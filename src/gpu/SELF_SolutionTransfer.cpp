#include "SELF_GPU_Macros.h"

// Device-side AMR solution transfer (Stage 6a).
//
// Applies a prebuilt transfer plan on the device, replacing the host round trip that the
// adaptive loop previously performed every epoch (solution UpdateHost -> host
// ApplyTransferPlanRange -> solution UpdateDevice). The mathematics mirror the portable
// reference implementation exactly: SELF_TransferPlan_2D.f90 (ApplyTransferPlanRange) driving
// SELF_SolutionTransfer_2D.f90 (ProlongToChildren / RestrictFromChildren).
//
// Each new element is classified by the plan as
//
//   COPY     - the element is unchanged; take the old element's nodal values verbatim.
//   PROLONG  - the element descends from a surviving coarser ancestor; interpolate down.
//   RESTRICT - the element replaces a family of four children; L2-project them up.
//
// and may then require `depth` further prolongation steps down a recorded quadrant `path`.
// Both element-local operators are tensor products of the 1-D mortar matrices, applied as an
// xi pass followed by an eta pass:
//
//   prolong:  tmp(i,j)  = sum_ii mortarR(ii,i,kx) * parent(ii,j)
//             child(i,j)= sum_jj mortarR(jj,j,ky) * tmp(i,jj)
//   restrict: tmp(i,jj) = sum_ii mortarP(ii,i,kx) * child(ii,jj)
//             parent   += sum_jj mortarP(jj,j,ky) * tmp(i,jj)     (accumulated over children)
//
// mortarP carries the sub-edge Jacobian, which is what makes restriction conservative.
//
// Deliberate divergence from the host reference: the host descent calls ProlongToChildren,
// which forms all four children and discards three (SELF_TransferPlan_2D.f90:283-286). This
// kernel applies only the operator pair of the child actually on the path, doing a quarter of
// the work per descent step. What is dropped is the three discarded children - independent
// computations writing disjoint slices - and not any term of the retained child, whose
// contractions sum the same products against the same mortar column in the same ascending index
// order as the host. (An earlier revision of this comment claimed the operation order across a
// descent differs; it does not. See the 3-D header below for the same argument spelled out.)
// Device and host nonetheless agree to round-off rather than bitwise, because the device
// compiler contracts these multiply-accumulates into FMAs. Conservation of the Jacobian-weighted
// integral is exact either way and is what the tests assert.
//
// Thread mapping: one workgroup per new element, one thread per (i,j) node of that element,
// variables handled sequentially. The Np x Np working buffers therefore live in shared memory
// instead of per-thread scratch, which a thread-per-element mapping would spill to scratch
// memory at realistic degree and variable counts.

// The uOld contract, shared by both kernels: uOld holds a contiguous WINDOW of old elements,
// laid out (node, element, variable) with `nOld` as its per-variable element stride. `oldFirst0`
// is the 0-based global old index of the window's first element, so a plan's global source index
// is rebased as `src - oldFirst0`. The whole-field case is exactly oldFirst0 == 0 with nOld the
// global old count, which is what the single-rank path passes, so its indexing is unchanged.
//
// The window form is what lets a multi-rank adaptation run this kernel at all: the old elements a
// rank's new range reads are migrated to it point-to-point (SELF_SolutionMigration.f90) rather
// than allgathered, and on a GPU build they are received directly into device memory.
//
// The caller guarantees oldFirst0 <= src < oldFirst0 + nOld for every element it launches, by
// scanning the plan over [eFirst,eLast] before the launch (ApplyTransferPlan_DGModel{2,3}D) - the
// same three conditions ApplyTransferPlanWindow enforces on the host. That check is deliberately
// NOT repeated here: an out-of-range early return would not be block-uniform and would strand the
// __syncthreads() calls below.

// Maximum supported (N+1). Sizes the shared working buffers below; the Fortran caller guards
// the degree so an unsupported N fails loudly rather than overrunning (the same literal appears
// in ApplyTransferPlan_DGModel2D and its error message - keep the three in step). Matches the
// bound used by the modal indicator in SELF_Refinement.cpp.
#define AMR2D_MAXNP 16

// Transfer-kind tags, shared by the 2-D and 3-D kernels. Must match the SELF_TRANSFER_*
// parameters in BOTH SELF_TransferPlan_2D.f90 and SELF_TransferPlan_3D.f90, which agree.
#define SELF_TRANSFER_COPY 0
#define SELF_TRANSFER_PROLONG 1
#define SELF_TRANSFER_RESTRICT 2

// Child quadrant (1=SW,2=SE,3=NE,4=NW) -> (x-half, y-half); half h selects mortar sub-edge
// k = h. Mirrors transferAxc / transferAyc in SELF_SolutionTransfer_2D.f90, less the Fortran
// 1-based offset.
__device__ __constant__ int d_transferAxc[4] = {0, 1, 1, 0};
__device__ __constant__ int d_transferAyc[4] = {0, 0, 1, 1};

__global__ void TransferSolution_2D_gpukernel(real *uOld, real *uNew,
                                             int *sourceKind, int *sourceElem, int *family,
                                             int *depth, int *path,
                                             real *mortarR, real *mortarP,
                                             int pathStride, int eFirst0, int oldFirst0,
                                             int N, int nvar, int nOld, int nNew){

  // One block per new (rank-local) element.
  int lo = blockIdx.x;              // 0-based index into the rank-local new field
  int gi = eFirst0 + lo;            // 0-based index into the plan's global new-leaf arrays

  int Np = N+1;
  int tid = threadIdx.x;            // blockDim.x == Np*Np, so every thread is active
  int i = tid % Np;
  int j = tid / Np;

  __shared__ real buf[AMR2D_MAXNP*AMR2D_MAXNP];
  __shared__ real tmp[AMR2D_MAXNP*AMR2D_MAXNP];
  __shared__ real acc[AMR2D_MAXNP*AMR2D_MAXNP];

  int kind = sourceKind[gi];
  int d = depth[gi];

  for(int v=0; v<nvar; v++){

    if( kind == SELF_TRANSFER_RESTRICT ){

      // L2-project the four children onto their parent, accumulating over children.
      acc[i + Np*j] = 0.0;
      __syncthreads();

      for(int c=0; c<4; c++){
        // Fortran 1-based global old index, rebased onto uOld's window (oldFirst0 == 0 for a
        // whole-field uOld, so the single-rank case is unchanged arithmetic).
        int src = family[c + 4*gi] - 1 - oldFirst0;
        int kx = d_transferAxc[c];
        int ky = d_transferAyc[c];

        buf[i + Np*j] = uOld[SC_2D_INDEX(i,j,src,v,N,nOld)];
        __syncthreads();

        // xi pass: tmp(i,jj) = sum_ii mortarP(ii,i,kx) * child(ii,jj)
        real s = 0.0;
        for(int ii=0; ii<Np; ii++){
          s += mortarP[ii + Np*(i + Np*kx)]*buf[ii + Np*j];
        }
        tmp[i + Np*j] = s;
        __syncthreads();

        // eta pass, accumulated over the four children
        s = 0.0;
        for(int jj=0; jj<Np; jj++){
          s += mortarP[jj + Np*(j + Np*ky)]*tmp[i + Np*jj];
        }
        acc[i + Np*j] += s;
        __syncthreads();
      }

      buf[i + Np*j] = acc[i + Np*j];
      __syncthreads();

    } else {

      // COPY (d == 0) and PROLONG (d > 0) both start from a single surviving old element.
      int src = sourceElem[gi] - 1 - oldFirst0;
      buf[i + Np*j] = uOld[SC_2D_INDEX(i,j,src,v,N,nOld)];
      __syncthreads();

    }

    // Descend the recorded path, prolonging one level per step onto the child that is on the
    // path (rather than onto all four).
    for(int step=0; step<d; step++){
      // pathStride is the ALLOCATED leading dimension of plan%path, i.e.
      // max(forest%MaxLevel(),1) - not plan%maxDepth, which is the largest depth actually
      // encountered and is generally smaller.
      int c = path[step + pathStride*gi] - 1;  // quadrant 1..4 -> 0..3
      int kx = d_transferAxc[c];
      int ky = d_transferAyc[c];

      // xi pass: tmp(i,j) = sum_ii mortarR(ii,i,kx) * parent(ii,j)
      real s = 0.0;
      for(int ii=0; ii<Np; ii++){
        s += mortarR[ii + Np*(i + Np*kx)]*buf[ii + Np*j];
      }
      tmp[i + Np*j] = s;
      __syncthreads();

      // eta pass: child(i,j) = sum_jj mortarR(jj,j,ky) * tmp(i,jj)
      s = 0.0;
      for(int jj=0; jj<Np; jj++){
        s += mortarR[jj + Np*(j + Np*ky)]*tmp[i + Np*jj];
      }
      buf[i + Np*j] = s;
      __syncthreads();
    }

    uNew[SC_2D_INDEX(i,j,lo,v,N,nNew)] = buf[i + Np*j];
    __syncthreads();
  }

}

extern "C"
{
  void TransferSolution_2D_gpu(real *uOld, real *uNew,
                               int *sourceKind, int *sourceElem, int *family,
                               int *depth, int *path,
                               real *mortarR, real *mortarP,
                               int pathStride, int eFirst0, int oldFirst0,
                               int N, int nvar, int nOld, int nNew, int nLocal)
  {
    if( nLocal <= 0 ){ return; }
    int Np = N+1;
    TransferSolution_2D_gpukernel<<<dim3(nLocal,1,1), dim3(Np*Np,1,1), 0, 0>>>(
      uOld, uNew, sourceKind, sourceElem, family, depth, path,
      mortarR, mortarP, pathStride, eFirst0, oldFirst0, N, nvar, nOld, nNew);
  }
}

// ------------------------------------------------------------------------------------------- //
// 3-D device-side AMR solution transfer.
//
// The 3-D analogue of the kernel above, with the same contract: it applies a prebuilt transfer
// plan on the device, so the per-element interpolation runs on the GPU instead of on one CPU
// core and an adapting single-GPU run moves no solution data across the host link. Measured on
// one B300 at three refinement depths, that is worth ~10-15% off an adaptation, and essentially
// all of it is the interpolation rather than the link: at 20,784 elements the eliminated round
// trip is ~250 MB, some 5 ms against a 2.2 s adaptation. See section 6.4 of
// docs/Learning/AdaptiveMeshRefinement.md. The mathematics mirror the portable reference
// exactly - SELF_TransferPlan_3D.f90
// (ApplyTransferPlanRange) driving SELF_SolutionTransfer_3D.f90 (ProlongToChildren /
// RestrictFromChildren) - with the element-local operators now triple tensor products applied as
// an xi pass, an eta pass and a zeta pass:
//
//   prolong:  t1(i,j,k) = sum_ii mortarR(ii,i,kx) * parent(ii,j,k)
//             t2(i,j,k) = sum_jj mortarR(jj,j,ky) * t1(i,jj,k)
//             child     = sum_kk mortarR(kk,k,kz) * t2(i,j,kk)
//   restrict: t1(i,j,k) = sum_ii mortarP(ii,i,kx) * child(ii,j,k)
//             t2(i,j,k) = sum_jj mortarP(jj,j,ky) * t1(i,jj,k)
//             parent   += sum_kk mortarP(kk,k,kz) * t2(i,j,kk)     (accumulated over children)
//
// Each 1-D mortarP carries the half-interval Jacobian, so the triple product carries 1/8, which
// is what makes restriction conservative in 3-D.
//
// Thread mapping: one workgroup per new element, (N+1)*(N+1) threads, and thread (a,b) owns one
// (a,b) pencil and loops the remaining index in each of the three passes - the decomposition the
// 3-D modal indicator already uses (RefinementIndicator_3D_gpukernel in SELF_Refinement.cpp).
// An (N+1)^3 thread block is not an option: it is 4096 threads at N=15, past the block limit,
// and it would put the working buffers in per-thread scratch. Variables are handled sequentially.
//
// Deliberate divergence from the host reference, inherited from the 2-D kernel: the host descent
// calls ProlongToChildren, which forms all EIGHT children and discards seven
// (SELF_TransferPlan_3D.f90:288-291). This kernel applies only the operator triple of the child
// actually on the recorded path, doing an eighth of the work per descent step (a quarter in 2-D).
// The work that is dropped is the seven children that are thrown away, not any part of the
// retained value: each contraction below sums the same terms against the same mortar column, in
// the same ascending index order as the host loop, so the reduction order is preserved as
// CLAUDE.md requires.
//
// Device and host are nonetheless expected to agree only to round-off, not bitwise, because the
// device compiler contracts these multiply-accumulates into FMAs while the host build need not.
// That is the same situation as every other kernel in this directory. Conservation of the
// Jacobian-weighted integral is exact either way - it follows from sum_c P_kz P_ky P_kx being a
// partition of the parent's quadrature, not from the summation order - and that is what the AMR
// regression asserts per epoch; test/solution_transfer_3d_device.f90 pins the transferred values
// themselves against the host reference to the DOUBLE_PRECISION tolerance.

// Maximum supported (N+1) for the 3-D transfer kernel, and the size of the three (N+1)^3 working
// buffers, which live in per-block __shared__ memory. At 12 that is 3*12^3*8 = 41,472 B, inside
// both the 48 KB CUDA static shared-memory limit and the 64 KB gfx90a/gfx942 LDS budget. The
// Fortran caller guards the degree (ApplyTransferPlan_DGModel3D carries the same literal, in the
// test and in its error message - keep the three in step), so an unsupported N fails loudly
// rather than overrunning. Matches the bound the 3-D modal indicator uses (AMR3D_MAXNP in
// SELF_Refinement.cpp).
//
// Note the static allocation is AMR3D_MAXNP^3 whatever the degree in use, not (N+1)^3: 41,472 B
// on every launch, where N = 7 would only need 12,288 B. On an NVIDIA device that is not the
// limiter - a dynamic-shared variant sized 3*(N+1)^3*sizeof(real), which is what
// SELF_MatrixMultiply.cpp does, measured 1-3% SLOWER on a B300 across three runs at two
// refinement depths, so it was not adopted. It may still be worth revisiting on a 64 KB-LDS AMD
// device, where the static footprint admits only one workgroup per CU rather than five; that has
// not been measured, and since the alternative on every platform is the host round trip this
// kernel replaces, the static form cannot regress anything as it stands.
#define AMR3D_MAXNP 12

// Child octant (SELF/CGNS corner order) -> (x-half, y-half, z-half); half h selects mortar
// sub-interval k = h. Mirrors transferAxc / transferAyc / transferAzc in
// SELF_SolutionTransfer_3D.f90, less the Fortran 1-based offset.
__device__ __constant__ int d_transferAxc3D[8] = {0, 1, 1, 0, 0, 1, 1, 0};
__device__ __constant__ int d_transferAyc3D[8] = {0, 0, 1, 1, 0, 0, 1, 1};
__device__ __constant__ int d_transferAzc3D[8] = {0, 0, 0, 0, 1, 1, 1, 1};

__global__ void TransferSolution_3D_gpukernel(real *uOld, real *uNew,
                                              int *sourceKind, int *sourceElem, int *family,
                                              int *depth, int *path,
                                              real *mortarR, real *mortarP,
                                              int pathStride, int eFirst0, int oldFirst0,
                                              int N, int nvar, int nOld, int nNew){

  // One block per new (rank-local) element.
  int lo = blockIdx.x;              // 0-based index into the rank-local new field
  int gi = eFirst0 + lo;            // 0-based index into the plan's global new-leaf arrays

  int Np = N+1;
  int tid = threadIdx.x;            // blockDim.x == Np*Np, so every thread is active
  int a = tid % Np;                 // this thread's first pencil index
  int b = tid / Np;                 // this thread's second pencil index

  __shared__ real buf[AMR3D_MAXNP*AMR3D_MAXNP*AMR3D_MAXNP];
  __shared__ real tmp[AMR3D_MAXNP*AMR3D_MAXNP*AMR3D_MAXNP];
  __shared__ real acc[AMR3D_MAXNP*AMR3D_MAXNP*AMR3D_MAXNP];

  int kind = sourceKind[gi];
  int d = depth[gi];

  // Both branch conditions below are uniform across the block (they depend only on gi), so every
  // __syncthreads() is reached by every thread.
  for(int v=0; v<nvar; v++){

    if( kind == SELF_TRANSFER_RESTRICT ){

      // L2-project the eight children onto their parent, accumulating over children.
      for(int k=0; k<Np; k++){ acc[a + Np*(b + Np*k)] = 0.0; }
      __syncthreads();

      for(int c=0; c<8; c++){
        // Fortran 1-based global old index, rebased onto uOld's window (oldFirst0 == 0 for a
        // whole-field uOld, so the single-rank case is unchanged arithmetic).
        int src = family[c + 8*gi] - 1 - oldFirst0;
        int kx = d_transferAxc3D[c];
        int ky = d_transferAyc3D[c];
        int kz = d_transferAzc3D[c];

        // Stage the child into shared memory first. The xi pass contracts over the FIRST node
        // index, so reading uOld inside it would make every thread's address independent of a -
        // uncoalesced, and each child value re-read Np times over. Loading it here is coalesced
        // in a and reads each value once, as the 2-D kernel above does.
        for(int k=0; k<Np; k++){
          buf[a + Np*(b + Np*k)] = uOld[SC_3D_INDEX(a,b,k,src,v,N,nOld)];
        }
        __syncthreads();

        // xi pass: tmp(a,b,k) = sum_ii mortarP(ii,a,kx) * child(ii,b,k)
        for(int k=0; k<Np; k++){
          real s = 0.0;
          for(int ii=0; ii<Np; ii++){
            s += mortarP[ii + Np*(a + Np*kx)]*buf[ii + Np*(b + Np*k)];
          }
          tmp[a + Np*(b + Np*k)] = s;
        }
        // buf is overwritten by the eta pass below, and was just read across threads by the xi
        // pass, so this barrier closes both that read and the tmp write.
        __syncthreads();

        // eta pass: buf(a,b,k) = sum_jj mortarP(jj,b,ky) * tmp(a,jj,k)
        for(int k=0; k<Np; k++){
          real s = 0.0;
          for(int jj=0; jj<Np; jj++){
            s += mortarP[jj + Np*(b + Np*ky)]*tmp[a + Np*(jj + Np*k)];
          }
          buf[a + Np*(b + Np*k)] = s;
        }
        // tmp was just read across threads and the next child's xi pass rewrites it.
        __syncthreads();

        // zeta pass (own-pencil in both operands), accumulated over the eight children:
        //   acc(a,b,r) += sum_kk mortarP(kk,r,kz) * buf(a,b,kk)
        for(int r=0; r<Np; r++){
          real s = 0.0;
          for(int kk=0; kk<Np; kk++){
            s += mortarP[kk + Np*(r + Np*kz)]*buf[a + Np*(b + Np*kk)];
          }
          acc[a + Np*(b + Np*r)] += s;
        }
        // The next child's staging pass overwrites buf, which this child's zeta pass has just
        // read; that read is own-pencil, but the staging write is too, so this barrier is what
        // keeps the eight children ordered as a block rather than per thread.
        __syncthreads();
      }

      for(int k=0; k<Np; k++){ buf[a + Np*(b + Np*k)] = acc[a + Np*(b + Np*k)]; }
      __syncthreads();

    } else {

      // COPY (d == 0) and PROLONG (d > 0) both start from a single surviving old element.
      int src = sourceElem[gi] - 1 - oldFirst0;
      for(int k=0; k<Np; k++){
        buf[a + Np*(b + Np*k)] = uOld[SC_3D_INDEX(a,b,k,src,v,N,nOld)];
      }
      __syncthreads();

    }

    // Descend the recorded path, prolonging one level per step onto the child that is on the
    // path (rather than onto all eight).
    for(int step=0; step<d; step++){
      // pathStride is the ALLOCATED leading dimension of plan%path, i.e. max(forest%MaxLevel(),1)
      // - not plan%maxDepth, which is the largest depth actually encountered and is generally
      // smaller.
      int c = path[step + pathStride*gi] - 1;  // octant 1..8 -> 0..7
      int kx = d_transferAxc3D[c];
      int ky = d_transferAyc3D[c];
      int kz = d_transferAzc3D[c];

      // xi pass: tmp(a,b,k) = sum_ii mortarR(ii,a,kx) * parent(ii,b,k)
      for(int k=0; k<Np; k++){
        real s = 0.0;
        for(int ii=0; ii<Np; ii++){
          s += mortarR[ii + Np*(a + Np*kx)]*buf[ii + Np*(b + Np*k)];
        }
        tmp[a + Np*(b + Np*k)] = s;
      }
      __syncthreads();

      // eta pass: acc(a,b,k) = sum_jj mortarR(jj,b,ky) * tmp(a,jj,k)
      for(int k=0; k<Np; k++){
        real s = 0.0;
        for(int jj=0; jj<Np; jj++){
          s += mortarR[jj + Np*(b + Np*ky)]*tmp[a + Np*(jj + Np*k)];
        }
        acc[a + Np*(b + Np*k)] = s;
      }
      // acc is written and read within one thread's own pencil, so strictly this barrier is not
      // required by the zeta pass below; it is kept because the pass structure here is otherwise
      // identical to the cross-thread xi/eta pair above, and a future indexing change that made
      // the zeta contraction cross-thread would silently race without it.
      __syncthreads();

      // zeta pass: child(a,b,r) = sum_kk mortarR(kk,r,kz) * acc(a,b,kk)
      for(int r=0; r<Np; r++){
        real s = 0.0;
        for(int kk=0; kk<Np; kk++){
          s += mortarR[kk + Np*(r + Np*kz)]*acc[a + Np*(b + Np*kk)];
        }
        buf[a + Np*(b + Np*r)] = s;
      }
      __syncthreads();
    }

    for(int k=0; k<Np; k++){
      uNew[SC_3D_INDEX(a,b,k,lo,v,N,nNew)] = buf[a + Np*(b + Np*k)];
    }
    // The next variable's first pass overwrites buf (and, on the RESTRICT branch, acc), so the
    // stores above must complete for the whole block before it begins.
    __syncthreads();
  }

}

extern "C"
{
  void TransferSolution_3D_gpu(real *uOld, real *uNew,
                               int *sourceKind, int *sourceElem, int *family,
                               int *depth, int *path,
                               real *mortarR, real *mortarP,
                               int pathStride, int eFirst0, int oldFirst0,
                               int N, int nvar, int nOld, int nNew, int nLocal)
  {
    if( nLocal <= 0 ){ return; }
    int Np = N+1;
    TransferSolution_3D_gpukernel<<<dim3(nLocal,1,1), dim3(Np*Np,1,1), 0, 0>>>(
      uOld, uNew, sourceKind, sourceElem, family, depth, path,
      mortarR, mortarP, pathStride, eFirst0, oldFirst0, N, nvar, nOld, nNew);
  }
}
