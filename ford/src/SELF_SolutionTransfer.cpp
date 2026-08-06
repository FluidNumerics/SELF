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
// the work per descent step. The arithmetic per retained value is identical in form but the
// operation order across a descent differs, so device and host results agree to round-off
// rather than bitwise; conservation of the Jacobian-weighted integral is exact either way and
// is what the tests assert.
//
// Thread mapping: one workgroup per new element, one thread per (i,j) node of that element,
// variables handled sequentially. The Np x Np working buffers therefore live in shared memory
// instead of per-thread scratch, which a thread-per-element mapping would spill to scratch
// memory at realistic degree and variable counts.

// Maximum supported (N+1). Sizes the shared working buffers below; the Fortran caller guards
// the degree so an unsupported N fails loudly rather than overrunning. Matches the bound used
// by the modal indicator in SELF_Refinement.cpp.
#define AMR2D_MAXNP 16

// Transfer-kind tags. Must match the SELF_TRANSFER_* parameters in SELF_TransferPlan_2D.f90.
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
                                             int pathStride, int eFirst0,
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
        int src = family[c + 4*gi] - 1;   // Fortran 1-based element index
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
      int src = sourceElem[gi] - 1;
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
                               int pathStride, int eFirst0,
                               int N, int nvar, int nOld, int nNew, int nLocal)
  {
    if( nLocal <= 0 ){ return; }
    int Np = N+1;
    TransferSolution_2D_gpukernel<<<dim3(nLocal,1,1), dim3(Np*Np,1,1), 0, 0>>>(
      uOld, uNew, sourceKind, sourceElem, family, depth, path,
      mortarR, mortarP, pathStride, eFirst0, N, nvar, nOld, nNew);
  }
}
