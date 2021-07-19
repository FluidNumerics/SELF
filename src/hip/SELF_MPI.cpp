
__global__ void ApplyFlip_MappedScalar2D_gpu(real *extBoundary, real *boundary, int *elemInfo, int *sideInfo, int *elemToRank, int rankId, int N, int nVar, int nEl){

  size_t ivar = blockIdx.x;
  size_t s1 = blockIdx.y;
  size_t e1 = blockIdx.z;
  size_t i1 = threadIdx.x;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;
  int i2 = N-i1;

  __shared__ real extBuff[16];

  extBuff[i1] = extBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0 && neighborRank /= rankId){
    if(flip == 1){
      extBoundary[SCB_2D_INDEX(i1,ivar,s1,e1,N,nVar)] = extBuff[i2];
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedScalar2D_gpu_wrapper(int **sideInfo, int **elemToRank, real **extBoundary, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedScalar2D_gpu<<<dim3(nVar,4,nEl), dim3(N+1,1,1), 0, 0>>>(*sideInfo, *elemToRank, *extBoundary, rankId, N, nVar);
  }

}

__global__ void ApplyFlip_MappedVector2D_gpu(real *extBoundary, real *boundary, int *elemInfo, int *sideInfo, int *elemToRank, int rankId, int N, int nVar, int nEl){

  size_t ivar = blockIdx.x;
  size_t s1 = blockIdx.y;
  size_t e1 = blockIdx.z;
  size_t dir = threadIdx.x;
  size_t i1 = threadIdx.y;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;
  int i2 = N-i1;

  __shared__ real extBuff[32];

  extBuff[dir+2*i1] = extBoundary[VEB_2D_INDEX(dir+1,i1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0 && neighborRank /= rankId){
    if(flip == 1){
      extBoundary[VEB_2D_INDEX(dir+1,i1,ivar,s1,e1,N,nVar)] = extBuff[dir+2*i2];
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedVector2D_gpu_wrapper(int **sideInfo, int **elemToRank, real **extBoundary, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedVector2D_gpu<<<dim3(nVar,4,nEl), dim3(2,N+1,1), 0, 0>>>(*sideInfo, *elemToRank, *extBoundary, rankId, N, nVar);
  }

}

__global__ void ApplyFlip_MappedTensor2D_gpu(real *extBoundary, real *boundary, int *elemInfo, int *sideInfo, int *elemToRank, int rankId, int N, int nVar, int nEl){

  size_t ivar = blockIdx.x;
  size_t s1 = blockIdx.y;
  size_t e1 = blockIdx.z;
  size_t row = threadIdx.x;
  size_t col = threadIdx.y;
  size_t i1 = threadIdx.z;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;
  int i2 = N-i1;

  __shared__ real extBuff[64];

  extBuff[row+2*(col+2*i1)] = extBoundary[TEB_2D_INDEX(row+1,col+1,i1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0 && neighborRank /= rankId){
    if(flip == 1){
      extBoundary[TEB_2D_INDEX(row+1,col+1,i1,ivar,s1,e1,N,nVar)] = extBuff[row+2*(col+2*i2)];
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedTensor2D_gpu_wrapper(int **sideInfo, int **elemToRank, real **extBoundary, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedTensor2D_gpu<<<dim3(nVar,4,nEl), dim3(2,2,N+1), 0, 0>>>(*sideInfo, *elemToRank, *extBoundary, rankId, N, nVar);
  }

}

__global__ void ApplyFlip_MappedScalar3D_gpu(real *extBoundary, real *boundary, int *elemInfo, int *sideInfo, int *elemToRank, int rankId, int N, int nVar, int nEl){

  size_t ivar = blockIdx.x;
  size_t s1 = blockIdx.y;
  size_t e1 = blockIdx.z;
  size_t i1 = threadIdx.x;
  size_t j1 = threadIdx.x;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;

  __shared__ real extBuff[256];

  extBuff[i1+(N+1)*j1] = extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0 && neighborRank /= rankId){
    if(flip == 2){
      int i2 = N-j1;
      int j2 = i1;
      extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = extBuff[i2+(N+1)*j2];
    }
    else if(flip == 3){
      int i2 = N-i1;
      int j2 = N-j1;
      extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = extBuff[i2+(N+1)*j2];
    }
    else if(flip == 4){
      int i2 = j1;
      int j2 = N-i1;
      extBoundary[SCB_3D_INDEX(i1,j1,ivar,s1,e1,N,nVar)] = extBuff[i2+(N+1)*j2];
    }
  }
}

extern "C"
{
  void ApplyFlip_MappedScalar3D_gpu_wrapper(int **sideInfo, int **elemToRank, real **extBoundary, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedScalar3D_gpu<<<dim3(nVar,6,nEl), dim3(N+1,N+1,1), 0, 0>>>(*sideInfo, *elemToRank, *extBoundary, rankId, N, nVar);
  }

}

__global__ void ApplyFlip_MappedVector3D_gpu(real *extBoundary, real *boundary, int *elemInfo, int *sideInfo, int *elemToRank, int rankId, int N, int nVar, int nEl){

  size_t ivar = blockIdx.x;
  size_t s1 = blockIdx.y;
  size_t e1 = blockIdx.z;
  size_t dir = threadIdx.x;
  size_t i1 = threadIdx.y;
  size_t j1 = threadIdx.z;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;

  __shared__ real extBuff[768];

  extBuff[dir+3*(i1+(N+1)*j1)] = extBoundary[VEB_3D_INDEX(dir+1,i1,j1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0 && neighborRank /= rankId){
    if(flip == 2){
      int i2 = N-j1;
      int j2 = i1;
      extBoundary[VEB_3D_INDEX(dir,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[dir+3*(i2+(N+1)*j2)];
    }
    else if(flip == 3){
      int i2 = N-i1;
      int j2 = N-j1;
      extBoundary[VEB_3D_INDEX(dir,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[dir+3*(i2+(N+1)*j2)];
    }
    else if(flip == 4){
      int i2 = j1;
      int j2 = N-i1;
      extBoundary[VEB_3D_INDEX(dir,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[dir+3*(i2+(N+1)*j2)];
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedVector3D_gpu_wrapper(int **sideInfo, int **elemToRank, real **extBoundary, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedVector3D_gpu<<<dim3(nVar,6,nEl), dim3(3,N+1,N+1), 0, 0>>>(*sideInfo, *elemToRank, *extBoundary, rankId, N, nVar);
  }

}

__global__ void ApplyFlip_MappedTensor3D_gpu(real *extBoundary, real *boundary, int *elemInfo, int *sideInfo, int *elemToRank, int rankId, int N, int nVar, int nEl){

  size_t ivar = blockIdx.x;
  size_t s1 = blockIdx.y;
  size_t e1 = blockIdx.z;
  size_t dir = threadIdx.x;
  size_t i1 = threadIdx.y;
  size_t j1 = threadIdx.z;
  size_t row = dir/3;
  size_t col = dir - dir*row;
  
  int e2 = sideInfo[INDEX3(2,s1,e1,5,4)];
  int s2 = sideInfo[INDEX3(3,s1,e1,5,4)]/10;
  int flip = sideInfo[INDEX3(3,s1,e1,5,4)]-s2*10;

  __shared__ real extBuff[2304];

  extBuff[row+3*(col+3*(i1+(N+1)*j1))] = extBoundary[TEB_3D_INDEX(row+1,col+1,i1,j1,ivar,s1,e1,N,nVar)];

  __syncthreads();

  if(bcid == 0 && neighborRank /= rankId){
    if(flip == 2){
      int i2 = N-j1;
      int j2 = i1;
      extBoundary[TEB_3D_INDEX(row+1,col+1,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[row+3*(col+3*(i1+(N+1)*j1))];
    }
    else if(flip == 3){
      int i2 = N-i1;
      int j2 = N-j1;
      extBoundary[TEB_3D_INDEX(row+1,col+1,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[row+3*(col+3*(i1+(N+1)*j1))];
    }
    else if(flip == 4){
      int i2 = j1;
      int j2 = N-i1;
      extBoundary[TEB_3D_INDEX(row+1,col+1,i1,j1,ivar,s1,e1,N,nVar)] = extBuff[row+3*(col+3*(i1+(N+1)*j1))];
    }
  }
  
}

extern "C"
{
  void ApplyFlip_MappedTensor3D_gpu_wrapper(int **sideInfo, int **elemToRank, real **extBoundary, int rankId, int N, int nVar, int nEl)
  {
    ApplyFlip_MappedTensor3D_gpu<<<dim3(nVar,6,nEl), dim3(9,N+1,N+1), 0, 0>>>(*sideInfo, *elemToRank, *extBoundary, rankId, N, nVar);
  }

}
