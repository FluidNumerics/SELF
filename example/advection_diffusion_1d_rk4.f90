program advection_diffusion_1d_rk4

  use self_data
  use self_advection_diffusion_1d

  implicit none
  character(SELF_INTEGRATOR_LENGTH), parameter :: integrator = 'rk4'
  integer, parameter :: nvar = 1
  integer, parameter :: nelem = 50
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  real(prec), parameter :: u = 1.0_prec ! velocity
  real(prec), parameter :: nu = 0.001_prec ! diffusivity
  real(prec), parameter :: dt = 1.0_prec*10.0_prec**(-4) ! time-step size
  real(prec), parameter :: endtime = 1.0_prec
  real(prec), parameter :: iointerval = 0.1_prec
  type(advection_diffusion_1d) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh1D),target :: mesh
  type(Geometry1D),target :: geometry
  type(MPILayer),target :: decomp

  ! Create a mesh using the built-in
  ! uniform mesh generator.
  ! The domain is set to x in [0,1]
  ! We use `nelem` elements
  call mesh % UniformBlockMesh(nGeo=1,&
                               nElem=nelem,&
                               x=(/0.0_prec,1.0_prec/))

  ! We create a domain decomposition.
  call decomp % Init(enableMPI=.false.)
  call decomp % GenerateDecomposition(nelem,nelem+1)

  ! Create an interpolant
  call interp % Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

  ! Generate geometry (metric terms) from the mesh elements
  call geometry % Init(interp,mesh % nElem)
  call geometry % GenerateFromMesh(mesh)

  ! Initialize the model
  call modelobj % Init(nvar,mesh,geometry,decomp)

  ! Set the velocity
  modelobj % u = u
  !Set the diffusivity
  modelobj % nu = nu

  ! Set the initial condition
  call modelobj % solution % SetEquation( 1, 'f = 1.0')
  call modelobj % solution % SetInteriorFromEquation( geometry, 0.0_prec ) 

  print*, "min, max (interior)", &
    minval(modelobj % solution % interior % hostdata), &
    maxval(modelobj % solution % interior % hostdata)

  ! Set the model's time integration method
  call modelobj % SetTimeIntegrator( integrator )

  ! forward step the model to `endtime` using a time step
  ! of `dt` and outputing model data every `iointerval`
  call modelobj % ForwardStep(endtime,dt,iointerval)

  print*, "min, max (interior)", &
  minval(modelobj % solution % interior % hostdata), &
  maxval(modelobj % solution % interior % hostdata)
  

  ! Clean up
  call modelobj % free()
  call decomp % free()
  call mesh % free()
  call geometry % free()
  call interp % free()
  
end program advection_diffusion_1d_rk4
