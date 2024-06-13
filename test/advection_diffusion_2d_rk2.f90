program advection_diffusion_2d_rk2

  use self_data
  use self_advection_diffusion_2d

  implicit none
  character(SELF_INTEGRATOR_LENGTH),parameter :: integrator = 'rk2'
  integer,parameter :: nvar = 1
  integer,parameter :: controlDegree = 7
  integer,parameter :: targetDegree = 16
  real(prec),parameter :: u = 0.25_prec ! velocity
  real(prec),parameter :: v = 0.25_prec
  real(prec),parameter :: nu = 0.001_prec ! diffusivity
  real(prec),parameter :: dt = 1.0_prec*10.0_prec**(-4) ! time-step size
  real(prec),parameter :: endtime = 0.2_prec
  real(prec),parameter :: iointerval = 0.1_prec
  type(advection_diffusion_2d) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh2D),target :: mesh
  type(SEMQuad),target :: geometry
  type(MPILayer),target :: decomp
  character(LEN=255) :: WORKSPACE

  ! We create a domain decomposition.
  call decomp%Init(enableMPI=.false.)

  ! Create a uniform block mesh
  call get_environment_variable("WORKSPACE",WORKSPACE)
  call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D_mesh.h5",decomp)

  ! Create an interpolant
  call interp%Init(N=controlDegree, &
                   controlNodeType=GAUSS, &
                   M=targetDegree, &
                   targetNodeType=UNIFORM)

  ! Generate geometry (metric terms) from the mesh elements
  call geometry%Init(interp,mesh%nElem)
  call geometry%GenerateFromMesh(mesh)

  ! Initialize the model
  call modelobj%Init(nvar,mesh,geometry,decomp)

  ! Set the velocity
  modelobj%u = u
  modelobj%v = v
  !Set the diffusivity
  modelobj%nu = nu

  ! Set the initial condition
  call modelobj%solution%SetEquation(1,'f = \exp( -( (x-0.5)^2 + (y-0.5)^2 )/0.005 )')
  call modelobj%solution%SetInteriorFromEquation(geometry,0.0_prec)

  print*,"min, max (interior)", &
    minval(modelobj%solution%interior), &
    maxval(modelobj%solution%interior)

  ! Set the model's time integration method
  call modelobj%SetTimeIntegrator(integrator)

  ! forward step the model to `endtime` using a time step
  ! of `dt` and outputing model data every `iointerval`
  call modelobj%ForwardStep(endtime,dt,iointerval)

  print*,"min, max (interior)", &
    minval(modelobj%solution%interior), &
    maxval(modelobj%solution%interior)

  ! Clean up
  call modelobj%free()
  call decomp%free()
  call mesh%free()
  call geometry%free()
  call interp%free()

endprogram advection_diffusion_2d_rk2
