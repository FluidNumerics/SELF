# Viscous Burgers Equation

The [`SELF_Burgers1D_t` module](../ford/module/self_burgers1d_t.html) defines the [Burgers1D_t class](ford/type/burgers1d_t.html) class. In SELF, models are posed in the form of a generic conservation law

\begin{equation}
\vec{s}_t + \nabla \cdot \overleftrightarrow{f} = \vec{q}
\end{equation}

where $\vec{s}$ is a vector of solution variables, $\overleftrightarrow{f}$ is the conservative flux, and $\vec{q}$ are non-conservative source terms. 

For Burgers equation in 1-D

\begin{equation}
\vec{s} = s
\end{equation}

\begin{equation}
\overleftrightarrow{f} = \frac{s^2}{2} \hat{x}
\end{equation}

\begin{equation}
\vec{q} = 0
\end{equation}

To track stability of the Burgers equation in 1-D, the total entropy function is

\begin{equation}
e = \int_x \frac{s^2}{2} \hspace{1mm} dx
\end{equation}

## Implementation
The viscous Burgers equation model is implemented as a type extension of the [`DGModel1D` class](../ford/type/dgmodel1d_t.html). The [Burgers1D_t class](ford/type/burgers1d_t.html) adds a parameter for the viscosity and overrides the `SetMetadata`, `entropy_func`, `flux1d`, and `riemannflux1d` type-bound procedures.

## Riemann Solver
The `Burgers1D` class is defined using the conservative form of the conservation law. The Riemman solver for the hyperbolic part of Burgers equation is the local Lax Friedrichs upwind riemann solver

\begin{equation}
f_h^* = \frac{1}{2}( f_L + f_R + c_{max}(s_L - s_R))
\end{equation}

where 
\begin{equation}
f_L = \frac{1}{2}s_L^2
\end{equation}

\begin{equation}
f_R = \frac{1}{2}s_R^2
\end{equation}

and 
\begin{equation}
c_{max} = max( |s_L|, |s_R| )
\end{equation}

The parabolic part of the flux (the viscous flux) is computed using the Bassi-Rebay flux, which computes the flux using an average of the gradient values on either side of the shared edge.

\begin{equation}
f_p^* = -\frac{\nu}{2}\left( \frac{∂ s_L}{∂ x}+ \frac{∂ s_R}{∂ x}\right)
\end{equation}

The details for this implementation can be found in [`self_burgers1d_t.f90`](../ford/sourcefile/self_burgers1d_t.f90.html)

### Boundary conditions
By default, the boundary conditions are periodic boundary conditions. When initializing the mesh for your Burgers equation solver, you can change the boundary conditions to `SELF_BC_Radiation` to set the external state on model boundaries to 0 in the Riemann solver

```fortran

type(Mesh1D),target :: mesh

  ! Create a mesh using the built-in
  ! uniform mesh generator.
  ! The domain is set to x in [0,1] with 10 elements
  call mesh%UniformBlockMesh(nGeo=1, &
                             nElem=10, &
                             x=(/0.0_prec,1.0_prec/))

  ! Reset the boundary conditions to radiation
  call mesh%ResetBoundaryConditionType(SELF_BC_RADIATION,SELF_BC_RADIATION)

```

If you need to explicitly set the boundary conditions as a function of position and time, you can create a type-extension of the `Burgers1D` class and override the [`hbc1d_Prescribed`](../ford/proc/hbc1d_prescribed_model.html) and [`pbc1d_Prescribed`](../ford/proc/pbc1d_prescribed_model.html) boundary condition procedures.

To make a type extension, you can first create a module that defines your model with the the new type-bound procedures for the boundary conditions.

```fortran
module my_burgers_model

use self_burgers1d

implicit none

  type,extends(Burgers1D) :: myModel
  contains
    procedure :: hbc1d_Prescribed => hbc1d_mymodel ! For the hyperbolic part
    procedure :: pbc1d_Prescribed => pbc1d_mymodel ! For the parabolic part
  end type myModel

  contains
  pure function hbc1d_mymodel(this,x,t) result(exts)
    !! Sets the external solution state at model boundaries
    class(myModel),intent(in) :: this
    real(prec),intent(in) :: x
    real(prec),intent(in) :: t
    real(prec) :: exts(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      exts(ivar) = ! To do : fill in the external state 
                   !          here as a function of space and time
    enddo

  endfunction hbc1d_mymodel

  pure function pbc1d_mymodel(this,x,t) result(extDsdx)
    !! Sets the external solution state derivative at model boundaries
    class(myModel),intent(in) :: this
    real(prec),intent(in) :: x
    real(prec),intent(in) :: t
    real(prec) :: extDsdx(1:this%nvar)
    ! Local
    integer :: ivar

    do ivar = 1,this%nvar
      extDsdx(ivar) = ! To do : fill in the external state 
                      !          here as a function of space and time
    enddo

  endfunction pbc1d_mymodel
  
end module my_burgers_model
```

In your program, you can use your new class. Your new class will inherit all of the features and other type bound procedures from the `Burgers1D` class but will enforce your boundary conditions. The snippet below shows the steps required to instantiate your model

```fortran
program run_my_model

use my_burgers_model

implicit none

  type(mymodel) :: modelobj
  type(Lagrange),target :: interp
  type(Mesh1D),target :: mesh
  type(Geometry1D),target :: geometry

  call mesh % UniformBlockMesh(nGeo=1, &
                               nElem=10, &
                               x=(/0.0_prec,1.0_prec/))

  ! Set the left and right boundary conditions to prescribed
  call mesh % ResetBoundaryConditionType(SELF_BC_PRESCRIBED,SELF_BC_PRESCRIBED)

  ! Create a 7th degree polynomial interpolant
  ! with Legendre-Gauss quadrature.
  ! The target grid for plotting is 12 uniformly spaced
  ! points within each element.
  call interp % Init(N=7, &
                     controlNodeType=GAUSS, &
                     M=12, &
                     targetNodeType=UNIFORM)

  ! Generate geometry (metric terms) from the mesh elements
  call geometry % Init(interp,mesh%nElem)
  call geometry % GenerateFromMesh(mesh)

  ! Initialize the model
  call modelobj % Init(nvar,mesh,geometry)
  ! Enable gradient calculations
  ! so that we can compute diffusive fluxes
  modelobj % gradient_enabled = .true.

  ! To do : Set the initial conditions
  ! To do : Forward step the model

end program run_my_model
```



## Example usage

For examples, see any of the following

* [`examples/burgers1d_shock.f90`](https://github.com/FluidNumerics/SELF/blob/main/examples/burgers1d_shock.f90) - Implements the travelling viscous shock problem