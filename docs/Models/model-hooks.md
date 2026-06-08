# Model Hooks

SELF models expose a small set of **overridable hook methods** that let you
inject custom work into the time-integration loop without modifying the
integrators or the core tendency calculation. Every hook has a default
no-op implementation on the base `Model` type, so overriding them is
entirely optional and existing models are unaffected.

You override a hook the same way you override any other type-bound
procedure: by extending a concrete model type and binding your own
implementation.

## Available Hooks

| Hook | When it fires | Cadence | Typical use |
|------|---------------|---------|-------------|
| `PreStepHook`     | Immediately **before** a time step is taken (before any Runge-Kutta stages) | Once per **time step** | Per-step setup; refresh quantities that are constant across the step |
| `PreTendencyHook` | At the **start of every tendency evaluation** (`CalculateTendency`) | Once per **RK stage** (e.g. 3× per step for RK3) | Recompute derived fields the flux/source needs before each stage |
| `PostStepHook`    | Immediately **after** a completed time step (after `t` advances by `dt`) | Once per **time step** | Record receiver samples, accumulate diagnostics, on-the-fly output |

!!! note "Step vs. stage cadence"
    `PreStepHook` and `PostStepHook` fire **once per time step**, while
    `PreTendencyHook` fires **once per Runge-Kutta stage** — three times per
    step for the default `rk3` integrator, five for `rk4`. Choose the hook
    whose cadence matches your intent: use a step hook for work that should
    happen once per step, and `PreTendencyHook` only for work the tendency
    calculation depends on at every stage.

## Where the hooks fire

All four time integrators (`euler`, `rk2`, `rk3`, `rk4`) follow the same
pattern. Schematically, a single step of the low-storage RK3 integrator is:

```fortran
do while(this%t < tn)
  t0 = this%t
  call this%PreStepHook()              ! <-- once per step (start)
  do m = 1, 3
    call this%CalculateTendency()      !     PreTendencyHook() fires inside,
    call this%UpdateGRK3(m)            !     at the start of CalculateTendency
    this%t = t0 + rk3_b(m)*this%dt
  enddo
  this%t = t0 + this%dt
  call this%PostStepHook()             ! <-- once per step (end)
enddo
```

`PreTendencyHook` is invoked from within each model's `CalculateTendency`
(before the boundary exchange / flux / divergence work), so it runs once for
every stage `m`.

## Overriding a hook

Bind your implementation in the `contains` block of your extended type:

```fortran
module my_model_module

  use self_lineareuler2d

  implicit none

  type, extends(LinearEuler2D) :: my_model
    ! ... your model data (e.g. receiver buffers) ...
  contains
    procedure :: PostStepHook => PostStepHook_my_model
  endtype my_model

contains

  subroutine PostStepHook_my_model(this)
    implicit none
    class(my_model), intent(inout) :: this
    ! ... per-step work, e.g. sample the device-resident solution ...
  endsubroutine PostStepHook_my_model

endmodule my_model_module
```

The hook receives the model instance (`class(Model)`-compatible) and is
called automatically by the integrator — you do not call it yourself.

## Why hooks instead of stepping manually

A common reason to reach for a hook is **per-step receiver sampling or
output**. It is tempting to drive the integrator one step at a time from the
host:

```fortran
do k = 1, nsteps
  call modelobj%ForwardStep(t + dt, dt, dt)   ! NOT recommended for this
  ! ... sample here ...
enddo
```

but `ForwardStep` performs the entropy diagnostic, metric reporting, and
`WriteModel` file I/O on every interval — including device-to-host copies of
the full solution — which is wasteful when you only want a few sampled
values. Doing the work in `PostStepHook` keeps it **inside the native
time-stepping loop**, so the bulk solution stays resident on the device and
you copy back only what you need.

```fortran
subroutine PostStepHook_my_model(this)
  implicit none
  class(my_model), intent(inout) :: this
  ! Interpolate the device-resident solution at receiver points (on device)
  ! and copy back only the sampled values, not the whole field.
  call this%receivers%EvaluateScalar(this%solution, this%rxvals_dev)
  ! ... copy rxvals_dev -> host, store into a trace buffer ...
endsubroutine PostStepHook_my_model
```

## Guidelines

- **Keep hooks lightweight.** They run inside the time loop; expensive work
  (especially host/device transfers) directly inflates runtime per step.
- **Prefer device-resident operations.** Operate on the `*_gpu` buffers and
  copy back only small, sampled quantities.
- **Match the cadence.** Use `PreStepHook`/`PostStepHook` for once-per-step
  work and `PreTendencyHook` for work the tendency depends on every stage.
- **Hooks are optional.** The base-class defaults are no-ops; override only
  the hooks you need.
