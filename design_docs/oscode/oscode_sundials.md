# Proposal For Including OSCODE Into Sundials

## Summary

This design document describes the motivation and architecture for adding the OSCODE numerical routine for the efficient solution of second order, linear, homogeneous Ordinary Differential Equations (ODEs) with rapidly oscillating solutions. The proposed solver will will operate with initial value problems (IVPs) of the form

$$ 
u(t_0) = u_0, \quad u(t_1) = u_1,
$$
$$
\ddot{u}(t) + 2 \gamma(t) \dot{u}(t) + \omega^2(t)u(t) = 0, \quad t \in [t_0, t_1], 
$$

where $\gamma$ and $\omega$ may or may not be expressed as a closed form function of time. We further require $\gamma$ and $\omega$ to be real-valued and $\omega > 0$ If $\omega$ is large, the components may become highly oscillatory, and standard (polynomial-based) numerical methods will require $\mathcal{O}(\omega)$ timesteps/discretization points. In regimes where $\omega$ is large and smooth, OSCODE exploits an asympotic approximation to reduce this computational cost to $\mathcal{O}(1)$ ($\omega$-independent). In other regimes of the solution interval, OSCODE behaves as a Runge--Kutta solver, and is thus robust to changes in the behavior in the solution from oscillatory to non-oscillatory.

## Motivation

Solutions to oscillatory problems are often bespoke and use analytic approximations such as the Wentzel-Kramers-Brillouin approximation (WKB) method used by the current OSCODE solver. Giving researchers access to a general purpose solver for solutions for ODEs of the form written above will allow them to explore more computationally complex models in their research and speed up existing implementations.

The second order ODE can be reposed as a set of two, first order ODEs to fit into the standard schema for sundials solvers,

$$
\dot{\mathbf{y}} = \begin{pmatrix}
\mathbf{y}[1]\\
-\mathbf{y}[0]\omega^2(t)-2\gamma(t)\mathbf{y}[1]
\end{pmatrix}
$$

Among the SUNDIALS modules, ARKODE --which solves equations with known explicit and implicit, stiff and nonstill timescale components-- is structured most similarly to OSCODE.

$$
M(t)\dot{\mathbf{y}} = f^E(t, \mathbf{y}) + f^I(t, \mathbf{y}), \quad \mathbf{y}(t_0) = \mathbf{y}_0
$$

For the new OSCODE module we propose the right-hand-side will be a subset of what ARKODE can handle, specifically,

$$
M(t)\dot{\mathbf{y}} = f^O(t, \mathbf{y}), \quad \mathbf{y}|(t_0) = \mathbf{y}_0, \\
$$
where $f^O(t)$ takes the specific form
$$
f^O(t, y) = \begin{pmatrix}
\mathbf{y}[1]\\
-\mathbf{y}[0]\omega^2(t)-2\gamma(t)\mathbf{y}[1]
\end{pmatrix}.
$$

As the solution $\mathbf{y}(t)$ may vary between oscillatory and smoothly varying, OSCODE will dynamically switch between a Runge-Kutta solver and its asymptotic Adaptive Riccati Defect Correction method (ARDC). At each step the solver will choose the method which yields the larger step length (while keeping the local error within user defined tolerance).


## Guide-level Explanation

------------------------------------
Section Goals:

Explain the proposal as if it was already included in the project and you were teaching it to another Sundials programmer in the manual. That generally means:

- Introducing new named concepts.
- Explaining the feature largely in terms of examples.
- Explaining how Sundials programmers should *think* about the feature, and how it should impact the way they use the relevant package. It should explain the impact as concretely as possible.
- If applicable, provide sample error messages, deprecation warnings, or migration guidance.

------------------------------------

(Note: Much of this is taken from the ARKODE docs and modified for OSCODE)

The OSCODE infrastructure provides adaptive-step time integration modules for a class of nonstiff Ordinary Differential Equations (ODEs) where the components may be highly oscillatory. Current users of ARKODE will notice a similar structure in the C API of OSCODE.

OSCODE supports ODE systems posed in the form

$$
M(t)\dot{\mathbf{y}} = f^O(t, \mathbf{y}), \quad \mathbf{y}|(t_0) = \mathbf{y}_0, \\
$$
where $f^O(t)$ takes the specific form
$$
f^O(t, y) = \begin{pmatrix}
\mathbf{y}[1]\\
-\mathbf{y}[0]\omega^2(t)-2\gamma(t)\mathbf{y}[1]
\end{pmatrix}.
$$

Here, $t$ is the independent variable, $\mathbf{y}$ is a set of two (complex-valued) dependent variables with one being the derivative of the other, $\omega$ and $\gamma$ are user specified, real-valued callables. The solver has the ability to switch dynamically between using a Runge-Kutta method and the adaptive Riccati defect correction method to advance the solution, while adaptively updating its stepsize.

### Adaptive single-step methods

The OSCODE framework is designed to support single step, IVP integration methods. The stepsize is determined by the time-stepping method (based on user-provided accuracy requirements). However, users may place minimum/maximum bounds on if desired.

At each step, either the ARDC or RK solver is chosen for that step based on their estimated local error. For RK, the error is generated by the difference in the $N$th and $N-1$th order RK steps. The ARDC step's error is calculated as the residual of a transformed form of the original ODE, by substituting in ARDC's solution estimate.

OSCODE's time stepping modules may be run in a variety of "modes":

- NORMAL: The solver will take internal steps until it has just overtaken a user-specified end time, $t_1$, in the direction of integration, i.e. $t_{n-1} < t_1 \le t_n$ for forward integration, or $t_{n} < t_1 \le t_{n-1}$ for backward integration. It will then compute an approximation to the solution $\mathbf{y}(t_1)$ by interpolation (using one of the dense output routines described in the section on interpolation below).

- ONE-STEP: The solver will only take a single internal step $\mathbf{y}_{n-1} \rightarrow \mathbf{y}_n$ and then return control back to the calling program. If this step will overtake $t_1$ then the solver will again return an interpolated result; otherwise it will return a copy of the internal solution $\mathbf{y}_n$.

- NORMAL-TSTOP: The solver will take internal steps until the point where the _next_ step would overtake $t_1$. It will then limit this next step so that $t_n = t_{n-1} + h_n = t_1$, and once the step completes it will return a copy of the internal solution.

- ONE-STEP-TSTOP: The solver will check whether the next step would overtake $t_1$. If not, then this mode is identical to "one-step" above; otherwise it will limit this next step so that $t_n = t_{n-1} + h_n = t_1$. In either case, once the step completes it will return a copy of the internal solution $\mathbf{y}_n$.

### Interpolation

As mentioned above, the time-stepping modules in OSCODE support interpolation of solution $\mathbf{y}(t_{\mathrm{out}})$ and derivatives $\frac{\mathrm{d}\mathbf{y}}{\mathrm{d}t}|_{t = t_{\mathrm{out}}}$, where $t_{\mathrm{out}}$ occurs within a completed time step from $[t_{n-1},  t_n]$. To this end, OSCODE currently supports construction of polynomial interpolants $p_q(t)$ where the degree $q$ is determined by a parameter passed to the solver by the user.

OSCODE provides two complementary interpolation approaches, both of which are accessible from any of the time-stepping modules: "Hermite" and "Lagrange". The Hermite approach has been included with ARKODE since its inception, and is more suitable for non-stiff problems; Lagrange is designed to provide increased accuracy when integrating stiff problems. These two methods are described in detail within the ARKODE documentation (link).

### User callable functions

```c++
/**
 * @param fe the name of the C function (of type `OSRhsFn()`)
 *   defining the explicit portion of the right-hand side function in
 *   `M(t)\, y'(t) = f^E(t,y) + f^I(t,y)`.
 * @param fi the name of the C function (of type `OSRhsFn()`)
 *       defining the implicit portion of the right-hand side function in
 *       `M(t)\, y'(t) = f^E(t,y) + f^I(t,y)`.
 * @param t0 the initial value of :math:`t`.
 * @param y0 the initial condition vector :math:`y(t_0)`.
 * @param sunctx  the `SUNContext` object (see :numref:`SUNDIALS.SUNContext`)
 * @return If successful, a pointer to initialized problem memory
 * of type ``void*``, to be passed to all user-facing OSCStep routines
 * listed below.  If unsuccessful, a ``NULL`` pointer will be
 * returned, and an error message will be printed to ``stderr``.
 */
void* OSCStepCreate(OSCRhsFn, f, realtype t0, N_Vector y0, SunContext sunctx);

void OSCStepFree(void** oscode_mem);

// Tolerances
int OSCStepSStolerances(void* oscode_mem, realtype reltol, realtype abstol);

int OSCStepSVtolerances(void* oscode_mem, realtype reltol, N_Vector abstol);

int OSCStepWFtolerances(void* oscode_mem, OSEwtFn efun);

int OSCStepResStolerance(void* oscode_mem, realtype rabstol);

int OSCStepResVtolerance(void* oscode_mem, N_Vector rabstol);

int OSCStepResFtolerance(void* oscode_mem, OSRwtFn rfun);

// Main function
int OSCStepEvolve(void* oscode_mem, realtype tout, N_Vector yout, realtype *tret, int itask);

// Set Interpolator
int OSCStepSetInterpolantType(void* oscode_mem, int itype);

int OSCStepSetInterpolantDegree(void* oscode_mem, int degree);

// Error and diagnostic files
int OSCStepSetDiagnostics(void* oscode_mem, FILE* diagfp);

int OSCStepSetErrFile(void* oscode_mem, FILE* errfp);

int OSCStepSetErrHandlerFn(void* oscode_mem, OSErrHandlerFn ehfun, void* eh_data);

// initial size, and stopping rules
int OSCStepSetInitStep(void* oscode_mem, realtype hin);

int OSCStepSetMaxHnilWarns(void* oscode_mem, int mxhnil);

int OSCStepSetMaxNumSteps(void* oscode_mem, long int mxsteps);

int OSCStepSetMaxStep(void* oscode_mem, realtype hmax);

int OSCStepSetMinStep(void* oscode_mem, realtype hmin);

int OSCStepSetMinStep(void* oscode_mem, realtype hmin);

int OSCStepSetInterpolateStopTime(void* oscode_mem, booleantype interp);

int OSCStepClearStopTime(void* oscode_mem);

// User data passed into functions
int OSCStepSetUserData(void* oscode_mem, void* user_data);

// ARDC solver
int OSCStepARDCSetOrder(void* oscode_mem, int ord);

// Use only RK or ARDC (Mostly for debugging purposes probably)
// OSCODE_MODE is an enum of {RK, ARDC, BOTH} with default of BOTH
int OSCStepSetExplicit(void* oscode_mem, OSCODE_MODE mode);

// Set order and tables for ARDC solver's Chebyshev gradient method
int OSCStepARDCSetDerivativeMatrices(void* oscode_mem, int o, OSodeButcherTable Bi);

// If the user sets the ARDC tables and order they must also set the series and derivative functions.
int OSCStepARDCSetSeriesFuncs(void* oscode_mem, OSCddsFn dds, OSCdsiFn dsi, OSCdsfFn dsf);

// Set order and tables for RK solver 
int OSCStepRKSetTables(void* oscode_mem, int o, OSodeButcherTable Bi)

int OSCStepRKSetTableNum(void* oscode_mem, OSODE_ERKTableID etable)

typedef double (*fun_ptr)(double) OSCGradFun;
struct OSCGradFuns {
    OSCGradFun* first_deriv_;
    OSCGradFun* second_deriv_;
};
enum OSCODE_GRADIENT_MODE {Chebyshev, UserDefined};
// If user sets mode to `UserDefined` g must have pointers to functions returning nth gradient
OSCStepARDCGradientCalculation(OSCODE_GRADIENT_MODE mode, OSCGradFuns g);
```

Should we have methods for setting the time step adaptivity? Ex


| Optional input                                               | Function name                            | Default  |
|--------------------------------------------------------------|------------------------------------------|----------|
| Set a custom time step adaptivity function                   | `OSCStepSetAdaptivityFn()`      | internal |
| Choose an existing time step adaptivity method               | `OSCStepSetAdaptivityMethod()`  | 0        |
| Explicit stability safety factor                             | `OSCStepSetCFLFraction()`       | 0.5      |
| Time step error bias factor                                  | `OSCStepSetErrorBias()`         | 1.5      |
| Bounds determining no change in step size                    | `OSCStepSetFixedStepBounds()`   | 1.0 - 1.5|
| Maximum step growth factor on convergence fail               | `OSCStepSetMaxCFailGrowth()`    | 0.25     |
| Maximum step growth factor on error test fail                | `OSCStepSetMaxEFailGrowth()`    | 0.3      |
| Maximum first step growth factor                             | `OSCStepSetMaxFirstGrowth()`    | 10000.0  |
| Maximum allowed general step growth factor                   | `OSCStepSetMaxGrowth()`         | 20.0     |
| Minimum allowed step reduction factor on error test fail     | `OSCStepSetMinReduction()`      | 0.1      |
| Time step safety factor                                      | `OSCStepSetSafetyFactor()`      | 0.96     |
| Error fails before MaxEFailGrowth takes effect               | `OSCStepSetSmallNumEFails()`    | 2        |
| Explicit stability function                                  | `OSCStepSetStabilityFn()`       | none     |


## Reference-level explanation
------------------------------------
This is the technical portion of the RFC. Explain the design in sufficient detail that:

- Its interaction with other features is clear.
- It is reasonably clear how the feature would be implemented.
- Corner cases are dissected by example.

The section should return to the examples given in the previous section, and explain more fully how the detailed proposal makes those examples work.

------------------------------------

The main driver of `OSCODE` will be the function `oscEvolve` with a similar signature to `arkEvolve`

```C
/*---------------------------------------------------------------
  oscEvolve:

  This routine is the main driver of OSCODE-based integrators.

  It integrates over a time interval defined by the user, by
  calling the time step module to do internal time steps.

  The first time that oscEvolve is called for a successfully
  initialized problem, it computes a tentative initial step size.

  oscEvolve supports two modes as specified by itask: OSC_NORMAL and
  OSC_ONE_STEP.  In the OSC_NORMAL mode, the solver steps until
  it reaches or passes tout and then interpolates to obtain
  y(tout).  In the OSC_ONE_STEP mode, it takes one internal step
  and returns.  The behavior of both modes can be over-rided
  through user-specification of osc_tstop (through the
  *StepSetStopTime function), in which case if a solver step
  would pass tstop, the step is shortened so that it stops at
  exactly the specified stop time, and hence interpolation of
  y(tout) is not required.
  ---------------------------------------------------------------*/
int oscEvolve(OSCodeMem osc_mem, realtype tout, N_Vector yout,
              realtype *tret, int itask);
```

The `OSCodeMemRec` will be similar to `ARKodeMemRec` containing the information for the estimation using the rk and wkb solver

(TODO:) The example program `example_oscode.cpp` shows an example of the user API for the burst equation.

Users can set:

1. `w` and `g` functions that take in a timestep `t` along with a `void*` argument for their own data
2. `OSCODEOmegaFn` and `OSCODEGammaFn`, the functors $\omega(t)$ and $\gamma(t)$ which have pointer type
    - `complextype (OSCODEOmegaFn*)(realtype)`
3. Options for the rk solver
    - order
    - tolerances
    - butcher tables
4. Options for the ARDC solver
    - Order
    - tolerance(s)
    - User supplied gradient functions (default Chebyshev methods)
5. Options for solver
    - initial conditions for the ODE, $ x(t)  \frac{dx}{dt} $ evaluated at the start of the integration range
    - start of integration range
    - end of integration range
    - do_times timepoints at which dense output is to be produced. Must be sorted
    - Direction of integration
    - (local) relative tolerance
    - (local) absolute tolerance
    - initial stepsize suggestion (to be refined by the solver)

The implimentation can utilize the ARKODE explicit runga kutta solver, but we will have to write an implementation of the ARDC solver. The ARDC solver is currently only written in python, but will be ported to C for use in OSCODE.

Tests and example code will be written for the airy and burst examples.

## Drawbacks

The main question is the user base. As long as we think there enough researchers that would use this method then it most likely makes sense to add.

The other main question to this proposal is whether it's possible to make this fit in the ARKODE framework.

## Rationale and alternatives

------------------------------------

- Why is this design the best in the space of possible designs?
- What other designs have been considered and what is the rationale for not choosing them?
- What is the impact of not doing this?

------------------------------------
## Prior art

Discuss prior art, both the good and the bad, in relation to this proposal.
A few examples of what this can include are:

- For language, library, tools, and compiler proposals: Does this feature exist in other programming languages and what experience have their community had?

The OSCODE solver is currently available as a python package and standalone header files. the riccati solver is also available as it's own python package.

## Papers

Are there any published papers or great posts that discuss this? If you have some relevant papers to refer to, this can serve as a more detailed theoretical background.

- Original paper: ["Efficient method for solving highly oscillatory ordinary differential equations
with applications to physical systems"](https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.2.013030)
- JOSS Paper on Python / C++ implementation of OSCODE: [Link](https://joss.theoj.org/papers/10.21105/joss.02830)
- Riccati paper: ["An adaptive spectral method for oscillatory second-order linear ODEs with frequency-independent cost"](https://arxiv.org/abs/2212.06924)
- C++ Implementation: [Link](https://github.com/fruzsinaagocs/oscode)
- Python Wrapper Docs: [Link](https://oscode.readthedocs.io/en/latest/)
- Riccati python package [Link](https://github.com/fruzsinaagocs/riccati)

## Unresolved Questions

- What parts of the design do you expect to resolve through the RFC process before this gets merged?

Whether OSCODE should be its own seperate module or whether it is possible to fit it in the ARKODE framework.

- What parts of the design do you expect to resolve through the implementation of this feature before stabilization?

- What related issues do you consider out of scope for this RFC that could be addressed in the future independently of the solution that comes out of this RFC?

Any reworking of the algorithm itself would be out of scope of this document.
