# Title

## Summary

This design document describes the motivation and architecture for adding the OSCODE numerical routine for the efficient solution of one dimensional, second order, linear, homogeneous ordinary differential equations with rapidly oscillating solutions. The proposed solver will will operate with differential equations of the form

$$ \ddot{x} + 2 \gamma(t) \dot{x}(t) + \omega^2(t)x(t) = 0 $$

where $\gamma$ and $\omega$ may or may not be expressed as a closed form function of time. We further require $\gamma$ and $\omega$ to be real-valued. If $\omega$ is large, the solution may be highly oscillatory, and standard (polynomial-based) numerical methods will require $\mathcal{O}(\omega)$ timesteps/discretization points. In regimes where $\omega$ is large and smooth, OSCODE exploits an asympotic approximation to reduce this computational cost to $\mathcal{O}(1)$ ($\omega$-independent). In other regimes of the solution interval, OSCODE behaves as a Runge--Kutta solver, and is thus robust to changes in the behavior in the solution from oscillatory to non-oscillatory.

<!--- We should cite the oscode papers for details, and the JOSS paper! --->

This new module will follow the same C based API with an underlying C++ implementation.

## Motivation

Solutions to oscillatory problems are often bespoke and use analytic approximations such as the Adaptive Ricatti Defect Correction (ARDC) method used by OSCODE. Giving researchers access to a general purpose solver for solutions for ODEs of the form written above will allow them to explore more computationally complex models in their research and speed up existing implementations.

The second order ODE can be reposed as a set of two, first order ODEs to fit into the standard schema for sundials solvers,

$$
\dot{y} = \begin{pmatrix}
y[1]\\
-y[0]w^2(t)-2\gamma(t)y[1]
\end{pmatrix}.
$$

The most similar module is ARKODE, which solves equations with known explicit and implicit, nonstiff and stiff time scale components,

$$
M(t)\ddot{y} = f^E(t, y) + f^I(t, y), \quad y(t_0) = y_0.
$$

<!--- Are you sure the \ddot above is correct? I'd expect just \dot.--->

For the new OSCODE module we propose the right-hand-side will be a subset of what ARKODE can handle, specifically,

$$
M(t)\ddot{y} = f^E(t, y), \quad y(t_0) = y_0
$$


with $f^B(t, y)$ taking the form as the right-hand-side of a second-order, linear, homogeneous ODE as described above. OSCODE will then dynamically switch between a Runge-Kutta solver and its asymptotic ARDC-based solver, at each step choosing the method which yields the larger step length (while keeping the local error within user tolerance).

The OSCODE framework is packaged with a built in Butcher tableau for a 6-stage, 5th order, and 4-stage, 4th order Runge-Kutta method, though users can supply their own Buthcer tableaus of any order.

<!--- OSCODE also supports unevenly spaced samples of $y$. What do you mean by this? --->

Runge-Kutta Solver:

Augmented Riccati Defect Correction Solver:

<!--- What do you want me to put here? A quick math overview of each? --->

## Guide-level Explanation

------------------------------------
Section Goals:

Explain the proposal as if it was already included in the project and you were teaching it to another Sundials programmer in the manual. That generally means:

- Introducing new named concepts.
- Explaining the feature largely in terms of examples.
- Explaining how Sundials programmers should *think* about the feature, and how it should impact the way they use the relevant package. It should explain the impact as concretely as possible.
- If applicable, provide sample error messages, deprecation warnings, or migration guidance.
- If applicable, describe the differences between teaching this to existing Stan programmers and new Stan programmers.

For implementation-oriented RFCs (e.g. for compiler internals), this section should focus on how compiler contributors should think about the change, and give examples of its concrete impact. For policy RFCs, this section should provide an example-driven introduction to the policy, and explain its impact in concrete terms.

<!--- Should the above be deleted? --->
------------------------------------

(Note: Much of this is taken from the ARKODE docs and modified for OSCODE)

<!--- I removed/modified many of the statements because they didn't apply for OSCODE.--->

The OSCODE infrastructure provides adaptive-step time integration modules for a class of nonstiff ordinary differential equations (ODEs) where the components may be highly oscillatory.  

OSCODE supports ODE systems posed in the form below

$$
M(t)\ddot{y} = f^O(t, y), \quad y(t_0) = y_0
$$

with

$$
f^O(t, y) = \begin{pmatrix}
y[1]\\
-y[0]w^2(t)-2\gamma(t)y[1]
\end{pmatrix}.
$$

Here, $f^O$ is a function with varying sized oscillatory regions. $t$ is the independent variable, $y$ is a set of two (complex-valued) dependent variables with one being the derivative of the other, $\omega$ and $\gamma$ are user specified, real-valued callables. The solver has the ability to switch dynamically between using a Runge-Kutta method and a Wentzel-Kramers-Brillouin method to advance the solution, while adaptively updating its stepsize.

### Adaptive single-step methods

The OSCODE framework is designed to support single step, IVP integration methods

$$
\ddot{y} = f(t ,y), y(t_0) = y_0
$$

The choice of step size is determined by the time-stepping method (based on user-provided inputs, typically accuracy requirements). However, users may place minimum/maximum bounds on if desired.

At each step, either the ARDC or RK solver is chosen for that step based on the relative difference in error between the two methods. For RK, the error is generated by the difference in the $N$th and $N-1$th order RK steps. The ARDC steps error for the function evaluation and the first derivative of $x$ with respect to time is defined as

$$
\delta x_{ARDC} = A_+ \delta f_+ + A_- \delta f_-
$$
$$
\delta f_\pm = f_\pm \sum_{i=0}^n \delta[S_i]_t^{t + h}
$$

$$
\delta \dot{x}\_{ARDC} = B\_+ \delta \dot{f}\_+ + B\_- \delta \dot{f}\_-
$$

$$
\delta \dot{f}\_\pm = \delta f\_\pm \frac{\dot{f}\_\pm}{f\_\pm}
$$

OSCODE's time stepping modules may be run in a variety of "modes":

- NORMAL: The solver will take internal steps until it has just overtaken a user-specified output time, $t_{out}$, in the direction of integration, i.e. $t_{n-1} < t_{out} \le t_n$ for forward integration, or $t_{n} < t_{out} \le t_{n-1}$ for backward integration. It will then compute an approximation to the solution $y(t_{out})$ by interpolation (using one of the dense output routines described in the section TBD).

- ONE-STEP: The solver will only take a single internal step $y_{n-1} \rightarrow y_n$ and then return control back to the calling program. If this step will overtake $t_{out}$ then the solver will again return an interpolated result; otherwise it will return a copy of the internal solution $y_n$

- NORMAL-TSTOP: The solver will take internal steps until the next step will overtake $t_{out}$. It will then limit this next step so that $t_n = t_{n-1} + h_n = t_{out}$, and once the step completes it will return a copy of the internal solution

- ONE-STEP-TSTOP: The solver will check whether the next step will overtake $t_{out}$ if not then this mode is identical to "one-step" above; otherwise it will limit this next step so that $t_n = t_{n-1} + h_n = t_{out}$. In either case, once the step completes it will return a copy of the internal solution $y_n$.
.

We note that interpolated solutions may be slightly less accurate than the internal solutions produced by the solver. Hence, to ensure that the returned value has full method accuracy one of the "tstp" modes may be used.

### Interpolation

As mentioned above, the time-stepping modules in OSCODE support interpolation of solution $y(t_{out})$ and derivatives $y^{(d)}(t_{out})$, where $t_{out}$ occurs within a completed time step from $t_{n-1} \rightarrow t_n$. Additionally, this module supports extrapolation of solutions and derivatives for outside this interval (e.g. to construct predictors for iterative nonlinear and linear solvers). To this end, OSCODE currently supports construction of polynomial interpolants $p_q(t)$ of polynomial degree up to $q=5$, although users may select interpolants of lower degree.

OSCODE provides three complementary interpolation approaches, both of which are accessible from any of the time-stepping modules: "Hermite" and "Lagrange". The Hermite approach has been included with ARKODE since its inception, and is more suitable for non-stiff problems; Lagrange is designed to provide increased accuracy when integrating stiff problems. These two methods are described in detail within the ARKODE documentation (link).

### User callable functions

```c++
/**
 * @param fe -- the name of the C function (of type `OSRhsFn()`)
 *   defining the explicit portion of the right-hand side function in
 *   `M(t)\, y'(t) = f^E(t,y) + f^I(t,y)`.
 * @param fi -- the name of the C function (of type `OSRhsFn()`)
 *       defining the implicit portion of the right-hand side function in
 *       `M(t)\, y'(t) = f^E(t,y) + f^I(t,y)`.
 * @param t0 -- the initial value of :math:`t`.
 * @param y0 -- the initial condition vector :math:`y(t_0)`.
 * @param sunctx -- the :c:type:`SUNContext` object (see :numref:`SUNDIALS.SUNContext`)
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

// Set order and tables for ARDC solver 
int OSCStepARDCSetTables(void* oscode_mem, int o, OSodeButcherTable Bi);

// If the user sets the ARDC tables and order they must also set the series and derivative functions.
int OSCStepARDCSetSeriesFuncs(void* oscode_mem, OSCddsFn dds, OSCdsiFn dsi, OSCdsfFn dsf);

// Set order and tables for RK solver 
int OSCStepRKSetTables(void* oscode_mem, int o, OSodeButcherTable Bi)

int OSCStepRKSetTableNum(void* oscode_mem, OSODE_ERKTableID etable)
```

Should we have methods for setting the time step adaptivity? Ex

```C
========================================================   ======================================  ========
Optional input                                             Function name                           Default
========================================================   ======================================  ========
Set a custom time step adaptivity function                 :c:func:`OSCStepSetAdaptivityFn()`      internal
Choose an existing time step adaptivity method             :c:func:`OSCStepSetAdaptivityMethod()`  0
Explicit stability safety factor                           :c:func:`OSCStepSetCFLFraction()`       0.5
Time step error bias factor                                :c:func:`OSCStepSetErrorBias()`         1.5
Bounds determining no change in step size                  :c:func:`OSCStepSetFixedStepBounds()`   1.0  1.5
Maximum step growth factor on convergence fail             :c:func:`OSCStepSetMaxCFailGrowth()`    0.25
Maximum step growth factor on error test fail              :c:func:`OSCStepSetMaxEFailGrowth()`    0.3
Maximum first step growth factor                           :c:func:`OSCStepSetMaxFirstGrowth()`    10000.0
Maximum allowed general step growth factor                 :c:func:`OSCStepSetMaxGrowth()`         20.0
Minimum allowed step reduction factor on error test fail   :c:func:`OSCStepSetMinReduction()`      0.1
Time step safety factor                                    :c:func:`OSCStepSetSafetyFactor()`      0.96
Error fails before MaxEFailGrowth takes effect             :c:func:`OSCStepSetSmallNumEFails()`    2
Explicit stability function                                :c:func:`OSCStepSetStabilityFn()`       none
========================================================   ======================================  ========

```

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
    - `std::complex<double> (OSCODEOmegaFn*)(realtype)`
3. Options for the rk solver
    - order
    - tolerances
    - butcher tables
    - Gauss - Lobatto weights
    - Integral Function
4. Options for the wkb solver
    - Order
    - order related weight matrices
    - Derivative and series functions
    - Function for Hermite and Lagrange
    - Gauss-Lobatto weights
5. Options for solver
    - initial conditions for the ODE, $ x(t) \f$, \f$ \frac{dx}{dt} $ evaluated at the start of the integration range
    - start of integration range
    - end of integration range
    - do_times timepoints at which dense output is to be produced. Must be sorted
    - Direction
    - order of ARDC approximation to be used
    - (local) relative tolerance
    - (local) absolute tolerance
    - initial stepsize to use

The implimentation can utilize the ARKODE explicit runga kutta solver, but we will have to write an implementation of the ARDC solver. For a first implementation I think it would be fine to only have the 6th order ARDC solver already written in C++, then in another PR expose functions that allow the user to change the ARDC solvers order.

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

The OSCODE solver is currently available as a python package and standalone header files.

## Papers

Are there any published papers or great posts that discuss this? If you have some relevant papers to refer to, this can serve as a more detailed theoretical background.

- Original paper: ["Efficient method for solving highly oscillatory ordinary differential equations
with applications to physical systems"](https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.2.013030)
- JOSS Paper on Python / C++ implementation: [Link](https://joss.theoj.org/papers/10.21105/joss.02830)
- C++ Implementation: [Link](https://github.com/fruzsinaagocs/oscode)
- Python Wrapper Docs: [Link](https://oscode.readthedocs.io/en/latest/)

## Unresolved Questions

- What parts of the design do you expect to resolve through the RFC process before this gets merged?

Whether OSCODE should be its own seperate module or whether it is possible to fit it in the ARKODE framework.

- What parts of the design do you expect to resolve through the implementation of this feature before stabilization?

I'd like to resolve what the user should have access to modify, for example the butcher tables in the ARDC method.

- What related issues do you consider out of scope for this RFC that could be addressed in the future independently of the solution that comes out of this RFC?

Any reworking of the algorithm itself would be out of scope of this document.
