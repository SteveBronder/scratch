- Start Date: Aug 2nd 2023

# Summary
[summary]: #summary

This design document describes the motivation and architecture for adding the OSCODE numerical routine for the efficient solution of one dimensional, second order, ordinary differential equations with rapidly oscillating solutions. The proposed solver will will operate with differential equations of the form 

$$ \ddot{x} + 2 \gamma(t) \dot{x}(t) + \omega^2(t)x(t) = 0 $$

where $\gamma$ and $\omega$ may or may not be expressed as a closed form function of time.

Where equation expressed above can be reposed as a first order ODE

$$
\dot{y} = \begin{pmatrix}
\dfrac{\partial M_x}{\partial t}\\
\dfrac{\partial M_y}{\partial t}\\
\dfrac{\partial M_z}{\partial t}
\end{pmatrix}
$$

This new module will follow the same C based API with an underlying C++ implementation. 

# Motivation
[motivation]: #motivation

Oscillatory ODEs are commonly used in fields such as physics, finance, and engineering. Solutions to oscillatory problems are often bespoke and use anlytic approximations such as the Wentzel-Kramers-Brilouin (WKB) method. Giving researchers access to a general purpose solvers for solutions for ODEs of the form written above will allow them to do their research faster and with better accuracy.





# Guide-level explanation
[guide-level-explanation]: #guide-level-explanation



Explain the proposal as if it was already included in the project and you were teaching it to another Stan programmer in the manual. That generally means:

- Introducing new named concepts.
- Explaining the feature largely in terms of examples.
- Explaining how Stan programmers should *think* about the feature, and how it should impact the way they use the relevant package. It should explain the impact as concretely as possible.
- If applicable, provide sample error messages, deprecation warnings, or migration guidance.
- If applicable, describe the differences between teaching this to existing Stan programmers and new Stan programmers.

For implementation-oriented RFCs (e.g. for compiler internals), this section should focus on how compiler contributors should think about the change, and give examples of its concrete impact. For policy RFCs, this section should provide an example-driven introduction to the policy, and explain its impact in concrete terms.

# Reference-level explanation
[reference-level-explanation]: #reference-level-explanation

This is the technical portion of the RFC. Explain the design in sufficient detail that:

- Its interaction with other features is clear.
- It is reasonably clear how the feature would be implemented.
- Corner cases are dissected by example.

The section should return to the examples given in the previous section, and explain more fully how the detailed proposal makes those examples work.

The main driver of `OSCODE` will be the function `oscEvolve` with a similar signature to `arkEvolve`

```C
/*---------------------------------------------------------------
  oscEvolve:

  This routine is the main driver of OSCODE-based integrators.

  It integrates over a time interval defined by the user, by
  calling the time step module to do internal time steps.

  The first time that arkEvolve is called for a successfully
  initialized problem, it computes a tentative initial step size.

  arkEvolve supports two modes as specified by itask: OSC_NORMAL and
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
2. An interpolation schema for either or both of `w` and `g`
3. Options for the rk solver
  - tolerances
  - butcher tables
  - (P_dense?)
  - function `f` for going from second order ODE into a system of first order ODEs (?)
4. Options for the wkb solver
  - Order
  - Gauss-Lobatto weights
  - vandermonde iterpolation values (?) 
5. Options for solver
  - initial conditions for the ODE, \f$ x(t) \f$, \f$ \frac{dx}{dt} \f$ evaluated at the start of the integration range
  - start of integration range
  - end of integration range
  - do_times timepoints at which dense output is to be produced. Must be sorted
  - order of WKB approximation to be used
  - (local) relative tolerance
  - (local) absolute tolerance
  - initial stepsize to use





# Drawbacks
[drawbacks]: #drawbacks

Why should we *not* do this?

# Rationale and alternatives
[rationale-and-alternatives]: #rationale-and-alternatives

- Why is this design the best in the space of possible designs?
- What other designs have been considered and what is the rationale for not choosing them?
- What is the impact of not doing this?

# Prior art
[prior-art]: #prior-art

Discuss prior art, both the good and the bad, in relation to this proposal.
A few examples of what this can include are:

- For language, library, tools, and compiler proposals: Does this feature exist in other programming languages and what experience have their community had?
- For community proposals: Is this done by some other community and what were their experiences with it?
- For other teams: What lessons can we learn from what other communities have done here?
- Papers: Are there any published papers or great posts that discuss this? If you have some relevant papers to refer to, this can serve as a more detailed theoretical background.

This section is intended to encourage you as an author to think about the lessons from other languages, provide readers of your RFC with a fuller picture.
If there is no prior art, that is fine - your ideas are interesting to us whether they are brand new or if it is an adaptation from other languages.

Note that while precedent set by other languages is some motivation, it does not on its own motivate an RFC.
Please also take into consideration that rust sometimes intentionally diverges from common language features.

# Unresolved questions
[unresolved-questions]: #unresolved-questions

- What parts of the design do you expect to resolve through the RFC process before this gets merged?
- What parts of the design do you expect to resolve through the implementation of this feature before stabilization?
- What related issues do you consider out of scope for this RFC that could be addressed in the future independently of the solution that comes out of this RFC?
