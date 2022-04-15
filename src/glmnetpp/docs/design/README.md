# Design Documentation for `glmnetpp`

## Table of Content
- [Motivation](#motivation)
- [Current Design](#current-design)
  - [Library Structure](#library-structure)
  - [Elnet Driver](#elnet-driver)
  - [Elnet Path](#elnet-path)
    - [Why is there an `ElnetPathCRTPBase` and an `ElnetPathBase`?](#why-is-there-an-elnetpathcrtpbase-and-an-elnetpathbase)
    - [Is CRTP really needed?](#is-crtp-really-needed)
    - [What are these `...Pack` nested class templates?](#what-are-these-pack-nested-class-templates)
  - [Elnet Point](#elnet-point)
    - [Why do we need yet another policy for `ElnetPointInternalPolicy`?](#why-do-we-need-yet-another-policy-for-elnetpointinternalpolicy)
    - [Why does `ElnetPointCRTPBase` inherit the internal policy?](#why-does-elnetpointcrtpbase-inherit-the-internal-policy)
  - [Elnet Point Internal](#elnet-point-internal)
    - [What is `ElnetPointInternalStaticBase`?](#what-is-elnetpointinternalstaticbase)
- [Other Design Ideas](#other-design-ideas)
  - [Policy-Based Design](#policy-based-design)
  - [Dynamic Polymorphism (Virtual Mechanism)](#dynamic-polymorphism-virtual-mechanism)

This is a documentation for the design choices of the C++ implementation of `glmnet` original Fortran routines, `glmnetpp`.

## Motivation

For context, [`glmnet5dpclean.m`](../../../../inst/mortran/glmnet5dpclean.m) 
(or [`glmnet5dpclean.f`](../../src/legacy/glmnet5dpclean.f), the generated Fortran code)
is the MORTRAN code that was the workhorse for `glmnet` up until version `4.1-2`.
This code contains all the implementation of the IRLS algorithm for the various GLMs (Gaussian, Binomial, Poisson, etc.).
As seen in the MORTRAN code, the biggest drawback with the code is the difficulty in maintenance.
We believe that this backend code of `glmnet` 
should arguably be the _easiest_ to maintain, if anything,
as it is the most crucial part of the library.
Hence, the primary goal of this C++ implementation is to use better design principles 
to ease maintainability and extendability (see [Current Design](#current-design) for more detail).
We note in passing that in the process of converting the MORTRAN code, however, to our surprise, 
we actually found very tangible speed-ups as well
(at least 2x speed-up for dense versions and around 20% speed-up for sparse versions).

The biggest drawback with the MORTRAN code is the large code duplication.
All versions (for different GLMS) of the IRLS algorithm are _essentially_ the same algorithm,
yet the algorithm is copied and pasted as separate subroutines for each GLM and dense/sparse version,
leading to about 16 copies in total.
While this is tough for maintenance and future extendability, 
we should note that it may be more beneficial performance-wise.
A refactored code would potentially abstract a common block into a function
and reuse this function in the different versions, which would require a function call.
If the Fortran compiler cannot inline this function, we would pay a large cost in performance
if that function call were made frequently.
However, a copy-paste would be one way to manually force inlining of the function.
In C++, we have ways to strong-inline functions (force the compiler to inline a function if possible)
(see [macros.hpp](../../include/glmnetpp_bits/util/macros.hpp))
and we are unaware of similar capabilities of the Fortran compiler.

The second drawback of the MORTRAN code is that the linear algebra routines are not as optimized 
as those provided by some matrix libraries in C++.
Our implementation relies on the widely-used, industry-standard
matrix library, [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page).
We noticed that `Eigen` allows for more vectorized code, leading to performance boosts.

From the above, the problem is clear:
_how do we write a refactored library where we reuse abstracted code 
as much as possible to remove code duplications,
write cache-friendly code, and use more vectorized code
without incurring the cost of function calls and generally avoiding costs that come from abstraction?_

Let's dive into the problem a little bit more.
At first glance, it may seem like a simple task to refactor the MORTRAN code,
however, the abstraction is quite subtle.
Usually, when we refactor, we have the following problem:
```cpp
void foo1() 
{
    // some stuff..
    // BEGIN COMMON BLOCK
    // ...
    // END COMMON BLOCK
    // some stuff..
}

void foo2() 
{
    // some different stuff..
    // BEGIN COMMON BLOCK
    // ...
    // END COMMON BLOCK
    // some different stuff..
}
```
and we would like to abstract the common block into a separate function from two macroscopic functions like:
```cpp
void common_block() { /*...*/ }
void foo1() 
{
    // some stuff..
    common_block();
    // some stuff..
}

void foo2() 
{
    // some different stuff..
    common_block();
    // some different stuff..
}
```
Indeed, we do see these kinds of problems where the common block could be 
the formula for gradient computation,
residual updates, and beta updates.

This is, indeed, a simple task, however, this is not the __only__ problem we see in the MORTRAN code.
We also have the following kind of problem:
```cpp
void foo(...) 
{
    // 1. set-up quantities needed for the elastic net path-solver

    // 2. loop through each lambda in the regularization path
    for (int m = 0; m < max_n_lambdas; ++m) {
        // 3. compute the current lambda at iteration m
        // 4. solve the elastic net problem at the current lambda
        // 5. save the result into output
    }

    // 6. do any final computations after solving for the whole path
}
```
Essentially, there is a _common macroscopic function_ `foo` and some of the _internal details_ (numbered steps) are different
depending on the GLM and dense/sparse storage of the data matrix.
Moreover, step 4 also has further common code such as the weighted-least-squares (WLS) algorithm of the form:
```cpp
void wls(...)
{
    while (1) {
        double convergence_measure = 0.0;
        // 1. iterate through each feature and do coordinate descent.
        for (int k = 0; k < p; ++k) {
            // update beta...
            if (new_beta == old_beta) continue;
            // update convergence_measure...
            // update other quantities...
        }

        // 2. check for convergence
        if (convergence_measure < tolerance) {
            // check KKT...
            if (kkt_passed) return;
            continue; 
        }

        // 3. do coordinate descent on only the active variables until convergence
    }
}
```
All versions of the subroutines reuse this general algorithm.
Only the exact details of how we update `beta`, `convergence_measure`, and `other quantities` are different
across the different versions of the algorithm.

In the [next section](#current-design),
we explain our current design of the library and how it tries to solve this very problem.

## Current Design

In this section, we go over the current design of `glmnetpp`.
First, we summarize the library structure at a high-level.
Then, we delve into each of the key components of the library.
Of course, we cannot explain every minute detail,
so we only explain the big _key_ design choices. 
We do not discuss the correctness of the implementation,
as we did not change the algorithm itself (except bug fixes).

### Library Structure

The following file structure outlines the library structure:
```
glmnetpp
    \_ include
        \_ glmnetpp_bits
            \_ elnet_driver
            \_ elnet_path
            \_ elnet_point
                \_ internal
            \_ util
            \_ elnet_driver.hpp
            \_ elnet_path.hpp
            \_ elnet_point.hpp
        \_ glmnetpp
    \_ src
        \_ legacy
```
Recall that `glmnetpp` is a header-only template library, 
so the whole library is defined in `include`.
`src` exists purely for legacy reasons to store the old Fortran code.
`glmnetpp_bits` contains the bits and pieces of the library
and the full library is exported into the final header file `glmnetpp`.
Inside `glmnetpp_bits`, there are four major components:

- [**elnet_driver**](#elnet-driver): contains driver classes that are at the highest-level of the hierarchy of classes.
    The main job is to delegate to the correct elastic net path-solver depending on the input.
    The exported header file is `elnet_driver.hpp`.
- [**elnet_path**](#elnet-path): contains elastic net path-solver classes for the different GLMs (dense and sparse).
    These classes are responsible for solving the elastic net problem for the whole regularization path.
    The exported header file is `elnet_path.hpp`.
- [**elnet_point**](#elnet-point): contains generic elastic net point-solver classes for the different GLMs (dense and sparse).
    These classes are responsible for solving the elastic net problem for a _fixed_ regularization parameter.
    The exported header file is `elnet_point.hpp`.
- [**elnet_point/internal**](#elnet-point-internal): contains the internal details that define the GLM-specific routines
    for the point-solvers.

Finally, `util` contains some utility functions/macros that are not intended to be exported for direct use,
but rather as tools in the implementation of the library.

### Elnet Driver

All drivers have very similar implementation details, so we only discuss one of the drivers, `gaussian.hpp`.
These drivers mimic the Fortran functions such as `elnet, lognet, fishnet`, 
the highest-level functions that get exported to `R`.
The class declaration is of the form
```cpp
template <>
struct ElnetDriver<util::glm_type::gaussian>
    : ElnetDriverBase
{
private:
    static constexpr util::glm_type glm = util::glm_type::gaussian;
    using mode_t = util::mode_type<glm>;

public:
    template <...>
    void fit(...);
};
```
that is, we have a single class template `ElnetDriver`
which takes in an `enum class` value of type `util::glm_type`
defined in [types.hpp](../../include/glmnetpp_bits/util/types.hpp).
The idea is that we only specialize this class template only
for particular GLMs,
since logically, there is no "default" definition of a driver 
for a general GLM.
The same idea holds in `elnet_path` and `elnet_point`.
Note that for each `glm_type`, there is an associated `Mode` class
which defines flags to indicate further differences in the algorithm
for a given GLM.
As an example, we will see that Gaussian GLM 
allows for two different algorithms: naive and covariance.
These flags along with dense/sparse storage type of the data matrix
uniquely define each of the different versions of IRLS algorithm.

The `fit` function simply prepares the inputs 
for the path solvers such as 
standardizing features and normalizing weights if necessary.
Note that `fit` is heavily templatized! 
We could have easily written it non-templated and take in
specific `Eigen` matrix objects, since that is the only use-case currently.
However, looking ahead, we plan to export `glmnetpp` not just to `R`
but to other higher-level languages such as `Python` and `Matlab`.
It is not clear what objects we would need to pass to the `fit` function
in that case, as these languages may not support bindings for `Eigen` objects.
The only cost in templatizing is compile-time,
which is abundant for our use-case since we build this library once
when exporting in `R`.

`fit` also checks _at compile time_ whether `X` data matrix is dense or sparse.
Hence, at compile-time, we choose the correct version of the path solver.
The technical details of this "choosing" is slightly nuanced 
since we must be able to support C++14 standard and do not have the convenient `if constexpr` in C++17.
In `gaussian.hpp`, you will see `FitPathGaussian` class template in `namespace details`,
which is a metaprogramming tool to choose the dense or sparse version.
The problem is that if we write something like:
```cpp
template <...>
void fit(...) 
{
    constexpr bool is_dense = ...;
    if (is_dense) {
        call_dense_version(...);
    } else {
        call_sparse_version(...)
    }
}
```
the compiler has to be able to generate the code for _both_ sparse and dense versions 
_for the same given input_.
For example, if `X` were sparse, the compiler would still need to be able to generate code for
`call_dense_version` with `X`, which is not going to be case for us since
there are, naturally, operations that are not defined for sparse matrices
but are defined for dense matrices that get used in `call_dense_version` (such as `operator+=`).
Instead we add a layer of dispatch with `FitPathGaussian` like
```
FitPathGaussian<is_dense>::eval(...);
```
such that the compiler will first choose the specialization of `FitPathGaussian` depending on the boolean `is_dense`
and then generate the code for `eval` only for that specialization.

### Elnet Path

`ElnetPath` is a class template responsible for fitting an elastic net model for the full regularization path.
The declaration is as follows:
```cpp
template <util::glm_type glm
        , util::mode_type<glm> mode
        , class ElnetPointPolicy=ElnetPoint<glm, mode> >
struct ElnetPath;
```
As before, `glm` denotes the GLM enum value and `mode` is the `glm` corresponding mode value.
The third template parameter `ElnetPointPolicy` is the type of the point-solver.
A priori, `ElnetPath` does not depend on any specific implementation of a point-solver.
By default, we use our own implementation given in `ElnetPoint<glm, mode>`,
but we may substitute that type with a different one, so long as it satisfies a certain interface
(for now, read the code to see how `ElnetPoint` is used in `ElnetPath`, TODO: maybe write-up an interface).
The structure is hopefully clear: the most derived classes are the header files that don't end with `_base.hpp`
and the sparse versions are prefixed with `sp_`.
Each GLM has its own base class, and all such classes derive from the base class(es) in `base.hpp`.

#### Why is there an ElnetPathCRTPBase and an ElnetPathBase?

[Curiously Recurring Template Pattern (CRTP)](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern#:~:text=The%20curiously%20recurring%20template%20pattern,form%20of%20F%2Dbounded%20quantification.)
is a design pattern where a base class template takes in the derived type as a template parameter
and casts itself to that derived type to access its (public) interface.
A derived class `D` is expected to inherit from the CRTP base class `B<D>`, where `B` is the CRTP base class template.
As a result, for all derived classes `D_i` that inherit from `B<D_i>`, they all inherit from __different instances__ of `B`,
that is, `B<D_i>` are all __different class types__.
So, if there are any common functions defined in the class template `B` _that do not actually depend on the derived type_,
the compiler generates _copies_ of them in every `B<D_i>`, which leads to both increase in compile-time and code-bloat.
__Only__ the code that depends on the derived class, i.e. those that __must__ be differently generated per derived class,
should be defined inside `B`.
The rest should be in a non-CRTP base class (like `ElnetPathBase`) so that the compiler only generates one version of these member functions
for all derived classes.

#### Is CRTP really needed?

CRTP is not necessary for our purposes, but it does encapsulate and abstract very well so it is convenient.
All the CRTP base member functions could have been free function templates,
which would've also reduced compile-time,
but it does expose this function to the user technically, which isn't something we need to do.
We viewed these members as implementation detail and wanted to encapsulate them in a class.

#### What are these `...Pack` nested class templates?

All `ElnetPath` classes define three "pack" class templates: `FitPack, PathConfigPack, PointConfigPack`.
They all perform one job: pack a bunch of arguments into a single object.

- `FitPack` packs all the arguments from the user into a single object.
- `PathConfigPack`: packs any configuration variables needed for the path-solver.
- `PointConfigPack`: packs any configuration variables needed for the point-solver.

Not surprisingly, these packs have an inheritance hierarchy as there are many shared quantities across the different versions.
We decided to nest them rather than having them as separate classes because:

- Logically, they should not exist outside of `ElnetPath` -- they are structures that only have meaning in the context of `ElnetPath`.
- Their hierarchy structure directly follows that of `ElnetPath` and its base class(es), so it is convenient to nest.

But why create these packs in the first place? 
Let's consider the `FitPack` as an example.
All versions of the path solvers have different set of inputs.
For example, some take in an `offset` parameter whereas some don't.
Some take in a residual vector `r` where as others take instead a gradient vector `g`.
And then there are some that take both `r` and `g`.
It would be a massive headache to find a super-set of these arguments
and write our path-solvers assuming this one big set of arguments for the following reasons:

- We would have to dynamically check which arguments are used.
- It would make it difficult to see which version depends on what set of arguments.

Hence, we resolve this at compile-time by having each version define its own set of `FitPack`,
containing __only the arguments it needs__.
We now have a generic interface where the path-solvers simply work with one `FitPack`-like object.
See [base.hpp](../../include/glmnetpp_bits/elnet_path/base.hpp) to see how this generic `fit` function works
in `ElnetPathCRTPBase`.

Similarly, `PathConfigPack` and `PointConfigPack` define each version-specific set of values
that are needed for the path-solver and point-solver, respectively.

### Elnet Point

`ElnetPoint` is a class template responsible for fitting an elastic net model for a fixed regularization parameter.
The declaration is as follows:
```cpp
template <util::glm_type glm
        , util::mode_type<glm> mode
        , class ElnetPointInternalPolicy=ElnetPointInternal<glm, mode> >
struct ElnetPoint;
```
A priori, `ElnetPoint` does not depend on the internal implementation for a point solver.
By default, we use our policy `ElnetPointInternal` defined in [elnet_point/internal](../../include/glmnetpp_bits/elnet_point/internal),
but we may swap it for a different class so long as it satisfies the interface needed by `ElnetPoint`
(for now, read the code to see how `ElnetPoint` interacts with `ElnetPointInternal`.
Note that different specializations require a different interface).

#### Why do we need yet another policy for ElnetPointInternalPolicy?

It may seem like an over-generalization to have yet another policy for the internal implementation detail,
but once we explain what we mean by "internal", the design will make more sense.
As explained in [Motivation](#motivation),
there are two types of abstraction problems:

- Generalizing a small common block that is reused in many bigger functions.
- Generalizing the control-flow of many smaller functions, which can differ across different versions.

`ElnetPoint` aims to solve the second problem.
As an example, in [base.hpp](../../include/glmnetpp_bits/elnet_point/base.hpp),
`ElnetPointCRTPBase` contains some common routines used by all point solvers.
For example, the `fit` function is really the coordinate-descent algorithm.
Note that `fit` simply "manages" the smaller functions such as `increment_passes`,
`coord_desc_reset`, `for_each_with_skip`, `all_begin`, etc.
In particular, it defines the order and conditions under which these smaller functions are to be executed.
However, `fit` __does not__ rely on the __implementation-specific detail for the smaller functions__.
It is precisely these smaller functions that we refer to as "internal".
For each different GLM, we would obviously want a different internal policy 
that defines each of these smaller functions,
but for future development as well, we may want to try different internal policies for the same GLM.
This design grants us this flexibility.

#### Why does `ElnetPointCRTPBase` inherit the internal policy?

This was purely out of coding convenience.
The internal policies have constructors that are very cumbersome to type out,
so it was very convenient to simply inherit the policy and expose its constructor.
Because this doesn't break any logic, we stuck with this design.

### Elnet Point Internal

`ElnetPointInternal` is a class template responsible for defining the internal details for each `ElnetPoint` specialization.
The declaration is as follows:
```cpp
template <util::glm_type g
        , util::mode_type<g> mode
        , class ValueType = double
        , class IndexType = int
        , class BoolType = bool>
struct ElnetPointInternal;
```
Currently, the only special behavior that occurs is when `BoolType` is `bool`.
Internally, we keep an inclusion vector that contains which features may be included in the model
and which are always forbidden from being included.
If `BoolType` is `bool` we use `std::vector<bool>` to represent this vector,
as it is a more specialized and performant data structure.
Otherwise, we simply use an `Eigen::Matrix<BoolType, -1, 1>` (column vector).

This class simply tries to solve the first problem in the 
[previous section](#why-do-we-need-yet-another-policy-for-elnetpointinternalpolicy):
generalizing the smaller common blocks.
These classes contain the meat of the algorithm.
They were written with the intention that the members would be called in a certain order
by `ElnetPoint` classes, so when trying to read the code,
__read it together with the corresponding `ElnetPoint`__.

Currently, the most derived classes really only contain logic related to the `X` data matrix,
i.e. defining dense or sparse routines.
Otherwise, the two versions reuse the rest of the logic defined in their GLM-specific base class.
TODO: some work remains in abstracting the sparse logic across the different versions.

A lot of these base classes will reuse static members or further base class members.
Hopefully, this sheds some light into how much was being duplicated in the original MORTRAN code!

#### What is `ElnetPointInternalStaticBase`?

For encapsulation purposes and to associate routines under a certain name,
we place commonly used routines as static member functions in `ElnetPointInternalStaticBase`.
Other base classes may define more static members on top of these for the derived classes `ElnetPointInternal`.
The point of these static members is mainly to define formulas once.
A very common source of bugs stems from copying and pasting a common block multiple times.
Once an edit must be made in one, we have to manually update the copied-and-pasted places as well.
Our approach removes this problem from ever happening.
See the following issues on GitHub regarding such bugs:
- [Issue 39](https://github.com/trevorhastie/glmnet/issues/39)
- [Issue 40](https://github.com/trevorhastie/glmnet/issues/40)
- [Issue 41](https://github.com/trevorhastie/glmnet/issues/41)
- [Issue 42](https://github.com/trevorhastie/glmnet/issues/42)
- [Issue 43](https://github.com/trevorhastie/glmnet/issues/43)

## Other Design Ideas

The current design of `glmnetpp` is by no means a perfect design.
In this section, we go over some other ideas we thought about in the process.

### Policy-Based Design

Another idea we thought about was to use policy-based design.
We noticed that in [elnet_point/internal](../../include/glmnetpp_bits/elnet_point/internal),
we're borrowing from multiple base class members sometimes (such as Gaussian multi-response and Binomial multi-response grouped).
It is possible to remove this intricate dependency by defining a set of policies for each `ElnetPointInternal`.
We didn't go with this idea because it felt like the number of policies could get out-of-hand a priori.

### Dynamic Polymorphism (Virtual Mechanism)

In many places, we used CRTP and compile-time tricks to abstract away certain logic,
such as the `Pack` objects in `ElnetPack`.
However, one could have also used dynamic polymorphism without relying on template magic.
The big issue we saw was that dynamic polymorphism completely removes the benefit of inlining for us, which is detrimental!
Although it reduces binary size and compile-time,
once again, compile-time is abundant for us and code-bloat is not an issue
because this is not meant to be a generic library used by others for larger programs.
We only compile our export functions once using this library and create a dynamic library.
