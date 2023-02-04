---
layout: distill
title: Convex Solvers
description: A survey of the different classes of solvers for convex optimization problems
tags: math
giscus_comments: false
date: 2023-02-03

authors:
  - name: Govind Chari
    url: "https://govindchari.com/"
    affiliations:
      name: University of Washington Autonomous Controls Lab

bibliography: 2018-12-22-distill.bib

# Optionally, you can add a table of contents to your post.
# NOTES:
#   - make sure that TOC names match the actual section names
#     for hyperlinks within the post to work correctly.
#   - we may want to automate TOC generation in the future using
#     jekyll-toc plugin (https://github.com/toshimaru/jekyll-toc).
toc:
  - name: Introduction
    # if a section has subsections, you can add them as follows:
    # subsections:
    #   - name: Example Child Subsection 1
    #   - name: Example Child Subsection 2
  - name: Active Set Methods
  - name: Interior Point Methods (IPMs)
  - name: First Order Methods
  - name: Summary
  - name: References

# Below is an example of injecting additional post-specific styles.
# If you use this post as a template, delete this _styles block.
# _styles: >
#   .fake-img {
#     background: #bbb;
#     border: 1px solid rgba(0, 0, 0, 0.1);
#     box-shadow: 0 0px 4px rgba(0, 0, 0, 0.1);
#     margin-bottom: 12px;
#   }
#   .fake-img p {
#     font-family: monospace;
#     color: white;
#     text-align: left;
#     margin: 12px 0;
#     text-align: center;
#     font-size: 16px;
#   }

---

## Introduction

Convex optimization is a class of optimization concerned with minimizing a convex function over a convex set. 
One important feature of convex optimization is that any local minimum for a convex problem is the global minimum,
this means that the global minimum can be found very quickly. Mathamatically, a convex problem can be written as follows

$$
\begin{align*}
\min_{x} & \; f(x) \\
\textrm{s.t.} & \; h_{i}(x) = 0, \quad i = 1, \ldots, m \\
& \; g_{i}(x) \leq 0, \quad i = 1, \ldots, p
\end{align*}
$$

where $x \in \mathbb{R}^n$ is the optimization variable, $f(x)$ is a convex objective function, $h_{i}(x)=0$ are affine equality
constraints, and $g_{i}(x) \leq 0$ are convex inequality constraints.

This post focuses on stressing the intuition of different classes of convex solver and provides references for further reading at the end of each section.

***

## Active Set Methods
A constraint is said to be active at a point $x_0$ if $g_i(x_0) = 0$. We can define the optimal active set as the set of all constraints that are active at the optimal solution $x^*$. We can see that all the equality constraints will be in the optimal active set. 

Active set methods take advantage of the fact that equality constrainted problems are easier to solve than inequality constrainted problems. These methods start with a guess of the optimal active set and solve this equality constrained subproblem. Then it uses information from the solution of the subproblem, for example the sign of the dual variables, to add and remove constraints from the current guess of the active set.

If a good guess of the optimal active set is known, these methods can be be very fast and only take a handful of iterations. As a result these methods warmstart well and would be advantageous in applications like model predictive control since it is is unlikely that the optimal active set drastically changes between two solves.

The disadvantage of active set methods is that the theoretical worst-case runtime is exponential in the number of constraints since in the worst case, all combinations of constraints must be tested. 

One famous example of an active set algorithm is Simplex, which was invented by George Dantzig for Linear Programs (LPs). In Simplex, all iterates are vertices of the feasible set (which is a polytope), however this is not the case for quadratic programs (QPs) or any more complex optimization problem. Another active set solver is [qpOASES](https://github.com/coin-or/qpOASES).

Now we will see how equality constrained QPs can be solved in a simple way. For equality constrained QPs, the KKT conditions which are necessary and sufficient conditions for optimality are linear and thus they can be solved in one Newton step. We can write the equality constrained QP as follows where $Q>0$.

$$
\begin{align*}
\min_{x} & \; \frac{1}{2}x^\top Qx + q^\top x \\
\textrm{s.t.} & \; Ax=b\\
\end{align*}
$$

we can write the KKT conditions for this problem as follows, where $\lambda$ is a vector of dual variables

$$
\begin{align*}
Qx+q+A^\top \lambda &= 0 \\
Ax-b &= 0
\end{align*}
$$

We can see that this system of equations is linear in the primal and dual variables so we can find the optimal primal and dual solution by solving the following system
of equations

$$
\begin{bmatrix}
Q & A^\top \\
A & 0
\end{bmatrix}
\begin{bmatrix}
x \\
\lambda
\end{bmatrix} = 
\begin{bmatrix}
-q \\
b
\end{bmatrix}
$$

Thus we can see that solving the equality constrained quadratic program amounts to nothing more than solving a linear system.

To learn more about the details of active set methods reference Chapter 16 Section 5 of *Numerical Optimization* by Nocedal and Wright.

***

## Interior Point Methods (IPMs)

As the name suggests, interior point methods solve optimization problems in a way that the iterates lie in the interior of the feasible set. There are two main variants of IPMs: primal and primal-dual. In primal IPMs, we only compute iterates of the primal variables, and in primal-dual IPMs we compute iterates of both the primal and dual variables.

Primal IPMs make use of the fact that the following two problems are equivalent:

$$
\begin{align*}
\min_{x} & \; f(x) \\
\textrm{s.t.} & \; x \in \mathcal{D}\\
\end{align*}
$$

where $\mathcal{D}$ is some convex set and 

$$
\begin{align*}
\min_{x} & \; f(x) + \mathcal{I}_{\mathcal{D}}(x)\\
\end{align*}
$$


where $\mathcal{I}_{\mathcal{D}}(x)$ is the indicator function on $\mathcal{D}$ which is defined as

$$\mathcal{I}_{\mathcal{D}}(x) = 

\begin{cases}
0 & \text{if} \; x \in \mathcal{D} \\
\infty & \text{if} \; x \notin \mathcal{D}

\end{cases}$$

However, we cannot directly solve the minimization problem with the indicator function since it is nonsmooth at the boundary of the set $\mathcal{D}$, since it jumps from some finite value to infinity. This primal IPMs seek to replace the indicator function with a smooth approximation called a log-barrier function. This log-barrier function is roughly zero when $x \in \mathcal{D}$ and steeply approches infinity when you approach the boundary of $\mathcal{D}$. 

This steepness is controlled by a "barrier parameter." The steeper the barrier function is, the better it approximates the indicator function, but the less smooth and worse conditioned the minimization becomes. Initially the unconstrained problem is solved with a shallow barrier parameter and then is successively solved with steeper and steeper barrier parameters. This barrier can be thought of as a force-field that pushes the iterates away from the boundary of the feasible set and amount that this force field pushes the iterates is controlled by the barrier paramter.

Primal-Dual IPMs take a slightly different approach. They attempt to use Newton's method on the KKT conditions of the problem with some other fancy tricks such as taking a prediction step then a correction step which allows the algorithm to reuse the factorization of the KKT matrix. This famed trick is called Mehrotra's Predictor-Corrector. The Primal Dual IPM is famously used by SpaceX in their rocket landing algorithm [\[1,2\]](/blog/2023/optimization-algorithms/#references). 

One large drawback to using IPMs in real-time systems is the fact that they cannot be warmstarted which is a desirable property of real-time solvers.

One example of an IPM is [ECOS](https://github.com/embotech/ecos). 

To learn more about primal IPMs reference Chapter 11 in *Convex Optimization* by Boyd and Vanderberghe and to learn more about primal-dual IPMs reference Chapter 14 and Chapter 16 Section 6 in Nocedal and Wright.

***

## First-Order Methods

Both Active Set methods and IPMs typically rely on "second-order" information. Second order information is the information about the curvature of a function which is given by its second derivative, or for the multivariable case its Hessian. Using second order information allows these methods to converge quickly in few iterations, but each iteration requires the factorization of a large matrix which is a very expensive operation. 

The number of floating point operations for matrix factorization scales with $\mathcal{O}(n^3)$ and storing the Hessian scales with $\mathcal{O}(n^2)$ where $n$ is the number of variables in your problem. 

First Order methods on the other hand only use first-order information, which is information about the slope of a function which is given by its first derivative, or gradient in the multivariable case. First-order methods do not require matrix factorizations at each iteration and only require Matrix-vector multiplications. The number of floating point operations for matrix-vector multiplication scales with $\mathcal{O}(n^2)$, thus each iteration of a first-order method can be done quicker than an iteration of a second order method, but first order methods require more iterations to converge, since each iteration uses less information.

First order methods are more or less gradient descent algorithms with some modifications to handle constraints such as projections. For extremely large scale problems first order methods are preferable due to the high cost of factorizing and storing large matrices. However, even for medium sized problems we can gain a lot of performance from first-order methods by customizing the algorithm to the specific problem structure. For a detailed description of customization refer to [\[3\]](/blog/2023/optimization-algorithms/#references). 

All of this performance of first order methods does have some drawbacks. First order methods are extremely sensitive to ill conditioned objectives and badly scaled problem data. Thus an extrememly fast and robust implementation of a first-order method must scale the problem data, precondition the problem, and be customized to the problem structure. Without customization the algorithm will still be very fast (around the same speed as IPMs), but customization allows the full speed of the algorithm to be unlocked. To see some speed results of a particular first order algorithm called PIPG, refer to [\[3,4\]](/blog/2023/optimization-algorithms/#references). These papers have results that show that PIPG is faster than ECOS, SCS, MOSEK, Gurobi, and OSQP once customized, and preconditioned.

Some examples of first order solvers are [OSQP](https://osqp.org/), [SCS](https://github.com/cvxgrp/scs), and [PIPG]().

## Summary

Here we will quickly sumamrize the advantages and disadvantages of each method.

**Active Set**

Advantages: Easy to warmstart, fast if you have a good guess for the active set

Disadvantages: Worse case exponential runtime in the number of constraints, bad for large problems

**IPMs**

Advantages: Fast and robust for medium sized problems

Disadvantages: Bad for large problems, cannot warmstart, large code footprint

**First-Order**

Advantages: Small code footprint, good for large problems, can be good for medium sized problems with customization, ease of customization, easy to warmstart

Disadvantages: Highly sensitive to scaling and conditioning so they need scaling and preconditioning

## References

[1] *Blackmore, L. (2016). Autonomous precision landing of space rockets. In Frontiers of Engineering: Reports on Leading-Edge Engineering from the 2016 Symposium volume 46, 15–20*

[2] *J. Mattingley and S. Boyd. CVXGEN: A Code Generator for Embedded Convex Optimization. Optimization and Engineering, 13(1):1–27, 2012.*

[3] *Kamath, A.G., Elango, P., Yu, Y., Kim, T., Carson III, J.M., Mesbahi, M., and Açıkmeşe, B. (2023). Customized real-time first-order methods for onboard dual quaternion-based 6-dof powered-descent guidance. AIAA SciTech 2023 Forum.*

[4] *Yu, Y., Elango, P., Açıkmeşe, B., and Topcu, U., “Extrapolated Proportional-Integral Projected Gradient Method for
Conic Optimization,” arXiv preprint arXiv:2203.04188, 2022.*

***