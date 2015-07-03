# OBA
Author: [Nitish Shirish Keskar](http://users.iems.northwestern.edu/~nitish/)


**OBA** is a second-order method for convex L1-regularized optimization with active-set prediction. 
OBA belongs to the family of _Orthant-Based_ methods (such as OWL) and uses a selective-corrective mechanism which brings about increased efficiency and robustness. 

## Features
The OBA package

* allows for solving general convex L1-regularized problems including Logistic Regression and LASSO.
* is written in pure-MATLAB with minimal dependencies and emphasizes simplicity and cross-platform compatibility. 
* includes both Newton and quasi-Newton options for the proposed algorithm.


## Usage Guide

The algorithm can be run using the syntax 

```
x = OBA(funObj,lambda,[options]);
```

Here,
* `funObj` is an object with member functions for computing the function, gradient and Hessian-vector products at the iterates. Logistic Regression and LASSO classes are provided with the package. The file `funTemplate.m` can be used as a base for designing a custom function.
* `lambda` is the positive scalar for inducing sparsity in the solution.
* `options` is an optional argument for changing the default parameters used in OBA. For ease of use, the user can generate the default options struct using `options=GenOptions()` and change the parameters therein before passing it to OBA.

The parameters and their default values are
*      `options.optol`: termination tolerance
          (default: 1e-6)
*      `options.qn`: Quasi-Newton, 0 (Newton's Method), or 1 (quasi-Newton)
          (default: 0)
*      `options.mem_size`: quasi-Newton memory size
          (default: 20)
*      `options.maxiter`: max number of iterations
          (default: 1000)
*      `options.printlev`: print level, 0 (no printing) or 1
          (default: 1)
*      `options.CGtol`: CG termination tolerance (for Newton's Method)
          (default: 1e-1)
*      `options.maxCGiter`: max number of CG iterations (Newton's Method)
          (default: 1000).

For a detailed documentation of OBA and its associated functions, use `help OBA`.

## Citation
If you use OBA for your research, please cite the following pre-print

A Second-Order Method for Convex L1-Regularized Optimization with Active Set Prediction : http://arxiv.org/abs/1505.04315













