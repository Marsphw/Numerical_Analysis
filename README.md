# Numerical_Analysis
This is my repository for the course "Numerical Analysis" of Zhejiang University. 

To get the report and the results, you can input these commands in the terminal in the folder of every task:
1. ```make report``` to generate the report in the ```tex``` folder.
2. ```make run``` to compile and run the code in the ```code``` folder. (It should be mentioned that the ```code``` folder exists in the theoretical assignments.)
3. ```make clean``` to clean the temporary files, including the executable file and the report file.

In addition, there is a ```README.md``` in every folder of task, so you can get more information about it.
## Chapter 1: Solving Nonlinear Equations

This chapter is mainly about the numerical methods of solving nonlinear equations. The three methods are:
1. Bisection Method
2. Newton's Method
3. Secant Method

### Theoretical Problems

The answer can be seen in the [Chapter1_Theoretical](https://github.com/Marsphw/Numerical_Analysis/tree/main/Theroretical_Problem/Chapter1). This task concentrate on the theoretical analysis in the three methods, such as the convergence, the iteration formula and the number of steps required.

### Programming Assignments

The programming assignments are in the [Chapter1_Programming](https://github.com/Marsphw/Numerical_Analysis/tree/main/Programming_Assignments/Chapter1). They focus on the implementation of the three methods in C++, with the oriented object programming style.

## Chapter 2: Polynomial Interpolation 

This chapter focuses on the numerical methods of polynomial interpolation. *Interpolation* constructs new data points within the range of a discrete set of known data points. Given the excellent properties(including continuity, smoothness, and differentiability), we typically choose polynomial for interpolation. We have **the Lagrange formula, the Newton formula, the Hermite interpolation and the Bezier curve**.

### Theoretical Problems

The answer can be seen in the [Chapter2_Theoretical](https://github.com/Marsphw/Numerical_Analysis/tree/main/Theoretical_Problem/Chapter2). This task includes the theoretical analysis of the interpolation methods, such as the error analysis, the calculation of the interpolation polynomial, and some proof of the lemma from our textbook.

### Programming Assignments

The programming assignments are in the [Chapter2_Programming](https://github.com/Marsphw/Numerical_Analysis/tree/main/Programming_Assignments/Chapter2). The tasks are to implement the Newton formula, the Hermite interpolation, and the Bezier curve in C++, with necessary figures and explanations.

## Chapter 3: Splines

This chapter concentrates on the construction and prosperities of splines, including PP-form and B-form. **Splines** are piecewise-defined polynomial function that passes through a set of data points, with each segment between two points being a polynomial, and the entire curve is continuous.

### Theoretical Problems

The report can be seen in the [Chapter3_Theoretical](https://github.com/Marsphw/Numerical_Analysis/tree/main/Theoretical_Problem/Chapter3). The tasks include the calculation of splines, the proof of some equations and lemmas.

### Programming Project

The programming project needs us to achieve a package to implement the splines with PP-form and B-form, and give the following contents:
1. Implement the linear spline function $\mathbb{S}_1^0$;
2. Implement the cubic spline function $\mathbb{S}_3^2$ in PP-form with boundary conditions of periodic splines, natural splines and complete splines;
3. Implement the same function in B-form with the same boundary conditions;
4. Verify that we will get the same curve in PP-form and B-form when the knots and boundary conditions are the same;
5. Support the splines of arbitrary ranks and knots in B-form;
6. Implement the curve fitting in a plane plat;
7. Implement the curve fitting in a sphere.

## Chapter 4: Computer Arithmetic and Conditioning

This chapter helps us to understand the world of computer arithmetic, including the floating-point number systems, rounding error, accuracy and stability. In the computer, the number of numbers is limited by the restriction of precision.

### Theoretical Problems

The report can be found in the [Chapter4_Theoretical](https://github.com/Marsphw/Numerical_Analysis/tree/main/Theoretical_Problem/Chapter4). These tasks consider the convertion of FPNs, proof of some phenomena and the analysis for some functions.

## Chpater 5: Best Approximation and Least Square

This chapter deals with the best approximation problem in the 2-norm. We obtain the best approximation from Fourier expansions, and apply this abstract technique to continuous and discrete least squares.

### Theoretical Problems

The report can be found in the [Chapter5_Theoretical](https://github.com/Marsphw/Numerical_Analysis/tree/main/Theoretical_Problem/Chapter5). In this task, we finish the proof of some theroems and lemmas in detail and solve the specific problems of best approximation.

## Chapter 6: Numerical Integration and Differentiation

This chapter concentrates on the calculation of integration and differntiation in our computer. It introduces a lot of formulas to approximate the integration and differentiation with the Taylor expansion as foundation.

## Chapter 7: Numerical Solutions of Ordinary Differential Equations
