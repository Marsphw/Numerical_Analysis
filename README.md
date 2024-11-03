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
