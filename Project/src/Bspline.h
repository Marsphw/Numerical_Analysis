#ifndef BSPLINE_H
#define BSPLINE_H

#include <iostream>
#include "Spline.h"
#include "Polynomial.h"
#include "Equation_Solving.h"

class Bspline {

public:


private:
    int n, k, N; //s_n^k spline
    double boundary_conditions[5]; //boundary conditions, the number of elements equals to k
    std::vector<double> knots; //knot vector
    std::vector<Polynomial> basis_polynomials; //basis polynomials
    std::vector<double> coefficients; //coefficients of the spline
};

#endif