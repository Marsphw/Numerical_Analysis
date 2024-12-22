#ifndef BSPLINE_H
#define BSPLINE_H

#include <iostream>
#include "Spline.h"
#include "Polynomial.h"
#include "Equation_Solving.h"

class Bspline {

public:
    Bspline(int nn, int kk, int NN, double bc[]): n(nn), k(kk), N(NN) {
        for (int i = 0; i < k; i++) 
            boundary_conditions[i] = bc[i];
    }

    void construct_basis_polynomials(int i, int n) {
        
    }

private:
    int n, k, N; //S_n^k spline
    double boundary_conditions[5]; //boundary conditions, the number of elements equals to k
    std::vector<Polynomial> basis_polynomials; //basis polynomials
    std::vector<double> coefficients; //coefficients of the spline
};

#endif