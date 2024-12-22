#ifndef PPFORM_H
#define PPFORM_H

#include <iostream>
#include "Polynomial.h"
#include "Equation_Solving.h"
#include "Spline.h"

class PPForm {

public:
    PPForm(int nn, int kk, int NN, std::vector<knot> knots_in, double bc[] = 0): n(nn), k(kk), N(NN), knots(knots_in) {
        for (int i = 0; i < k; ++i) 
            boundary_conditions[i] = bc[i];
    }

    Spline linear_splines() {
        std::vector<Polynomial> S;
        for (int i = 0; i < N - 1; ++i) {
            double a = (knots[i + 1].value[1] - knots[i].value[1]) / (knots[i + 1].value[0] - knots[i].value[0]);
            double b = knots[i].value[1] - a * knots[i].value[0];
            std::vector<double> coeffs = {b, a};
            Polynomial P(1, coeffs);
            S.push_back(P);
        }
        Spline result(N, knots, S);
        return result; 
    }
    
private:
    int n, k, N; //S_n^k spline
    double boundary_conditions[5]; //boundary conditions, the number of elements equals to k
    std::vector<knot> knots; //knot vector
};

#endif