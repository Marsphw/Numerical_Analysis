#ifndef SOLVING_EQUATION_H
#define SOLVING_EQUATION_H

#include <iostream>
#include <cmath>
#include "lapacke.h"

class EquationSolver {

private:
    int N;
    double *A, *b, *x;

public:
    EquationSolver(int n, double *a, double *rhs): N(n), A(a), b(rhs), x(new double[n]) {}

    double* solve() {
        int IPIV[N];
        for (int i = 0; i < N; i++)
            IPIV[i] = 0;
        int INFO = 0;
        INFO = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, A, N, IPIV);
        if (INFO != 0) {
            std::cout << "Error in LU factorization" << std::endl;
            return NULL;
        }
        INFO = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', N, 1, A, N, IPIV, b, N);
        if (INFO != 0) {
            std::cout << "Error in solving equation" << std::endl;
            return NULL;
        }
        for (int i = 0; i < N; i++)
            x[i] = b[i];
        return x;
    }

};

#endif