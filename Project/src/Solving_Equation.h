#ifndef SOLVING_EQUATION_H
#define SOLVING_EQUATION_H

#include <iostream>
#include <cmath>
#include "lapacke.h"

class EquationSolver {

private:
    lapack_int N;
    double *A, *b;

public:
    EquationSolver(int n, double *a, double *rhs): N(n), A(a), b(rhs) {}

    double* solve() {
        lapack_int IPIV[N];
        lapack_int INFO = 0;
        /*INFO = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, A, N, IPIV);
        if (INFO != 0) {
            std::cout << "Error in LU factorization" << std::endl;
            return NULL;
        }
        INFO = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', N, 1, A, N, IPIV, b, N);
        if (INFO != 0) {
            std::cout << "Error in solving equation" << std::endl;
            return NULL;
        }*/
        INFO = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, 1, A, N, IPIV, b, N);
        if (INFO != 0) {
            std::cout << "Error in solving equation" << std::endl;
            return NULL;
        }
        return b;
    }

};

#endif