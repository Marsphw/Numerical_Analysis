#ifndef EQUATION_SOLVING_H
#define EQUATION_SOLVING_H

#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include "Matrix.h"

class Equation {

public:
    Equation(int nn, Matrix<double> MM, std::vector<double> bb): n(nn), M(MM), b(bb) {}

    std::vector<double> solve() {
        std::vector<int> indx(n);
        for (int i = 0; i < n; i++) indx[i] = i;
        LU_decomposition(M, indx);
        return LU_solve(M, indx, b);
    }

    void swapRow(Matrix<double>& A, int i, int j) {
        for (int k = 0; k < A.getCols(); k++) 
            std::swap(A(i, k), A(j, k));
    }

    void LU_decomposition(Matrix<double>& A, std::vector<int>& indx) {
        int n = A.getRows();
        for (int i = 0; i < n; ++i) {
            double maxEl = std::abs(A(i, i));
            int maxRow = i;
            for (int j = i + 1; j < n; ++j) {
                if (std::abs(A(j, i)) > maxEl) {
                    maxEl = std::abs(A(j, i));
                    maxRow = j;
                }
            }

            if (maxRow != i) {
                swapRow(A, i, maxRow);
                std::swap(indx[i], indx[maxRow]);
            }

            for (int j = i + 1; j < n; ++j) {
                double c = A(j, i) / A(i, i);
                for (int k = i; k < n; ++k) {
                    A(j, k) -= c * A(i, k);
                }
            }
        }
    }

    std::vector<double> LU_solve(Matrix<double>& A, std::vector<int>& indx, std::vector<double>& b) {
        int n = A.getRows();
        std::vector<double> x(n);

        //Ly = b
        for (int i = 0; i < n; ++i) {
            int ip = indx[i];
            double sum = b[ip];            
            for (int j = 0; j < i; ++j) 
                sum -= A(i, j) * x[j];
            x[i] = sum / A(i, i);
        }

        //Ux = y
        for (int i = n - 1; i >= 0; --i) {
            double sum = x[i];
            for (int j = i + 1; j < n; ++j) 
                sum -= A(i, j) * x[j];
            x[i] = sum / A(i, i);
        }
        return x;
    }

private:
    int n;
    Matrix<double> M;
    std::vector<double> b;

};

#endif