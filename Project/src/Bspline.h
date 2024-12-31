#ifndef BSPLINE_H
#define BSPLINE_H

#include <iostream>
#include "Spline.h"
#include "Polynomial.h"
#include "Solving_Equation.h"

class Bspline {

public:
    Bspline(int nn, int NN, std::vector<std::pair<double, double>> knots_in, double bc[] = 0): n(nn), N(NN), knots(knots_in) {
        for (int i = 0; i < n - 1; i++) 
            boundary_conditions[i] = bc[i];
    }

    void construct_basis_polynomials(int i, int ni) {
        if (ni == 0) {
            std::vector<double> coeffs{1.0};
            basis_polynomials[i][ni] = std::vector<Polynomial>{Polynomial(0, coeffs)};
            return ;
        }
        std::vector<Polynomial> result;
        if (i == 1 - n) 
            construct_basis_polynomials(i, ni - 1);
        construct_basis_polynomials(i + 1, ni - 1);

        double h1 = knots[1].first - knots[0].first;
        double h2 = knots[N - 1].first - knots[N - 2].first;

        std::vector<double> coeff1{-((i <= 0) ? knots[0].first + i * h1 : knots[i - 1].first) / (knots[i + ni].first - ((i <= 0) ? knots[0].first + i * h1 : knots[i - 1].first)), 1 / (knots[i + ni].first - ((i <= 0) ? knots[0].first + i * h1 : knots[i - 1].first))};
        Polynomial p1(1, coeff1);

        std::vector<double> coeff2{((i + ni + 1 >= N - 1) ? knots[N - 1].first + (i + ni + 1 - N + 1) * h2 : knots[i + ni + 1].first) / (knots[N - 1].first + (i + ni + 1 - N + 1) * h2 - knots[i].first)};
        Polynomial p2(1, coeff2);

        std::vector<double> temp{0.0};
        result.push_back(Polynomial(0, temp));
        for (int k = 0; k < ni; ++k) {
            result[k] = result[k] +  p1 * basis_polynomials[i][ni - 1][k];
            result.push_back(p2 * basis_polynomials[i + 1][ni - 1][k]);
        }

        basis_polynomials[i][ni] = result;
    }

    Spline complete_cubic_Bspline() {
        for (int i = 1 - n; i < N; ++i)
            construct_basis_polynomials(i, n);

        Matrix<double> coeff_matrix(N, N);

        coeff_matrix(0, 0) = 1.0/3 + (knots[2].first - knots[0].first) / (2 * (knots[2].first + knots[1].first - 2 * knots[0].first));
        coeff_matrix(0, 1) = (knots[1].first - knots[0].first) / (knots[2].first + knots[1].first - 2 * knots[0].first);

        for (int i = 1; i < N - 1; ++i) {
            coeff_matrix(i, i - 1) = basis_polynomials[i - 1][n][n - 1](knots[i].first);
            coeff_matrix(i, i) = basis_polynomials[i][n][n - 2](knots[i].first);
            coeff_matrix(i, i + 1) = basis_polynomials[i + 1][n][n - 3](knots[i].first);
        }

        coeff_matrix(N - 1, N - 2) = (knots[N - 1].first - knots[N - 2].first) / (2 * knots[N - 1].first - knots[N - 2].first - knots[N - 3].first);
        coeff_matrix(N - 1, N - 1) = 1.0/3 + (knots[N - 1].first - knots[N - 3].first) / (2 * (2 * knots[N - 1].first - knots[N - 2].first - knots[N - 3].first));

        std::vector<double> rhs(N);

        rhs[0] = knots[0].second + boundary_conditions[0] * (knots[1].first - knots[0].first) / 3;
        for (int i = 1; i < N - 1; ++i)
            rhs[i] = knots[i].second;
        rhs[N - 1] = knots[N - 1].second - boundary_conditions[1] * (knots[N - 1].first - knots[N - 2].first) / 3;

        Equation linear_system(N, coeff_matrix, rhs);
        std::vector<double> solution = linear_system.solve();

        std::vector<Polynomial> polynomials(N);
        for (int i = 0; i < N; ++i) {
            ;
        }
        Spline spline(N, knots, polynomials);
    }

    Spline natural_cubic_Bspline() {
        for (int i = 1 - n; i < N; ++i)
            construct_basis_polynomials(i, n);
        
    }

    Spline period_Bspline() {
        for (int i = 1 - n; i < N; ++i)
            construct_basis_polynomials(i, n);
    }

private:
    int n, N; //S_n^(n - 1) spline
    double boundary_conditions[5]; //boundary conditions, the number of elements equals to k
    std::vector<std::pair<double, double>> knots; //t_0 to t_{N - 1}
    std::vector<std::vector<std::vector<Polynomial>>> basis_polynomials; //basis polynomials
    std::vector<double> coefficients; //coefficients of the spline, a_i
};

#endif