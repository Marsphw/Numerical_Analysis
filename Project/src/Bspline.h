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
            basis_polynomials[i].push_back(std::vector<Polynomial>{Polynomial(0, coeffs)});
            return ;
        }
        std::vector<Polynomial> result;
        if (i == 1) 
            construct_basis_polynomials(i, ni - 1);
        construct_basis_polynomials(i + 1, ni - 1);

        std::vector<double> coeff1{-knots[i - 1].first / (knots[i + ni - 1].first -  knots[i - 1].first), 1.0 / (knots[i + ni - 1].first - knots[i - 1].first)};
        Polynomial p1(1, coeff1);
        
        std::vector<double> coeff2{knots[i + ni].first / (knots[i + ni].first - knots[i].first), -1.0 / (knots[i + ni].first - knots[i].first)};
        Polynomial p2(1, coeff2);
        
        std::vector<double> temp{0.0};
        result.push_back(Polynomial(0, temp));
        for (int k = 0; k < ni; ++k) {
            result[k] = result[k] +  p1 * basis_polynomials[i][ni - 1][k];
            result.push_back(p2 * basis_polynomials[i + 1][ni - 1][k]);
        }
        basis_polynomials[i].push_back(result);
    }

    Spline linear_Bspline() {
        knots.insert(knots.begin(), std::make_pair(2 * knots[0].first - knots[1].first, 0.0));
        knots.push_back(std::make_pair(2 * knots[N].first - knots[N - 1].first, 0.0));
        basis_polynomials.reserve(N + 3);
        for (auto& second_element: basis_polynomials) 
            second_element.reserve(n + 1);
        for (int i = 1; i <= N; ++i)
            construct_basis_polynomials(i, n);
        std::vector<Polynomial> polynomials;
        polynomials.push_back(Polynomial(0, {0.0}));
        polynomials.push_back(Polynomial(0, {0.0}));
        for (int i = 1; i <= N; ++i) {
            Polynomial temp(0, {knots[i].second});
            polynomials[i] = polynomials[i] + temp * basis_polynomials[i][n][0];
            polynomials.push_back(temp * basis_polynomials[i][n][1]);
        }
        for (int i = 0; i <= n; ++i)
            polynomials.erase(polynomials.begin());
        for (int i = 0; i < n; ++i)
            knots.erase(knots.begin());
        return Spline(N, knots, polynomials);
    }

    Spline complete_cubic_Bspline() {
        double h1 = knots[1].first - knots[0].first;
        double h2 = knots[N - 1].first - knots[N - 2].first;
        for (int i = 1; i <= n; ++i) {
            knots.insert(knots.begin(), std::make_pair(knots[0].first - h1, 0.0));
            knots.push_back(std::make_pair(knots[N - 1 + i].first + i * h2, 0.0));
        }

        for (int i = 1; i <= N + n - 1; ++i)
            construct_basis_polynomials(i, n);

        double coeff_matrix[N * N];

        coeff_matrix[0] = 1.0/3 + (knots[4].first - knots[2].first) / (2 * (knots[4].first + knots[3].first - 2 * knots[2].first));
        coeff_matrix[N] = (knots[3].first - knots[2].first) / (knots[4].first + knots[3].first - 2 * knots[1].first);

        for (int i = 1; i < N - 1; ++i) {
            coeff_matrix[(i - 1)*N + i] = basis_polynomials[i + 1][n][2](knots[i + 2].first);
            coeff_matrix[i*N + i] = basis_polynomials[i + 2][n][1](knots[i + 2].first);
            coeff_matrix[i*N + i + 1] = basis_polynomials[i + 3][n][0](knots[i + 2].first);
        }

        coeff_matrix[N*N - N - 1] = (knots[N + 1].first - knots[N].first) / (2 * knots[N + 2].first - knots[N + 2].first - knots[N - 1].first);
        coeff_matrix[N*N - 1] = 1.0/3 + (knots[N + 1].first - knots[N - 1].first) / (2 * (2 * knots[N + 1].first - knots[N].first - knots[N - 1].first));

        double rhs[N];

        rhs[0] = knots[0].second + boundary_conditions[0] * (knots[3].first - knots[2].first) / 3;
        for (int i = 1; i < N - 1; ++i)
            rhs[i] = knots[i + 2].second;
        rhs[N - 1] = knots[N + 1].second - boundary_conditions[1] * (knots[N + 1].first - knots[N].first) / 3;

        EquationSolver linear_system(N, coeff_matrix, rhs);
        double *solution = linear_system.solve();

        std::vector<Polynomial> polynomials;
        polynomials.push_back(Polynomial(0, {0.0}));
        for (int i = 1; i <= N + 2; ++i) {
            if (i == 1) {
                double co = 3 * (knots[4].first - knots[2].first) * solution[1] / (knots[5].first + knots[4].first - 2 * knots[3].first) - 2 * (knots[4].first - knots[3].first) * boundary_conditions[0];
                Polynomial temp(0, {co});
                polynomials.push_back(temp * basis_polynomials[i][n][0]);
                polynomials.push_back(temp * basis_polynomials[i][n][1]);
                polynomials.push_back(temp * basis_polynomials[i][n][2]);
                continue;
            }
            if (i == N + 2) {
                double co = 3 * (knots[N + 2].first - knots[N + 1].first) * solution[N - 2] / (2 * knots[N + 2].first - knots[N + 1].first - knots[N].first) + 2 * (knots[N + 2].first - knots[N + 1].first) * boundary_conditions[1];
                Polynomial temp(0, {co});
                polynomials[N + 2] = polynomials[N + 2] + temp * basis_polynomials[i][n][0];
                continue;
            }
            Polynomial temp(0, {solution[i - 2]});
            polynomials[i] = polynomials[i] + temp * basis_polynomials[i][n][0];
            polynomials[i + 1] = polynomials[i + 1] + temp * basis_polynomials[i][n][1];
            polynomials.push_back(temp * basis_polynomials[i][n][2]); 
        }
        for (int i = 0; i <= n; ++i)
            polynomials.erase(polynomials.begin());
        for (int i = 0; i < n; ++i)
            knots.erase(knots.begin());
        return Spline(N, knots, polynomials);
    }

    Spline natural_cubic_Bspline() {
        double h1 = knots[1].first - knots[0].first;
        double h2 = knots[N - 1].first - knots[N - 2].first;
        for (int i = 1; i <= n; ++i) {
            knots.insert(knots.begin(), std::make_pair(knots[0].first - h1, 0.0));
            knots.push_back(std::make_pair(knots[N - 1 + i].first + i * h2, 0.0));
        }

        for (int i = 1; i <= N + n - 1; ++i)
            construct_basis_polynomials(i, n);

        double coeff_matrix[N * N];

        coeff_matrix[0] = 1.0/3 + (knots[4].first - knots[2].first) / (2 * (knots[4].first + knots[3].first - 2 * knots[2].first));

        for (int i = 1; i < N - 1; ++i) {
            coeff_matrix[(i - 1)*N + i] = basis_polynomials[i + 1][n][2](knots[i + 2].first);
            coeff_matrix[i*N + i] = basis_polynomials[i + 2][n][1](knots[i + 2].first);
            coeff_matrix[i*N + i + 1] = basis_polynomials[i + 3][n][0](knots[i + 2].first);
        }

        coeff_matrix[N*N - 1] = 1.0/3 + (knots[N + 1].first - knots[N - 1].first) / (2 * (2 * knots[N + 1].first - knots[N].first - knots[N - 1].first));

        double rhs[N];

        rhs[0] = knots[0].second + boundary_conditions[0] * (knots[3].first - knots[2].first) / 3;
        for (int i = 1; i < N - 1; ++i)
            rhs[i] = knots[i + 2].second;
        rhs[N - 1] = knots[N + 1].second - boundary_conditions[1] * (knots[N + 1].first - knots[N].first) / 3;

        EquationSolver linear_system(N, coeff_matrix, rhs);
        double *solution = linear_system.solve();

        std::vector<Polynomial> polynomials;
        polynomials.push_back(Polynomial(0, {0.0}));
        for (int i = 1; i <= N + 2; ++i) {
            if (i == 1) {
                double co = -3 * (knots[4].first - knots[2].first) * solution[1] / (knots[5].first + knots[4].first - 2 * knots[3].first);
                Polynomial temp(0, {co});
                polynomials.push_back(temp * basis_polynomials[i][n][0]);
                polynomials.push_back(temp * basis_polynomials[i][n][1]);
                polynomials.push_back(temp * basis_polynomials[i][n][2]);
                continue;
            }
            if (i == N + 2) {
                double co = -3 * (knots[N + 2].first - knots[N + 1].first) * solution[N - 2] / (2 * knots[N + 2].first - knots[N + 1].first - knots[N].first);
                Polynomial temp(0, {co});
                polynomials[N + 2] = polynomials[N + 2] + temp * basis_polynomials[i][n][0];
                continue;
            }
            Polynomial temp(0, {solution[i - 2]});
            polynomials[i] = polynomials[i] + temp * basis_polynomials[i][n][0];
            polynomials[i + 1] = polynomials[i + 1] + temp * basis_polynomials[i][n][1];
            polynomials.push_back(temp * basis_polynomials[i][n][2]); 
        }
        for (int i = 0; i <= n; ++i)
            polynomials.erase(polynomials.begin());
        for (int i = 0; i < n; ++i)
            knots.erase(knots.begin());
        return Spline(N, knots, polynomials);
    }

    Spline period_cubic_Bspline() {

        knots.insert(knots.begin(), std::make_pair(knots[0].first + knots[N - 2].first - knots[N - 1].first, knots[N - 2].second));
        knots.insert(knots.begin(), std::make_pair(knots[1].first + knots[N - 2].first - knots[N].first, knots[N - 2].second));
        knots.insert(knots.begin(), std::make_pair(knots[2].first + knots[N - 2].first - knots[N + 1].first, knots[N - 2].second));
        knots.push_back(std::make_pair(knots[N + 2].first + knots[4].first - knots[3].first, knots[4].second));
        knots.push_back(std::make_pair(knots[N + 2].first + knots[5].first - knots[3].first, knots[5].second));
        knots.push_back(std::make_pair(knots[N + 2].first + knots[6].first - knots[3].first, knots[6].second));

        for (int i = 1; i <= N + n - 1; ++i)
            construct_basis_polynomials(i, n);

        double coeff_matrix[N * N];

        coeff_matrix[0] = basis_polynomials[2][n][1](knots[3].first);
        coeff_matrix[N] = basis_polynomials[3][n][0](knots[3].first);
        coeff_matrix[N * N - 2 * N] = basis_polynomials[1][n][2](knots[3].first);

        for (int i = 1; i < N - 1; ++i) {
            coeff_matrix[(i - 1)*N + i] = basis_polynomials[i + 1][n][2](knots[i + 2].first);
            coeff_matrix[i*N + i] = basis_polynomials[i + 2][n][1](knots[i + 2].first);
            coeff_matrix[i*N + i + 1] = basis_polynomials[i + 3][n][0](knots[i + 2].first);
        }

        coeff_matrix[2 * N - 1] = basis_polynomials[N + 2][n][0](knots[N + 2].first);
        coeff_matrix[N * N - N - 1] = basis_polynomials[N][n][2](knots[N + 2].first);
        coeff_matrix[N * N - 1] = basis_polynomials[N + 1][n][1](knots[N + 2].first);

        double rhs[N];

        rhs[0] = knots[0].second + boundary_conditions[0] * (knots[3].first - knots[2].first) / 3;
        for (int i = 1; i < N - 1; ++i)
            rhs[i] = knots[i + 2].second;
        rhs[N - 1] = knots[N + 1].second - boundary_conditions[1] * (knots[N + 1].first - knots[N].first) / 3;

        EquationSolver linear_system(N, coeff_matrix, rhs);
        double *solution = linear_system.solve();

        std::vector<Polynomial> polynomials;
        polynomials.push_back(Polynomial(0, {0.0}));
        for (int i = 1; i <= N + 2; ++i) {
            if (i == 1) {
                double co = solution[N - 2];
                Polynomial temp(0, {co});
                polynomials.push_back(temp * basis_polynomials[i][n][0]);
                polynomials.push_back(temp * basis_polynomials[i][n][1]);
                polynomials.push_back(temp * basis_polynomials[i][n][2]);
                continue;
            }
            if (i == N + 2) {
                double co = solution[1];
                Polynomial temp(0, {co});
                polynomials[N + 2] = polynomials[N + 2] + temp * basis_polynomials[i][n][0];
                continue;
            }
            Polynomial temp(0, {solution[i - 2]});
            polynomials[i] = polynomials[i] + temp * basis_polynomials[i][n][0];
            polynomials[i + 1] = polynomials[i + 1] + temp * basis_polynomials[i][n][1];
            polynomials.push_back(temp * basis_polynomials[i][n][2]); 
        }
        for (int i = 0; i <= n; ++i)
            polynomials.erase(polynomials.begin());
        for (int i = 0; i < n; ++i)
            knots.erase(knots.begin());
        return Spline(N, knots, polynomials);
    }

    Spline arbitrary_order_Bspline(double *a) {
        double h1 = knots[1].first - knots[0].first;
        double h2 = knots[N - 1].first - knots[N - 2].first;
        for (int i = 1; i <= n; ++i) {
            knots.insert(knots.begin(), std::make_pair(knots[0].first - h1, 0.0));
            knots.push_back(std::make_pair(knots[N - 1 + i].first + i * h2, 0.0));
        }

        for (int i = 1; i <= N + n - 1; ++i)
            construct_basis_polynomials(i, n);
        
        std::vector<Polynomial> polynomials;
        polynomials.push_back(Polynomial(0, {0.0}));
        for (int i = 1; i <= N + n - 1; ++i) {
            if (i == 1) {
                Polynomial temp(0, {a[0]});
                polynomials.push_back(temp * basis_polynomials[i][n][0]);
                polynomials.push_back(temp * basis_polynomials[i][n][1]);
                polynomials.push_back(temp * basis_polynomials[i][n][2]);
                continue;
            }
            if (i == N + n - 1) {
                Polynomial temp(0, {a[N + n - 2]});
                polynomials[N + n - 1] = polynomials[N + n - 1] + temp * basis_polynomials[i][n][0];
                continue;
            }
            Polynomial temp(0, {a[i - 1]});
            polynomials[i] = polynomials[i] + temp * basis_polynomials[i][n][0];
            polynomials[i + 1] = polynomials[i + 1] + temp * basis_polynomials[i][n][1];
            polynomials.push_back(temp * basis_polynomials[i][n][2]); 
        }
        for (int i = 0; i <= n; ++i)
            polynomials.erase(polynomials.begin());
        for (int i = 0; i < n; ++i)
            knots.erase(knots.begin());
        return Spline(N, knots, polynomials);
    }

private:
    int n, N; //S_n^(n - 1) spline
    double boundary_conditions[5]; //boundary conditions, the number of elements equals to k
    std::vector<std::pair<double, double>> knots; //t_0 to t_{N - 1}
    std::vector<std::vector<std::vector<Polynomial>>> basis_polynomials; //basis polynomials
    std::vector<double> coefficients; //coefficients of the spline, a_i
};

#endif