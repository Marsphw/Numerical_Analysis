#ifndef PPFORM_H
#define PPFORM_H

#include <iostream>
#include <cstring>
#include "Polynomial.h"
#include "Solving_Equation.h"
#include "Spline.h"

class PPForm {

public:
    PPForm(int nn, int kk, int NN, std::vector<std::pair<double, double>> knots_in, double *bc = nullptr): n(nn), k(kk), N(NN), knots(knots_in), boundary_conditions(bc) {}

    Spline linear_splines() {
        std::vector<Polynomial> S;
        for (int i = 0; i < N - 1; ++i) {
            double a = (knots[i + 1].second - knots[i].second) / (knots[i + 1].first - knots[i].first);
            double b = knots[i].second - a * knots[i].first;
            std::vector<double> coeffs = {b, a};
            Polynomial P(1, coeffs);
            S.push_back(P);
        }
        Spline result(N, knots, S);
        return result; 
    }

    double get_lambda(int i) {
        return (knots[i + 1].first - knots[i].first) / (knots[i + 1].first - knots[i - 1].first);
    }

    double get_mu(int i) {
        return (knots[i].first - knots[i - 1].first) / (knots[i + 1].first - knots[i - 1].first);
    }

    double get_K(int i) {
        return (knots[i + 1].second - knots[i].second) / (knots[i + 1].first - knots[i].first);
    } 

    std::vector<Polynomial> construct_polynomials(double *solution) {
        std::vector<Polynomial> S;
        for (int i = 0; i < N - 1; ++i) {
            double c[4];
            c[0] = knots[i].second;
            c[1] = solution[i];
            c[2] = (3 * get_K(i) - 2 * solution[i] - solution[i + 1]) / (knots[i + 1].first - knots[i].first);
            c[3] = (solution[i] + solution[i + 1] - 2 * get_K(i)) / ((knots[i + 1].first - knots[i].first) * (knots[i + 1].first - knots[i].first));
            std::vector<double> coeffs(4);
            coeffs[3] = c[3];
            coeffs[2] = c[2] - 3 * c[3] * knots[i].first;
            coeffs[1] = c[1] - 2 * c[2] * knots[i].first + 3 * c[3] * knots[i].first * knots[i].first;
            coeffs[0] = c[0] - c[1] * knots[i].first + c[2] * knots[i].first * knots[i].first - c[3] * knots[i].first * knots[i].first * knots[i].first;
            Polynomial P(3, coeffs);
            S.push_back(P);
        }
        return S;
    }

    Spline complete_cubic_ppform() {
        double coeff_matrix[N * N];
        memset(coeff_matrix, 0.0, sizeof(double) * N * N);
        coeff_matrix[0] = 1;
        coeff_matrix[N * N - 1] = 1;
        for (int i = 1; i < N - 1; ++i) {
            coeff_matrix[(i - 1)*N + i] = get_lambda(i);
            coeff_matrix[i*N + i] = 2;
            coeff_matrix[(i + 1)*N + i] = get_mu(i);
        }
        double rhs[N];
        rhs[0] = boundary_conditions[0];
        for (int i = 1; i < N - 1; ++i) 
            rhs[i] = 3 * (get_mu(i) * get_K(i) + get_lambda(i) * get_K(i - 1));
        rhs[N - 1] = boundary_conditions[1];

        EquationSolver eq(N, coeff_matrix, rhs);
        double *solution = eq.solve();

        return Spline(N, knots, construct_polynomials(solution));
    }

    Spline natural_cubic_ppform() {
        boundary_conditions[0] = 0;
        boundary_conditions[1] = 0;
        return complete_cubic_ppform();
    }

    Spline periodic_cubic_ppform() {
        double coeff_matrix[N * N] = {};
        coeff_matrix[0] = 1;
        coeff_matrix[N * N - N] = -1;
        for (int i = 1; i < N - 1; ++i) {
            coeff_matrix[(i - 1)*N + i] = get_lambda(i);
            coeff_matrix[i*N + i] = 2;
            coeff_matrix[(i + 1)*N + i] = get_mu(i);
        }
        coeff_matrix[N - 1] = (knots[N - 1].first - knots[N - 2].first) / (knots[N - 1].first - knots[N - 2].first + knots[1].first - knots[0].first);
        coeff_matrix[N * N - N - 1] = (knots[1].first - knots[0].first) / (knots[N - 1].first - knots[N - 2].first + knots[1].first - knots[0].first);
        coeff_matrix[N * N - 1] = 2;

        double rhs[N];
        rhs[0] = 0;
        for (int i = 1; i < N - 1; ++i) 
            rhs[i] = 3 * (get_mu(i) * get_K(i) + get_lambda(i) * get_K(i - 1));
        rhs[N - 1] = 3 * (get_mu(N - 1) * get_K(0) + get_lambda(N - 1) * (knots[0].second - knots[N - 2].second) / (knots[N - 1].first - knots[N - 2].first));

        EquationSolver eq(N, coeff_matrix, rhs);
        double *solution = eq.solve();

        return Spline(N, knots, construct_polynomials(solution));
    }
    
private:
    int n, k, N; //S_n^k spline
    double* boundary_conditions; //boundary conditions, the number of elements equals to k
    std::vector<std::pair<double, double>> knots; //knot vector
};

#endif