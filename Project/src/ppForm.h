#ifndef PPFORM_H
#define PPFORM_H

#include <iostream>
#include "Polynomial.h"
#include "Equation_Solving.h"
#include "Spline.h"

class PPForm {

public:
    PPForm(int nn, int kk, int NN, std::vector<std::pair<double, double>> knots_in, double bc[] = 0): n(nn), k(kk), N(NN), knots(knots_in) {
        for (int i = 0; i < k; ++i) 
            boundary_conditions[i] = bc[i];
    }

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

    std::vector<Polynomial> construct_polynomials(std::vector<double> solution) {
        std::vector<Polynomial> S;
        for (int i = 0; i < N - 1; ++i) {
            double c[4];
            c[0] = knots[i].second;
            c[1] = solution[i];
            c[2] = (3 * get_K(i) - 2 * solution[i] - 2 * solution[i + 1]) / (knots[i + 1].first - knots[i].first);
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
        Matrix<double> coeff_matrix(N, N);
        coeff_matrix(0, 0) = 1;
        for (int i = 1; i < N - 1; ++i) {
            coeff_matrix(i, i - 1) = get_lambda(i);
            coeff_matrix(i, i) = 2;
            coeff_matrix(i, i + 1) = get_mu(i);
        }
        coeff_matrix(N - 1, N - 1) = 1;

        std::vector<double> rhs(N);
        rhs[0] = boundary_conditions[0];
        for (int i = 1; i < N - 1; ++i) 
            rhs[i] = 3 * (get_mu(i) * get_K(i) + get_lambda(i) * get_K(i - 1));
        rhs[N - 1] = boundary_conditions[1];

        Equation linear_system(N, coeff_matrix, rhs);
        std::vector<double> solution = linear_system.solve();

        return Spline(N, knots, construct_polynomials(solution));
    }

    Spline natural_cubic_ppform() {
        boundary_conditions[0] = 0;
        boundary_conditions[1] = 0;
        return complete_cubic_ppform();
    }

    Spline periodic_cubic_ppform() {
        Matrix<double> coeff_matrix(N, N);
        coeff_matrix(0, 0) = 1;
        coeff_matrix(0, N - 1) = -1;
        for (int i = 1; i < N - 1; ++i) {
            coeff_matrix(i, i - 1) = get_lambda(i);
            coeff_matrix(i, i) = 2;
            coeff_matrix(i, i + 1) = get_mu(i);
        }
        coeff_matrix(N - 1, 0) = (knots[N - 1].first - knots[N - 2].first) / (knots[N - 1].first - knots[N - 2].first + knots[1].first - knots[0].first);
        coeff_matrix(N - 1, N - 2) = (knots[1].first - knots[0].first) / (knots[N - 1].first - knots[N - 2].first + knots[1].first - knots[0].first);
        coeff_matrix(N - 1, N - 1) = 2;

        std::vector<double> rhs(N);
        rhs[0] = 0;
        for (int i = 1; i < N - 1; ++i) 
            rhs[i] = 3 * (get_mu(i) * get_K(i) + get_lambda(i) * get_K(i - 1));
        rhs[N - 1] = 3 * (get_mu(N - 1) * get_K(0) + get_lambda(N - 1) * (knots[0].second - knots[N - 2].second) / (knots[N - 1].first - knots[N - 2].first));

        Equation linear_system(N, coeff_matrix, rhs);
        std::vector<double> solution = linear_system.solve();

        return Spline(N, knots, construct_polynomials(solution));
    }
    
private:
    int n, k, N; //S_n^k spline
    double boundary_conditions[5]; //boundary conditions, the number of elements equals to k
    std::vector<std::pair<double, double>> knots; //knot vector
};

#endif