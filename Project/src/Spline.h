#ifndef SPLINE_H
#define SPLINE_H

#include <vector>
#include <iostream>
#include <cstdio>
#include <utility>
#include "Polynomial.h"

class Spline {

public:
    Spline(int NN, std::vector<std::pair<double, double>> knots_in, std::vector<Polynomial> s_in): N(NN), knots(knots_in), s(s_in) {}

    double operator()(double x) const {
        for (int i = 0; i < N; ++i) {
            if (x >= knots[i].first && x <= knots[i + 1].first)
                return s[i](x);
        }
        return 0.0;
    }

    void generate_data() {
        printf("%d\n", N);
        for (int i = 0; i < N; ++i) 
            printf("%.6f\n%.6f\n", knots[i].first, knots[i].second);
        double span = knots[N - 1].first - knots[0].first;
        double step = span / 1000;
        double x = knots[0].first;
        for (int i = 0; i < 1000; ++i) 
            printf("%.6f\n", (*this)(x + i * step));
    }

    void show_polynomials() {
        for (int i = 0; i < N - 1; ++i) 
            s[i].show();
    }
    
private:
    int N; //number of knots
    std::vector<std::pair<double, double>> knots; //knot vector
    std::vector<Polynomial> s;
};

#endif 