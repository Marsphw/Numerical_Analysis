#ifndef SPLINE_H
#define SPLINE_H

#include <vector>
#include <iostream>
#include <cstdio>
#include "Polynomial.h"

class Spline {

public:
    Spline(int NN, std::vector<double> knots_in, std::vector<Polynomial> s_in): N(NN), knots(knots_in), s(s_in) {}
    double operator()(double x) const {
        for (int i = 0; i < N; ++i) {
            if (x >= knots[i] && x <= knots[i + 1])
                return s[i](x);
        }
    }
    void generate_data() {
        double span = knots[N - 1] - knots[0];
        double step = span / 1000;
        for (double x = knots[0]; x <= knots[N - 1]; x += step) 
            printf("%.6f %.6f\n", x, (*this)(x));
    }
    
private:
    int N; //number of knots
    std::vector<double> knots; //knot vector
    std::vector<Polynomial> s;
};

#endif 