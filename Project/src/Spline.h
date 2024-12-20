#ifndef SPLINE_H
#define SPLINE_H

#include <vector>
#include <iostream>
#include <cstdio>
#include "Polynomial.h"

struct knot {
    int order;
    double value[5];
};

class Spline {

public:
    Spline(int nn, std::vector<knot> knots_in, std::vector<Polynomial> s_in): N(nn), knots(knots_in), s(s_in) {}
    double operator()(double x) const {
        for (int i = 0; i < N; ++i) {
            if (x >= knots[i].value[0] && x <= knots[i + 1].value[0])
                return s[i](x);
        }
    }
    
private:
    int N; //number of knots
    std::vector<knot> knots; //knot vector
    std::vector<Polynomial> s;
};

#endif 