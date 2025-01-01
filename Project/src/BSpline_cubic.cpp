#include <cstdio>
#include <iostream>
#include "Bspline.h"

double func(double x) {
    return 1.0 / (1.0 + x * x);
}

double der_func(double x){
    return -2.0 * x * func(x) * func(x);
}

int main() {
    double boundary[2];
    boundary[0] = der_func(-5.0);
    boundary[1] = der_func(5.0);
    printf("2\n");
    std::vector<std::pair<double, double>> knots;
    for (int i = 0; i <= 10; ++i) {
        knots.push_back(std::make_pair(-5.0 + i, func(-5.0 + i)));
    }
    Bspline form1(3, 11, knots, boundary);
    Spline sp1 = form1.complete_cubic_Bspline();
    printf("Cubic B-spline\n");
    sp1.generate_data();

    knots.clear();
    for (int i = 1; i <= 10; ++i) 
        knots.push_back(std::make_pair(-5.5 + i, func(-5.5 + i)));
}