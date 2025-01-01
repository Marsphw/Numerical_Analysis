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
    std::vector<double> boundary;
    boundary.clear();
    boundary.push_back(der_func(-5.0));
    boundary.push_back(der_func(5.0));
    printf("2\n");
    std::vector<std::pair<double, double>> knots;
    for (int i = 0; i <= 10; ++i) 
        knots.push_back(std::make_pair(-5.0 + i, func(-5.0 + i)));
    Bspline form1(3, 11, knots, boundary);
    Spline sp1 = form1.complete_cubic_Bspline();
    printf("Cubic B-spline\n");
    sp1.generate_data();
    //sp1.show_polynomials();

    knots.clear();
    boundary.clear();
    for (int i = 1; i <= 10; ++i) 
        knots.push_back(std::make_pair(-5.5 + i, func(-5.5 + i)));
    boundary.push_back(func(-5.0));
    boundary.push_back(func(5.0));
    Bspline form2(2, 10, knots, boundary);
    Spline sp2 = form2.cardinal_quad_Bspline();
    printf("Cardinal quadratic B-spline\n");
    sp2.generate_data();
    //sp2.show_polynomials();

    return 0;
}