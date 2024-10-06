#include "EquationSolver.h"
#include <iostream>
#include <cmath>
#include <cstdio>

double V = 12.4, L = 10.0, r = 1.0;

int main() {
    freopen("../result/output.txt", "a", stdout);
    puts("\nProblem E:");
    auto func = [](double h) {
        return V - L*(0.5*pi*r*r - r*r*std::asin(h/r) - h*std::sqrt(r*r - h*h));
    };
    auto de = [](double h) {
        return 2*L*std::sqrt(r*r - h*h);
    };
    
    Bisection_Method<double> eq1(func);
    puts("Result by bisection method:");
    eq1.solve(1e-4, 1e-2, 0, 1);
    
    Newton_Method<double> eq2(func, de);
    puts("Result by Newton's method:");
    eq2.solve(1e-4, 1e-2, 0.1);    

    Secant_Method<double> eq3(func);
    puts("Result by secant method:");
    eq3.solve(1e-4, 1e-2, 0, 1);

    fclose(stdout);
    return 0;
}