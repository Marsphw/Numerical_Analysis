#include "EquationSolver.h"
#include <iostream>
#include <cmath>
#include <cstdio>

int main() {
    freopen("../result/output.txt", "w", stdout);
    puts("Problem B:");

    //Q1
    auto func1 = [](double x) {
        return 1/x - tan(x);
    };
    Bisection_Method<double> eq1(func1);
    puts("Question1: 1/x - tan(x) on [0, pi/2]");
    eq1.solve(1e-6, 1e-4, 0.1, pi / 2);

    //Q2
    auto func2 = [](double x) {
        return 1/x - std::pow(2, x);
    };
    Bisection_Method<double> eq2(func2);
    puts("Question2: 1/x - 2^x on [0, 1]");
    eq2.solve(1e-6, 1e-4, 0.1, 1);

    //Q3
    auto func3 = [](double x) {
        return std::pow(2, -x) + std::exp(x) + 2*cos(x) - 6.0;
    };
    Bisection_Method<double> eq3(func3);
    puts("Question3: 2^(-x) + e^x + 2cos(x) - 6 on [1, 3]");
    eq3.solve(1e-6, 1e-4, 1, 3);   

    //Q4
    auto func4 = [](double x) {
        return (std::pow(x, 3) + 4*std::pow(x, 2) + 3*x + 5.0)/(2*std::pow(x, 3) - 9*std::pow(x, 2) + 18*x - 2.0);
    };
    Bisection_Method<double> eq4(func4);
    puts("Question4: (x^3 + 4x^2 + 3x^2 + 5)/(2x^3 - 9x^2 + 18x - 2) on [0, 4]");
    eq4.solve(1e-6, 1e-4, 1e-6, 4);

    fclose(stdout);
    return 0;
}