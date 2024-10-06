#include "EquationSolver.h"
#include <iostream>
#include <cmath>
#include <cstdio>

int main() {
    freopen("../result/output.txt", "a", stdout);
    puts("\nProblem D:");

    auto func1 = [](double x) {return sin(x / 2) - 1.0;};
    Secant_Method<double> eq1(func1);
    puts("Question1.1: sin(x/2) - 1 with x0 = 0 and x1 = pi/2");
    eq1.solve(1e-6, 1e-4, 0, pi / 2);
    puts("Question1.2: sin(x/2) - 1 with x0 = 7pi/2 and x1 = 4pi");
    eq1.solve(1e-6, 1e-4, 7*pi/2, 4*pi);

    auto func2 = [](double x) {return exp(x) - tan(x);};
    Secant_Method<double> eq2(func2);
    puts("Question2.1: e^x - tan(x) with x0 = 1 and x1 = 1.4");
    eq2.solve(1e-6, 1e-4, 1, 1.4);
    puts("Question2.2: e^x - tan(x) with x0 = -3 and x1 = -3.4");
    eq2.solve(1e-6, 1e-4, -3, -3.4);

    auto func3 = [](double x) {return std::pow(x, 3) - 12*std::pow(x, 2) + 3*x + 1.0;};
    Secant_Method<double> eq3(func3);
    puts("Question3.1: x^3 - 12x^2 + 3x + 1 with x0 = 0 and x1 = -0.5");
    eq3.solve(1e-6, 1e-4, 0, -0.5);
    puts("Question3.2: x^3 - 12x^2 + 3x + 1 with x0 = 0 and x1 = 0.5");
    eq3.solve(1e-6, 1e-4, 0, 0.5);

    fclose(stdout);
    return 0;
}