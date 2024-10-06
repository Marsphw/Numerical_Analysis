#include "EquationSolver.h"
#include <iostream>
#include <cmath>
#include <cstdio>

int main() {
    freopen("../result/output.txt", "a", stdout);
    puts("\nProblem C:");
    auto func = [](double x) {return x - tan(x);};
    auto de = [](double x) {return 1 - 1/std::pow(cos(x), 2);};
    Newton_Method<double> eq(func, de);
    eq.solve(1e-6, 1e-4, 4.5);
    eq.solve(1e-6, 1e-4, 7.7);
    fclose(stdout);
    return 0;
}