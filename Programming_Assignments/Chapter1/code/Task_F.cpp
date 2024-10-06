#include "EquationSolver.h"
#include <iostream>
#include <cmath>
#include <cstdio>

double l = 89.0, h = 49.0, A, B, C, D, E, beta = 11.5*pi/180.0;


int main() {
    freopen("../result/output.txt", "a", stdout);
    puts("\nProblem F:");

    A = l*std::sin(beta);
    B = l*std::cos(beta);

    auto func = [](double alpha) {
        return A*std::sin(2*alpha)/2 + B*std::sin(alpha)*std::sin(alpha) - C*std::cos(alpha) - E*std::sin(alpha);
    };
    auto de = [](double alpha) {
        return A*std::cos(2*alpha) + B*std::sin(2*alpha) + C*std::sin(alpha) - E*std::cos(alpha);
    };

    //Task a
    D = 55.0;
    C = (h + 0.5*D)*std::sin(beta) - 0.5*D*std::tan(beta);
    E = (h + 0.5*D)*std::cos(beta) - 0.5*D;
    Newton_Method<double> eq1(func, de);
    puts("Task a:");
    eq1.solve(1e-6, 1e-4, 33*pi/180.0);

    //Task b
    D = 30.0;
    C = (h + 0.5*D)*std::sin(beta) - 0.5*D*std::tan(beta);
    E = (h + 0.5*D)*std::cos(beta) - 0.5*D;
    Newton_Method<double> eq2(func, de);
    puts("Task b:");
    eq2.solve(1e-6, 1e-4, 33*pi/180.0);

    //Task c
    D = 55.0;
    C = (h + 0.5*D)*std::sin(beta) - 0.5*D*std::tan(beta);
    E = (h + 0.5*D)*std::cos(beta) - 0.5*D;
    Secant_Method<double> eq3(func);
    puts("Task c:");
    printf("The another initial guess is 93 degrees: ");
    eq3.solve(1e-6, 1e-4, 93*pi/180.0, 33*pi/180.0);
    printf("The another initial guess is 153 degrees: ");
    eq3.solve(1e-6, 1e-4, 153*pi/180.0, 33*pi/180.0);
    printf("The another initial guess is 213 degrees: ");
    eq3.solve(1e-6, 1e-4, 213*pi/180.0, 33*pi/180.0);

    fclose(stdout);
    return 0;
}