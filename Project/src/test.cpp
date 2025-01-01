#include "Solving_Equation.h"
#include <iostream>

int main() {
    int N = 4;
    double A[16] = {  1,  2,  3,  1,
                      4,  2,  0,  2,
                     -2,  0, -1,  2,
                      3,  4,  2, -3};
    double B[4] = {  6,  2,  1,  8};
    EquationSolver solver(2, A, B);
    double* solution = solver.solve();
    for (int i = 0; i < 4; i++) {
        std::cout << "x" << i+1 << " = " << solution[i] << std::endl;
    }
    return 0;
}