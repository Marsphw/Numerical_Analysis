#include <cstdio>
#include <iostream>
#include "ppForm.h"

double func(double x){
    return 1.0 / (1 + 25 * x * x);
}

double der_func(double x){
    return -50 * x * func(x) * func(x);
}

int main() {
    int N[5] = {6, 11, 21, 41, 81};
    std::vector<double> boundary;
    boundary.push_back(der_func(-1.0));
    boundary.push_back(der_func(1.0));
    printf("5\n");
    for (int i = 0; i < 5; ++i) {
        double h = 2.0 / (N[i] - 1);
        std::vector<std::pair<double, double>> knots(N[i]);
        for (int j = 0; j < N[i]; ++j) 
            knots[j] = std::make_pair(-1.0 + j * h, func(-1.0 + j * h));
        PPForm form(3, 2, N[i], knots, boundary);
        Spline sp = form.complete_cubic_ppform();
        printf("Cubic ppForm with %d knots\n", N[i]);
        sp.generate_data();
        //sp.show_polynomials();
    }
    return 0;
}