#include <iostream>
#include <cstdio>
#include "Bspline.h"

int main() {
    int N;
    scanf("%d", &N);
    std::vector<std::pair<double, double>> knots(N + 1);
    for (int i = 0; i < N; i++) 
        scanf("%lf %lf", &knots[i].first, &knots[i].second);
    Bspline form(1, N, knots);
    Spline sp = form.linear_Bspline();
    printf("1\n");
    printf("Linear B-spline:\n");
    sp.generate_data();
    //sp.show_polynomials();
    return 0;
}