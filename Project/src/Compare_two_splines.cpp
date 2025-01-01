#include <iostream>
#include <cstdio>
#include "Bspline.h"
#include "ppForm.h"

int main() {
    int N;
    scanf("%d", &N);
    std::vector<std::pair<double, double>> knots(N + 1);
    for (int i = 0; i < N; ++i) 
        scanf("%lf %lf", &knots[i].first, &knots[i].second);
    PPForm form1(3, 2, N, knots);
    Spline sp1 = form1.natural_cubic_ppform(); 
    printf("2\n");
    printf("Cubic Natural ppform-Spline\n");
    sp1.generate_data();
    Bspline form2(3, N, knots);
    Spline sp2 = form2.natural_cubic_Bspline();
    printf("Cubic Natural B-Spline\n");
    sp2.generate_data();
    return 0;
}
