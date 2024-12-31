#include <iostream>
#include <cstdio>
#include "ppForm.h"

int main() {
    int N;
    scanf("%d", &N);
    std::vector<std::pair<double, double>> knots(N + 1);
    for (int i = 0; i < N; i++) 
        scanf("%lf %lf", &knots[i].first, &knots[i].second);
    PPForm form(1, 0, N, knots);
    Spline sp = form.linear_splines();
    printf("1\n");
    printf("Linear ppForm\n");
    sp.generate_data();
    return 0;
}