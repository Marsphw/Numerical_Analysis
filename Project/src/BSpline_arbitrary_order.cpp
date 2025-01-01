#include <cstdio>
#include <iostream>
#include "Bspline.h"

int main() {
    int N;
    std::cin >> N;
    std::vector<std::pair<double, double>> knots(N + 1);
    for (int i = 0; i < N; ++i) {
        scanf("%lf", &knots[i].first);
        knots[i].second = 0;
    }
    int n;
    std::cin >> n;
    Bspline form(n, N, knots);
    double a[n + N + 1];
    for (int i = 0; i < n + N - 1; ++i) 
        std::cin >> a[i];
    Spline sp = form.arbitrary_order_Bspline(a);
    std::cout << "1" << std::endl;
    printf("B-spline of order %d\n", n);
    sp.generate_data();
    return 0;
}