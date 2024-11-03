#include "Interpolation.h"
#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>

double f(double x) {
    return 1 / (1 + 25 * x * x);
}

double Chebyshev(double x, int n) {
    return std::cos(n * std::acos(x));
}

int main() {
    int n[5] = {5, 10, 15, 20};
    for (int i = 0; i < 4; ++i) {
        std::vector<Point> points;
        for (int k = 1; k <= n[i]; ++k) {
            double x = std::cos((2*k - 1)*pi/(2*n[i]));
            double y = f(x);
            Point temp;
            temp.x = x;
            temp.order = 0;
            temp.f = new double[1];
            temp.f[0] = y;
            points.push_back(temp);
        }
        Newton_Interpolation interpolation(points);
        Polynomial p = interpolation.interpolate();
        std::cout << n[i] << std::endl;
        for (int j = 0; j < 1000; ++j) {
            double x = -1 + j/500.0;
            double y = p(x);
            std::cout << y << std::endl;
        }
    }
}