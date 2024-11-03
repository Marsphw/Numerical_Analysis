#include "Interpolation.h"
#include <iostream>
#include <cstdio>
#include <vector>

double f(double x) {
    return 1 / (1 + x * x);
}

int main() {
    int n[5] = {2, 4, 6, 8};
    for (int i = 0; i < 4; ++i) {
        std::vector<Point> points = {};
        for (int j = 0; j <= n[i]; ++j) {
            double x = -5 + double(10*j)/n[i];
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
            double x = -5 + j/100.0;
            double y = p(x);
            std::cout << y << std::endl;
        }
    }
    return 0;
}