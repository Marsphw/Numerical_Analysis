#include "Interpolation.h"
#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>

int main() {
    std::vector<Point> points;
    for (int i = 0; i < 5; ++i) {
        double x, f, f1;
        std::cin >> x >> f >> f1;
        Point temp;
        temp.x = x;
        temp.order = 1;
        temp.f = new double[2];
        temp.f[0] = f;
        temp.f[1] = f1;
        points.push_back(temp);
    }

    Hermite_Interpolation interpolation(points);
    Polynomial p = interpolation.interpolate();
    printf("The position of the car for t = 10s is %.4f.\n", p(10));
    printf("The velocity of the car for t = 10s is %.4f.\n", p.derivative(10));

    std::vector<double> values;
    double a = 0.0, b = 13.0;
    double h = (b - a) / 1000;
    for (double x = a; x <= b; x += h)
        values.push_back(p.derivative(x));
    auto maxn = std::max_element(values.begin(), values.end());
    double max_value = *maxn;
    printf("The maximum value of the velocity is %.4f.\n", max_value);

    return 0;
}