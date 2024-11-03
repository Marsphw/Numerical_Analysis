#include "Interpolation.h"
#include <iostream>
#include <cstdio>
#include <vector>

int main() {

    freopen("../results/Task_E/Task_E_scatter.txt", "w", stdout);
    std::vector<Point> points_Sp1, points_Sp2;

    for (int i = 0; i < 7; ++i) {
        double x, f;
        std::cin >> x >> f;
        Point temp;
        temp.x = x;
        temp.order = 0;
        temp.f = new double[1];
        temp.f[0] = f;
        points_Sp1.push_back(temp);
    }
    Newton_Interpolation interpolation1(points_Sp1);
    Polynomial p1 = interpolation1.interpolate();
    for (int j = 0; j < 1000; ++j) {
        double x = 28.0*j / 1000.0;
        double y = p1(x);
        std::cout << y << std::endl;
    }

    for (int i = 0; i < 7; ++i) {
        double x, f;
        std::cin >> x >> f;
        Point temp;
        temp.x = x;
        temp.order = 0;
        temp.f = new double[1];
        temp.f[0] = f;
        points_Sp2.push_back(temp);
    }
    Newton_Interpolation interpolation2(points_Sp2);
    Polynomial p2 = interpolation2.interpolate();
    for (int j = 0; j < 1000; ++j) {
        double x = 28.0*j / 1000.0;
        double y = p2(x);
        std::cout << y << std::endl;
    }

    fclose(stdout);

    freopen("../results/Task_E/Task_E_result.txt", "w", stdout);
    printf("The average weight predicted of Sp1 is %.4f\n", p1(43.0));
    printf("The average weight predicted of Sp2 is %.4f\n", p2(43.0));
    fclose(stdout);

    return 0;
}