#include "Interpolation.h"
#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

double f1(double x) {
    return 2.0 * (std::sqrt(std::abs(x)) + std::sqrt(3 - x*x)) / 3;
}

double f2(double x) {
    return 2.0 * (std::sqrt(std::abs(x)) - std::sqrt(3 - x*x)) / 3;
}

Point f1_derivative(double x) {
    Point res;
    double derivative;
    res.order = 0;
    res.f = new double[1];
    if (x == 0) {
        derivative = (f1(1e-4) - f1(0.0)) / 1e-4;
        res.x = 1.0 / std::sqrt(1 + derivative * derivative);
        res.f[0] = std::abs(derivative) / std::sqrt(1 + derivative * derivative);
        return res;
    }
    if (std::abs(x) - std::sqrt(3) < 1e-3) {
        res.x = 0.0;
        res.f[0] = 1.0;
        return res;
    }
    derivative = 2.0/3 * (1.0/(2 * std::sqrt(std::abs(x))) - x / std::sqrt(3 - x*x));
    res.x = 1.0 / std::sqrt(1 + derivative * derivative);
    res.f[0] = std::abs(derivative) / std::sqrt(1 + derivative * derivative);
    return res;
}

Point f2_derivative(double x) {
    Point res;
    double derivative;
    res.order = 0;
    res.f = new double[1];
    if (x == 0) {
        derivative = (f2(1e-4) - f2(0.0)) / 1e-4;
        res.x = 1.0 / std::sqrt(1 + derivative * derivative);
        res.f[0] = std::abs(derivative) / std::sqrt(1 + derivative * derivative);
        return res;
    }
    derivative = 2.0/3 * (1.0/(2 * std::sqrt(std::abs(x))) + x / std::sqrt(3 - x*x));
    res.x = 1.0 / std::sqrt(1 + derivative * derivative);
    res.f[0] = std::abs(derivative) / std::sqrt(1 + derivative * derivative);
    return res;
}

inline int sign(double x){
    return (x >= 0) - (x < 0);
}

int main() {
    int m[4] = {10, 40, 160};
    double p[7][3];
    p[0][0] = 0;
    p[0][1] = 2 / std::sqrt(3);
    p[1][0] = std::sqrt(3);
    p[1][1] = 2 * std::sqrt(std::sqrt(3)) / 3;
    p[2][0] = (std::sqrt(13) - 1) / 2;
    p[2][1] = 0;
    p[3][0] = 0;
    p[3][1] = -2 / std::sqrt(3);
    p[4][0] = -p[2][0];
    p[4][1] = 0;
    p[5][0] = -std::sqrt(3);
    p[5][1] = 2 * std::sqrt(std::sqrt(3)) / 3;
    std::vector<Point> marker_points;
    for (int i = 0; i < 3; ++i) {
        std::string filename = "../results/Task_F/scatter_plot_" + std::to_string(m[i]) + ".txt";
        std::ofstream outfile(filename);
        outfile << m[i] << std::endl;
        marker_points.clear();
        for (int j = 0; j < 6; ++j) {
            Point temp;
            temp.x = p[j][0];
            temp.order = 0;
            temp.f = new double[1];
            temp.f[0] = p[j][1];
            marker_points.push_back(temp);
            double h = (p[(j + 1) % 6][0] - p[j][0]) / ((m[i] - 5 + j)/6 + 1);
            for (int k = 1; k <= (m[i] - 5 + j) / 6; ++k) {
                Point t;
                t.x = p[j][0] + k * h;
                t.order = 0;
                t.f = new double[1];
                if (j == 0 || j == 5)
                    t.f[0] = f1(t.x);
                else
                    t.f[0] = f2(t.x);
                marker_points.push_back(t);
            }
        }
        for (int j = 0; j < m[i]; ++j) {
            std::vector<Point> points;
            points.clear();
            Point q[5];
            q[0].x = marker_points[j].x;
            q[0].order = 0;
            q[0].f = new double[1];
            q[0].f[0] = marker_points[j].f[0];
            q[3].x = marker_points[(j + 1)%m[i]].x;
            q[3].order = 0 ;
            q[3].f = new double[1];
            q[3].f[0] = marker_points[(j + 1)%m[i]].f[0];
            if (q[0].f[0] >= 2 * std::sqrt(std::sqrt(3)) / 3) {
                Point tangent = f1_derivative(q[0].x);
                q[1].x = q[0].x + 1.0/3 * tangent.x * sign(q[3].x - q[0].x);
                q[1].order = 0;
                q[1].f = new double[1];
                q[1].f[0] = q[0].f[0] + 1.0/3 * tangent.f[0] * sign(q[3].f[0] - q[0].f[0]);
                q[2].x = q[3].x - 1.0/3 * tangent.x * sign(q[3].x - q[0].x);
                q[2].order = 0;
                q[2].f = new double[1];
                q[2].f[0] = q[3].f[0] - 1.0/3 * tangent.f[0] * sign(q[3].f[0] - q[0].f[0]);
            } else {
                Point tangent = f2_derivative(q[0].x);
                q[1].x = q[0].x + 1.0/3 * tangent.x * sign(q[3].x - q[0].x);
                q[1].order = 0;
                q[1].f = new double[1];
                q[1].f[0] = q[0].f[0] + 1.0/3 * tangent.f[0] * sign(q[3].f[0] - q[0].f[0]);
                q[2].x = q[3].x - 1.0/3 * tangent.x * sign(q[3].x - q[0].x);
                q[2].order = 0;
                q[2].f = new double[1];
                q[2].f[0] = q[3].f[0] - 1.0/3 * tangent.f[0] * sign(q[3].f[0] - q[0].f[0]);
            } 
            for (int k = 0; k < 4; ++k) 
                points.push_back(q[k]);
            Bezier_Curve curve(points);
            Curve c = curve.interpolate();
            for (double t = 0.0; t < 1.01; t += 0.01) {
                Point p = c(t);
                outfile << p.x << std::endl << p.f[0] << std::endl;
            }
        }
        outfile.close();
    }
    return 0;
}