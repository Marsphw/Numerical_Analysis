#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <functional>
#include <vector>
#include <cmath>

const double pi = 3.14159265358979323846;

long long factorial(int n) {
    if (n == 0)
        return 1;
    return n * factorial(n - 1);
}

long long combination(int n, int k) {
    if (k == 0 || k == n)
        return 1;
    long long res = 1;
    int smaller_k = (k < n - k) ? k : (n - k);
    for (int i = 1; i <= smaller_k; ++i) 
        res *= (n - i + 1) / i;
    return res;
}

struct Point {
    double x;
    int order;
    double* f;
};

class Polynomial {
public:
    Polynomial(std::vector<double> coeff, std::vector<Point> points): coeff(coeff), points(points) {}
    double operator() (double x) const {
        double res = 0.0;
        int n = points.size(), N = n;
        for (int i = 0; i < n; ++i)
            N += points[i].order;
        int line = 0;
        double x_value[N + 1];
        for (int m = 0; m < n; ++m) {
            for (int j = 0; j <= points[m].order; ++j)
                x_value[line] = points[m].x, ++line;
        }
        for (int i = 0; i < N; ++i) {
            double temp = 1.0;
            for (int j = 0; j < i; ++j)
                temp *= (x - x_value[j]);
            temp *= coeff[i];
            res += temp;
        }
        return res;
    }

    double derivative(double x) const {
        double res = 0.0;
        int n = points.size(), N = n;
        for (int i = 0; i < n; ++i)
            N += points[i].order;
        int line = 0;
        double x_value[N + 1];
        for (int m = 0; m < n; ++m) {
            for (int j = 0; j <= points[m].order; ++j)
                x_value[line] = points[m].x, ++line;
        }
        for (int i = 0; i < N; ++i) {
            double temp1 = 0.0;
            for (int j = 0; j < i; ++j) {
                double temp2 = 1.0;
                for (int k = 0; k < i; ++k) {
                    if (j == k)
                        continue;
                    temp2 *= (x - x_value[k]);
                }
                temp1 += temp2;
            }
            temp1 *= coeff[i];
            res += temp1;
        }
        return res;
    }

private:
    std::vector<double> coeff;
    std::vector<Point> points;
};

class Curve {
public:
    Curve(std::vector<long long> coeff, std::vector<Point> points): coeff(coeff), points(points) {}
    Point operator() (double t) const {
        int n = coeff.size();
        double res_x = 0.0, res_y = 0.0;
        for (int i = 0; i < n; ++i) {
            res_x += coeff[i] * std::pow(t, i) * std::pow(1 - t, n - i - 1) * points[i].x;
            res_y += coeff[i] * std::pow(t, i) * std::pow(1 - t, n - i - 1) * points[i].f[0]; 
        }
        Point p;
        p.x = res_x;
        p.order = 0;
        p.f = new double[1];
        p.f[0] = res_y;
        return p;
    }

private:
    std::vector<long long> coeff;
    std::vector<Point> points;
};

class Interpolation {
public:
    Interpolation(std::vector<Point> points): points(points) {}

protected:
    std::vector<Point> points;
};

class Newton_Interpolation : public Interpolation {
public:
    Newton_Interpolation(std::vector<Point> points): Interpolation(points) {}
    Polynomial interpolate() {
        int n = points.size();
        double d[n + 1][n + 1];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j <= i; ++j) {
                if (j == 0)
                    d[i][j] = points[i].f[0];
                else 
                    d[i][j] = (d[i][j - 1] - d[i - 1][j - 1]) / (points[i].x - points[i - j].x);
            }
        }
        std::vector<double> coeff;
        for (int i = 0; i < n; ++i) 
            coeff.push_back(d[i][i]);
        return Polynomial(coeff, points);
    }
};

class Hermite_Interpolation : public Interpolation {
public:
    Hermite_Interpolation(std::vector<Point> points): Interpolation(points) {}
    Polynomial interpolate() {
        int n = points.size(), N = n;
        for (int i = 0; i < n; ++i)
            N += points[i].order;
        double d[N + 1][N + 1];
        int line = 0;
        double x_value[N + 1];
        for (int m = 0; m < n; ++m) {
            for (int j = 0; j <= points[m].order; ++j)
                x_value[line] = points[m].x, ++line;
        }
        line = 0;
        for (int m = 0; m < n; ++m) {
            for (int j = 0; j <= points[m].order; ++j) {
                for (int col = 0; col <= line; ++col) {
                    if (col <= j)
                        d[line][col] = points[m].f[col] / factorial(col);
                    else 
                        d[line][col] = (d[line][col - 1] - d[line - 1][col - 1]) / (points[m].x - x_value[line - col]);
                }
                ++line;
            }
        }
        std::vector<double> coeff;
        for (int i = 0; i < N; ++i)
            coeff.push_back(d[i][i]);
        return Polynomial(coeff, points);
    }
};

class Bezier_Curve: public Interpolation{
public:
    Bezier_Curve(std::vector<Point> points): Interpolation(points) {}
    Curve interpolate() {
        int n = points.size();
        std::vector<long long> coeff;
        for (int i = 0; i < n; ++i)
            coeff.push_back(combination(n - 1, i));
        return Curve(coeff, points);
    }
};

#endif