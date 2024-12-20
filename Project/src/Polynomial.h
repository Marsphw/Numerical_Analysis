#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>

class Polynomial {

public:
    Polynomial(int de, std::vector<double> co): degree(de), coeffs(co) {}
    double operator() (double x) const {
        double result = 0;
        for (int i = 0; i <= degree; ++i)
            result += coeffs[i] * pow(x, i);
        return result;
    }
    double derivative(double x) const {
        double result = 0;
        for (int i = 1; i <= degree; ++i) 
            result += i * coeffs[i] * pow(x, i - 1);
        return result;
    }

private:
    std::vector<double> coeffs;
    int degree;

};

#endif