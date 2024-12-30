#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include "Matrix.h"

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
    Polynomial operator* (const Polynomial& other) const {
        Polynomial result(degree + other.degree, std::vector<double>(degree + other.degree + 1));
        for (int i = 0; i <= degree; ++i) {
            for (int j = 0; j <= other.degree; ++j) {
                result.coeffs[i + j] += coeffs[i] * other.coeffs[j];
            }
        }
        return result;
    }
    Polynomial operator+ (const Polynomial& other) const {
        Polynomial result(std::max(degree, other.degree), std::vector<double>(std::max(degree, other.degree) + 1));
        for (int i = 0; i <= std::max(degree, other.degree); ++i) {
            result.coeffs[i] = coeffs[i] + other.coeffs[i];
        }
        return result;
    }

private:
    std::vector<double> coeffs;
    int degree;

};

#endif