#ifndef EQUATION_SOLVER_H
#define EQUATION_SOLVER_H

#include <algorithm>
#include <iostream>
#include <functional>
#include <cmath>

const double pi = 3.14159265358979323846;

inline bool sgn(double a, double b) {
    return std::signbit(a) == std::signbit(b);
}

template <typename T>
class EquationSolver {
public:
    using Equation = std::function<T(T)>;
    virtual void solve(double tol_func, double tol_x, T a, T b) = 0;
    EquationSolver(Equation eq, Equation de = 0): equation_(eq), derivative_(de) {}

protected:
    Equation equation_, derivative_;
};

template <typename T>
class Bisection_Method: public EquationSolver<T> {
public:
    using Equation = std::function<T(T)>;
    Bisection_Method(Equation eq): EquationSolver<T>(eq) {}
    void solve(double tol_func, double tol_x, T a, T b) {
        double u = this->equation_(a), w;
        if (sgn(this->equation_(b), u)) {
            printf("The precondition does not meet, so the root is not found.\n");
            exit(0);
        }
        int M = 1000;
        T h = b - a;
        T c;
        for (int k = 0; k < M; ++k) {
            h /= 2;
            c = a + h;
            if (std::abs(h) < tol_x) 
                break;
            w = this->equation_(c);
            if (std::abs(w) < tol_func)
                break;
            else if (sgn(w, u))
                a = c;
        }
        printf("The root approximates to %.8f, and the error is %.8f.\n", c, this->equation_(c));
    }
};

template <typename T>
class Newton_Method: public EquationSolver<T> {
public:
    using Equation = std::function<T(T)>;
    Newton_Method(Equation eq, Equation de): EquationSolver<T>(eq, de) {}
    void solve(double tol_func, double tol_x, T a, T b = 0) override {
        double u, x0 = a;
        int M = 1000;
        for (int k = 0; k < M; ++k) {
            u = this->equation_(a);
            if (std::abs(u) < tol_func)
                break;
            a -= u/this->derivative_(a);
        }
        printf("The root of near %.8f approximates to %.8f, and the error is %.8f.\n", x0, a, this->equation_(a));
    }
};

template <typename T>
class Secant_Method: public EquationSolver<T> {
public:
    using Equation = std::function<T(T)>;
    Secant_Method(Equation eq): EquationSolver<T>(eq) {}
    void solve(double tol_func, double tol_x, T a, T b) override {
        double s, v = this->equation_(a), u = this->equation_(b);
        int M = 100000;
        for (int k = 2; k < M; ++k) {
            s = (b - a) / (u - v);
            a = b;
            v = u;
            b -= u * s;
            if (std::abs(b - a) < tol_x)
                break;
            u = this->equation_(b);
            if (std::abs(u) < tol_func)
                break;
        }
        printf("The root approximates to %.8f, and the error is %.8f.\n", b, this->equation_(b));
    }
};

#endif