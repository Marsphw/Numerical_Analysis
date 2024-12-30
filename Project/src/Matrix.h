#ifndef MATRIX_H
#define MATRIX_H

#include <cstdio>
#include <vector>

template <typename T>
class Matrix {

public:
    Matrix(int r, int c): rows(r), cols(c), data(r, std::vector<T>(c, 0)) {}
    Matrix(const Matrix<T>& other): rows(other.rows), cols(other.cols), data(other.data) {}
    Matrix<T>& operator=(const Matrix<T>& other) {
        if (this == &other) 
            return *this;
        rows = other.rows;
        cols = other.cols;
        data = other.data;
        return *this;
    }
    int getRows() const { return rows; }
    int getCols() const { return cols; }
    T& operator()(int i, int j) { 
        return data[i][j]; 
    }
    void print() const {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                printf("%d ", data[i][j]);
            }
            printf("\n");
        }
    }

    Matrix<T> operator+(const Matrix<T>& other) const {
        if (rows != other.rows || cols != other.cols) {
            printf("Error: matrices must have the same size\n");
            return *this;
        }
        Matrix<T> result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result(i, j) = data[i][j] + other(i, j);
            }
        }
        return result;
    }

    Matrix<T> operator*(const Matrix<T>& other) const {
        if (cols != other.rows) {
            printf("Error: matrices must have compatible sizes\n");
            return *this;
        }
        Matrix<T> result(rows, other.cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < other.cols; j++) {
                for (int k = 0; k < cols; k++) {
                    result(i, j) += data[i][k] * other(k, j);
                }
            }
        return result;
        }
    }

private:
    std::vector<std::vector<T>> data;
    int rows;
    int cols;

};

#endif