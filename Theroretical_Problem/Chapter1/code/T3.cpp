#include <iostream>
#include <cmath>

int main() {
    double x[10];
    x[0] = -1.0;
    for (int i = 1; i < 10; i++) 
        x[i] = x[i - 1] - (4*pow(x[i - 1], 3) - 2*pow(x[i - 1], 2) + 3) / (12*pow(x[i - 1], 2) - 4*x[i - 1]);
    for (int i = 0; i < 10; i++) 
        std::cout << x[i] << " ";
    std::cout << std::endl;
    for (int i = 0; i < 10; i++)
        std::cout << 4*pow(x[i], 3) - 2*pow(x[i], 2) + 3 << " ";
    return 0;
}