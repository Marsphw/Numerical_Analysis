#include <iostream>
#include <cstdio>
#include "ppForm.h"
#include "Bspline.h"

int main() {
    int N;
    printf("Enter the number of control points: ");
    scanf("%d", &N);
    std::vector<knot> knots(N + 1);
    printf("Enter the x and y coordinates of the control points:\n");
    for (int i = 0; i < N; i++) {
        knots[i].order = 0;
        scanf("%lf %lf", &knots[i].value[0], &knots[i].value[1]);
    }
    PPForm form(1, 0, N, knots);
    
}