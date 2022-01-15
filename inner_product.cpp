#include <iostream>
#include <stdio.h>

int main()
{
    double numArr[3][4] = {    // 세로 크기 3, 가로 크기 4인 int형 2차원 배열 선언
        { 11.6, 22.2, -33.3, 44.0 },
        { 55.2, 66.3, -77, 88.3 },
        { 99.1, 110.9, -121.2, -132.1 }
    };

    double* A = &numArr[0][0];

    int nrow = 3;
    int ncol = 4;
    double** mAt = (double**) malloc(sizeof(double*) * ncol);
    for(int i = 0; i < ncol; i ++) mAt[i] = A + i * nrow;


    std::cout << A << "\n";
    std::cout << *A << "\n\n";

    std::cout << mAt[0] <<"\n";
    std::cout << mAt[1][1] <<"\n";

    //printf("%p\n", A);     // print A
    //printf("%f\n", *A);    // print *A


    return 0;
}
