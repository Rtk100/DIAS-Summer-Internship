#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <random>
#include <complex>
#include <vector>
#include "eigen/Eigen/Dense"

const int N = 4;

typedef std::complex<long double> R_or_C;
typedef Eigen:: Matrix<std::complex<long double>, N, 1> col_vector;
typedef Eigen:: Matrix<std::complex<long double>, 1, N> row_vector;
typedef Eigen:: Matrix<std::complex<long double>, N, N> matrix;

R_or_C I(0,1);
int main()
{

    col_vector col;
    row_vector row;
    
    matrix mat = col * row;

    R_or_C scalar = (row*col)[0];

    R_or_C el(1.0, 2.2);


    R_or_C X1_4 = 1.0;
    R_or_C X1_5 = 1.2;
    R_or_C coeff(4, 0);
    R_or_C phi1_1 = coeff * (X1_4+I*X1_5);

    std::complex<long double> I = {0.0, 1.0};

    std::cout << phi1_1 << '\n';


    std::cout << I;




    return 0;
}