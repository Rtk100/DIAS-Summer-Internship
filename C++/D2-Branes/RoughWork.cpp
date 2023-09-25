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

long double delta_t = 0.1;

R_or_C I(0,1);
int main()
{
    std:: fstream Z_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/Z.txt", std:: ios:: out);
    Z_Export << std::fixed << std::setprecision(15);

    std:: fstream Z_dot_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/Z_dot.txt", std:: ios:: out);
    Z_dot_Export << std::fixed << std::setprecision(15);

    std:: fstream A_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/A.txt", std:: ios:: out);
    A_Export << std::fixed << std::setprecision(15);

    R_or_C Z = {1.0, 2.0};
    R_or_C Z_dot = 0;
    R_or_C Z_dot_new, Z_new;
    R_or_C A, A_new;
    A = -Z;
    R_or_C coeff = 0.5;
    for (int j = 0; j < 100 / delta_t; ++j)
    {

        // velocity Verlet 1 to get new positions from old positions, momentums and rate of change of momentums

        
        Z_new = Z + Z_dot * delta_t + coeff * A * delta_t * delta_t;
        

        // Generate and store new A1, A2, A3, A4, A5, A6, A7, A8, and A9
       
        A_new = -Z_new;
        
        Z_dot_new = Z_dot + coeff * (A_new + A) * delta_t;

        

        // Copy elements from X_vector_new to X_vector
        Z = Z_new;
        Z_dot = Z_dot_new;
        A = A_new;


        Z_Export << Z << "\n";
        Z_dot_Export << Z_dot<< "\n";
        A_Export << A<< "\n";
        
    }
        R_or_C Z13 = R_or_C(0.5, 2);
        row_vector Z14;
        Z14(0,0) = R_or_C(1,1);
        Z14(0,1) = R_or_C(2,1);
        Z14(0,2) = R_or_C(3,1);
        Z14(0,3) = R_or_C(4,1);

        col_vector Z41;
        Z41(0,0) = R_or_C(1,1);
        Z41(1,0) = R_or_C(2,1);
        Z41(2,0) = R_or_C(3,1);
        Z41(3,0) = R_or_C(4,1);

        std::cout << Z14 << '\n' << (Z14 * Z14.adjoint());
        std::cout << Z41 << '\n' << 34 + (Z41.adjoint() * Z41)[0].real();
/*
        long double Z = ((Z14 * Z14.adjoint())).real();

        abs(Z13)*abs(Z13) + ((Z14 * Z14.adjoint())).real();
        abs(Z13)*abs(Z13) + ((Z41.adjoint() * Z41)).real();

        std::cout << Z14 << '\n' << (Z14 * Z14.adjoint());
*/



    return 0;
}