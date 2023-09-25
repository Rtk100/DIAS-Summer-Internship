//This is a starter to the full D2 branes Lagrangian. 
// Only using Zs as the generalised coords, and only having the potential term equal V_D


#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <random>
#include <complex>
#include <vector>
#include "eigen/Eigen/Dense"
double start = std::time(nullptr);

long double c1 = 2.5;
long double c2 = 3.3;
long double c3 = 1.2;
long double c4 = -c1-c2-c3;

long double c12 = 1.1;
long double c13 = 2.2;
long double c23 = 2.0;
long double c14 = 0.7;
long double c24 = 3.1;
long double c34 = 0.9;

// Define timestep
const long double delta_t = 1e-4;
const long double seconds_thermalised = 100;

// Repeat simulation for 1000 seconds.
const int simulation_repetitions = seconds_thermalised / delta_t;
// Number of D0-Branes
const int N = 3;
const long double g = 1.0;
const int rows = N;
const int cols = N;

// Dimension of space
const int dim = 9;


//typedef double Complex;
//typedef Eigen:: Matrix<double, N, N> matrix;


//To go from real matrices to complex matrices delete the above typedefs and use these typedefs.
//And change the commented out code in generateHermitianMatrix()
typedef long double l_or_d;
typedef std::complex<long double> Complex;
typedef Eigen:: Matrix<std::complex<long double>, N, N> matrix;
typedef Eigen:: Matrix<std::complex<l_or_d>, 1, N> row_vector;
typedef Eigen:: Matrix<std::complex<l_or_d>, N, 1> col_vector;

// Initialise all Z generalised coords
static std::random_device rd;
double sigma = 1.0;
static std::mt19937 rng(std::time(nullptr)); 
std::normal_distribution<double> dist(0.0, sigma);

Complex Z12 = Complex((dist(rng),dist(rng)));
Complex Z13 = Complex((dist(rng),dist(rng)));

Complex Z21 = Complex((dist(rng),dist(rng)));
Complex Z23 = Complex((dist(rng),dist(rng)));

Complex Z31 = Complex((dist(rng),dist(rng)));
Complex Z32 = Complex((dist(rng),dist(rng)));

// Random below gives the same Z vectors every time unfortunately.
row_vector Z14 = row_vector::Random();
row_vector Z24 = row_vector::Random();
row_vector Z34 = row_vector::Random();

col_vector Z41 = col_vector::Random();
col_vector Z42 = col_vector::Random();
col_vector Z43 = col_vector::Random();


// Create a zero matrix in order to populate the V_vector with it.
matrix zero_matrix = matrix::Zero();
row_vector zero_row = row_vector::Zero();
col_vector zero_col = col_vector::Zero();

Complex Z12_dot = Complex(0.0, 0.0);
Complex Z13_dot = Complex(0.0, 0.0);
Complex Z21_dot = Complex(0.0, 0.0);
Complex Z23_dot = Complex(0.0, 0.0);
Complex Z31_dot = Complex(0.0, 0.0);
Complex Z32_dot = Complex(0.0, 0.0);

row_vector Z14_dot = zero_row;
row_vector Z24_dot = zero_row;
row_vector Z34_dot = zero_row;
col_vector Z41_dot = zero_col;
col_vector Z42_dot = zero_col;
col_vector Z43_dot = zero_col;

// Define H = K + V
Complex H(
    const long double g, 
    Complex Z12, Complex Z13, Complex Z21, Complex Z23, Complex Z31, Complex Z32,
    row_vector Z14, row_vector Z24, row_vector Z34,
    col_vector Z41, col_vector Z42, col_vector Z43,
    Complex Z12_dot, Complex Z13_dot, Complex Z21_dot, Complex Z23_dot, Complex Z31_dot, Complex Z32_dot,
    row_vector Z14_dot, row_vector Z24_dot, row_vector Z34_dot, 
    col_vector Z41_dot, col_vector Z42_dot, col_vector Z43_dot)
{   

    Complex K = Complex(0.5, 0) * (std::conj(Z12_dot)*Z12_dot + std::conj(Z13_dot)*Z13_dot +
                                   std::conj(Z21_dot)*Z21_dot + std::conj(Z23_dot)*Z23_dot + 
                                   std::conj(Z31_dot)*Z31_dot + std::conj(Z32_dot)*Z32_dot + 
                                   (Z14_dot.adjoint()*Z14_dot).trace() + (Z24_dot.adjoint()*Z24_dot).trace() + 
                                   (Z34_dot.adjoint()*Z34_dot).trace() + (Z41_dot.adjoint()*Z41_dot).trace() + 
                                   (Z42_dot.adjoint()*Z42_dot).trace() + (Z43_dot.adjoint()*Z43_dot).trace() );

    Complex V_D_arg_k1 = (Z12 * std::conj(Z12) + Z13 * std::conj(Z13) + (Z14 * Z14.adjoint()).trace() -
                         (std::conj(Z21) * Z21 + std::conj(Z31) * Z31 + (Z41.adjoint() * Z41).trace() )
                         - c1/(g*g));

    Complex V_D_arg_k2 = (Z21 * std::conj(Z21) + Z23 * std::conj(Z23) + (Z24 * Z24.adjoint()).trace() - 
                         (std::conj(Z12) * Z12 + std::conj(Z32) * Z32 + (Z42.adjoint() * Z42).trace() )
                         - c2/(g*g) );

    Complex V_D_arg_k3 = (Z31 * std::conj(Z31) + Z32 * std::conj(Z32) + (Z34 * Z34.adjoint()).trace() - 
                         (std::conj(Z13) * Z13 + std::conj(Z23) * Z23 + (Z43.adjoint() * Z43).trace() ) 
                         - c3/(g*g) );

    matrix V_D_arg_k4 = (Z41 * Z41.adjoint() + Z42 * Z42.adjoint() + Z43 * Z43.adjoint() + 
                        (Z14.adjoint() * Z14 + Z24.adjoint() * Z24 + Z34.adjoint() * Z34 ) 
                        - c4/(g*g) * matrix::Identity() );

    Complex V_D = Complex(0.5, 0) * ((V_D_arg_k1 * V_D_arg_k1) + (V_D_arg_k2 * V_D_arg_k2) + 
             (V_D_arg_k3 * V_D_arg_k3) + (V_D_arg_k4 * V_D_arg_k4).trace() );

    Complex V = g*g*V_D;

    return K + V;
}

int main() 
{
/*
// For testing reproducibility use these X values
    std::ifstream inputX("initial_X.txt");
    if (!inputX.is_open()) {
        std::cerr << "Failed to open X initial file." << std::endl;
        return 1;
    }

    // Read the values from the file and store them in the matrices
    for (int i = 0; i < dim; ++i) 
    {

        for (int row = 0; row < rows; ++row) 
        {
            for (int col = 0; col < cols; ++col) 
            {
                 inputX >> X_vector[i](row, col);
            }
        }
    }

    // Close the input file
    inputX.close();



    // Generate and store A1, A2, A3, A4, A5, A6, A7, A8, and A9
    for (int i = 0; i < dim; ++i) 
    {
        A_vector[i] = Acceleration(i, X_vector, rows, cols, g);
    }  

*/
    // Define the initial accelerations from the equations of motion

/* 
There are terms in the equations of motion that appear a bunch of times. These are the bracketed terms that are most easily spotted by
looking for the c^k. All bracketed terms with c^1 are the same, and that stands for all k. So I define them here to more easily write out
the equations of motion.
*/ 
    Complex c1_term = abs(Z12)*abs(Z12) + abs(Z13)*abs(Z13) + (Z14 * Z14.adjoint())[0] - abs(Z21)*abs(Z21) - abs(Z31)*abs(Z31) - (Z41.adjoint() * Z41)[0] - c1/(g*g);
    Complex c2_term = abs(Z21)*abs(Z21) + abs(Z23)*abs(Z23) + (Z24 * Z24.adjoint())[0] - abs(Z12)*abs(Z12) - abs(Z32)*abs(Z32) - (Z42.adjoint() * Z42)[0] - c2/(g*g);
    Complex c3_term = abs(Z31)*abs(Z31) + abs(Z32)*abs(Z32) + (Z34 * Z34.adjoint())[0] - abs(Z13)*abs(Z13) - abs(Z23)*abs(Z23) - (Z43.adjoint() * Z43)[0] - c3/(g*g);
    matrix c4_term = (Z41 * Z41.adjoint()) + (Z42 * Z42.adjoint()) + (Z43 * Z43.adjoint()) - (Z14.adjoint() * Z14) - (Z24.adjoint() * Z24) - (Z34.adjoint() * Z34) - c4/(g*g) * matrix::Identity();
    

    // Define the accelerations of Z

    Complex Z12_double_dot = -Z12 * c1_term + Z12 * c2_term;
    Complex Z21_double_dot = +Z21 * c1_term - Z21 * c2_term;
    Complex Z13_double_dot = -Z13 * c1_term + Z13 * c3_term;
    Complex Z31_double_dot = +Z31 * c1_term - Z31 * c3_term;
    Complex Z23_double_dot = -Z23 * c2_term + Z23 * c3_term;
    Complex Z32_double_dot = +Z32 * c2_term - Z32 * c3_term;

    row_vector Z14_double_dot = -Z14 * c1_term + Z14 * c4_term;
    row_vector Z41_double_dot =  Z41 * c1_term - c4_term * Z41;
    row_vector Z24_double_dot = -Z24 * c2_term + Z24 * c4_term;
    col_vector Z42_double_dot =  Z42 * c1_term - c4_term * Z42;
    col_vector Z34_double_dot = -Z34 * c3_term + Z34 * c4_term;
    col_vector Z43_double_dot =  Z43 * c3_term - c4_term * Z43;



    Complex scalar_Zs[6] = {Z12, Z13, Z21, Z23, Z31, Z32};
    row_vector row_Zs[3] = {Z14, Z24, Z34};
    col_vector col_Zs[3] = {Z41, Z42, Z43};

    Complex scalar_Zs_new[6] = {Complex(0.0, 0.0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)};
    row_vector row_Zs_new[3] = {zero_row, zero_row, zero_row};
    col_vector col_Zs_new[3] = {zero_col, zero_col, zero_col};

    Complex scalar_Z_dots[6] = {Z12_dot, Z13_dot, Z21_dot, Z23_dot, Z31_dot, Z32_dot};
    row_vector row_Z_dots[3] = {Z14_dot, Z24_dot, Z34_dot};
    col_vector col_Z_dots[3] = {Z41_dot, Z42_dot, Z43_dot};

    Complex scalar_Z_dots_new[6] = {Complex(0.0, 0.0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0), Complex(0, 0)};
    row_vector row_Z_dots_new[3] = {zero_row, zero_row, zero_row};
    col_vector col_Z_dots_new[3] = {zero_col, zero_col, zero_col};

    Complex scalar_Z_double_dots[6] = {Z12_double_dot, Z13_double_dot, Z21_double_dot, Z23_double_dot, Z31_double_dot, Z32_double_dot};
    row_vector row_Z_double_dots[3] = {Z14_double_dot, Z24_double_dot, Z34_double_dot};
    col_vector col_Z_double_dots[3] = {Z41_double_dot, Z42_double_dot, Z43_double_dot};

    Complex scalar_Z_double_dots_new[6] = {Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)};
    row_vector row_Z_double_dots_new[3] = {zero_row, zero_row, zero_row};
    col_vector col_Z_double_dots_new[3] = {zero_col, zero_col, zero_col};


    Complex c1_term_new = abs(scalar_Zs_new[0])*abs(scalar_Zs_new[0]) + abs(scalar_Zs_new[1])*abs(scalar_Zs_new[1]) + (row_Zs_new[0] * row_Zs_new[0].adjoint())[0] - abs(scalar_Zs_new[2])*abs(scalar_Zs_new[2]) - abs(scalar_Zs_new[4])*abs(scalar_Zs_new[4]) - (col_Zs_new[0].adjoint() * col_Zs_new[0])[0] - c1/(g*g);
    Complex c2_term_new = abs(scalar_Zs_new[2])*abs(scalar_Zs_new[2]) + abs(scalar_Zs_new[3])*abs(scalar_Zs_new[3]) + (row_Zs_new[1] * row_Zs_new[1].adjoint())[0] - abs(scalar_Zs_new[0])*abs(scalar_Zs_new[0]) - abs(scalar_Zs_new[5])*abs(scalar_Zs_new[5]) - (col_Zs_new[1].adjoint() * col_Zs_new[1])[0] - c2/(g*g);
    Complex c3_term_new = abs(scalar_Zs_new[4])*abs(scalar_Zs_new[4]) + abs(scalar_Zs_new[5])*abs(scalar_Zs_new[5]) + (row_Zs_new[2] * row_Zs_new[2].adjoint())[0] - abs(scalar_Zs_new[1])*abs(scalar_Zs_new[1]) - abs(scalar_Zs_new[3])*abs(scalar_Zs_new[3]) - (col_Zs_new[2].adjoint() * col_Zs_new[2])[0] - c3/(g*g);
    matrix c4_term_new = (col_Zs_new[0] * col_Zs_new[0].adjoint()) + (col_Zs_new[1] * col_Zs_new[1].adjoint()) + (col_Zs_new[2] * col_Zs_new[2].adjoint()) - (row_Zs_new[0].adjoint() * row_Zs_new[0]) - (row_Zs_new[1].adjoint() * row_Zs_new[1]) - (row_Zs_new[2].adjoint() * row_Zs_new[2]) - c4/(g*g) * matrix::Identity();

    // Write simulation to thermalise system
    for (int j = 0; j < seconds_thermalised / delta_t; ++j)
    {

        // velocity Verlet 1 to get new positions from old positions, momentums and rate of change of momentums

        for (int i = 0; i < 6; ++i)
        {
            scalar_Zs_new[i] = scalar_Zs[i] + scalar_Z_dots[i] * delta_t + Complex(0.5, 0) * scalar_Z_double_dots[i] * delta_t * delta_t;
        }
        for (int i = 0; i < 3; ++i)
        {
            row_Zs_new[i] = row_Zs[i] + row_Z_dots[i] * delta_t + Complex(0.5, 0) * row_Z_double_dots[i] * delta_t * delta_t;
        }
        for (int i = 0; i < 3; ++i)
        {
            col_Zs_new[i] = col_Zs[i] + col_Z_dots[i] * delta_t + Complex(0.5, 0) * col_Z_double_dots[i] * delta_t * delta_t;
        }

        // Generate and store new A1, A2, A3, A4, A5, A6, A7, A8, and A9
/*
    matrix Z14_double_dot = -Z14 * c1_term + Z14 * c4_term;
    matrix Z41_double_dot =  Z41 * c1_term - c4_term * Z41;
    matrix Z24_double_dot = -Z24 * c2_term + Z24 * c4_term;
    matrix Z42_double_dot =  Z42 * c1_term - c4_term * Z42;
    matrix Z34_double_dot = -Z34 * c3_term + Z34 * c4_term;
    matrix Z43_double_dot =  Z43 * c3_term - c4_term * Z43;

    Complex scalar_Zs[6] = {Z12, Z13, Z21, Z23, Z31, Z32};
    row_vector row_Zs[3] = {Z14, Z24, Z34};
    col_vector col_Zs[3] = {Z41, Z42, Z43};
*/
        scalar_Z_double_dots_new[0] = -scalar_Zs_new[0] * c1_term_new + scalar_Zs_new[0] * c2_term_new;
        scalar_Z_double_dots_new[1] = +scalar_Zs_new[2] * c1_term_new - scalar_Zs_new[2] * c2_term_new;
        scalar_Z_double_dots_new[2] = -scalar_Zs_new[1] * c1_term_new + scalar_Zs_new[1] * c3_term_new;
        scalar_Z_double_dots_new[3] = +scalar_Zs_new[4] * c1_term_new - scalar_Zs_new[4] * c3_term_new;
        scalar_Z_double_dots_new[4] = -scalar_Zs_new[3] * c2_term_new + scalar_Zs_new[3] * c3_term_new;
        scalar_Z_double_dots_new[5] = +scalar_Zs_new[5] * c2_term_new - scalar_Zs_new[5] * c3_term_new;
        row_Z_double_dots_new[0] = -row_Zs_new[0] * c1_term_new + row_Zs_new[0] * c4_term_new;
        row_Z_double_dots_new[1] =  col_Zs_new[0] * c1_term_new - c4_term_new * col_Zs_new[0];
        row_Z_double_dots_new[2] = -row_Zs_new[1] * c2_term_new + row_Zs_new[1] * c4_term_new;
        col_Z_double_dots_new[0] =  col_Zs_new[1] * c1_term_new - c4_term_new * col_Zs_new[1];
        col_Z_double_dots_new[1] = -row_Zs_new[2] * c3_term_new + row_Zs_new[2] * c4_term_new;
        col_Z_double_dots_new[2] =  col_Zs_new[2] * c3_term_new - c4_term_new * col_Zs_new[2];


        // Use Velocity Verlet 2 to get new momentums
        for (int i = 0; i < 6; ++i)
        {
            scalar_Z_dots_new[i] = scalar_Z_dots[i] + Complex(0.5, 0) * (scalar_Z_double_dots_new[i] + scalar_Z_double_dots[i]) * delta_t;
        }
        for (int i = 0; i < 3; ++i)
        {
            row_Z_dots_new[i] = row_Z_dots[i] + 0.5 * (row_Z_double_dots_new[i] + row_Z_double_dots[i]) * delta_t;
        }
        for (int i = 0; i < 3; ++i)
        {
            col_Z_dots_new[i] = col_Z_dots[i] + 0.5 * (col_Z_double_dots_new[i] + col_Z_double_dots[i]) * delta_t;
        }
        // Copy elements from X_vector_new to X_vector
        std::memcpy(scalar_Zs, scalar_Zs_new, sizeof(scalar_Zs_new));  
        std::memcpy(row_Zs, row_Zs_new, sizeof(row_Zs_new));  
        std::memcpy(col_Zs, col_Zs_new, sizeof(col_Zs_new));  

        // Copy elements from V_vector_new to V_vector
        std::memcpy(scalar_Z_dots, scalar_Z_dots_new, sizeof(scalar_Z_dots_new));  
        std::memcpy(row_Z_dots, row_Z_dots_new, sizeof(row_Z_dots_new));  
        std::memcpy(col_Z_dots, col_Z_dots_new, sizeof(col_Z_dots_new));  

        // Copy elements from A_vector_new to A_vector
        std::memcpy(scalar_Z_double_dots, scalar_Z_double_dots_new, sizeof(scalar_Z_double_dots_new));  
        std::memcpy(row_Z_double_dots, row_Z_double_dots_new, sizeof(row_Z_double_dots_new));  
        std::memcpy(col_Z_double_dots, col_Z_double_dots_new, sizeof(col_Z_double_dots_new));  

        if (j % 1000 == 0)
        {
            //for (matrix el : V_vector)
            //{
            //    std::cout <<"Ideal " << el << std::endl;
            //}
            //std::cout  << std::endl << gauss_law(X_vector_new[0], X_vector_new[1], X_vector_new[2], X_vector_new[3], X_vector_new[4], X_vector_new[5], X_vector_new[6], X_vector_new[7], X_vector_new[8],
            //                       V_vector_new[0], V_vector_new[1], V_vector_new[2], V_vector_new[3], V_vector_new[4], V_vector_new[5], V_vector_new[6], V_vector_new[7], V_vector_new[8]);

            std::cout << std::endl;
            std::cout << "H" << std::setprecision(15) << H(g, 
                            scalar_Zs_new[0], scalar_Zs_new[1], scalar_Zs_new[2], scalar_Zs_new[3], scalar_Zs_new[4], scalar_Zs_new[5],
                            row_Zs_new[0], row_Zs_new[1], row_Zs_new[2],
                            col_Zs_new[0], col_Zs_new[1], col_Zs_new[2],
                            scalar_Z_dots_new[0], scalar_Z_dots_new[1], scalar_Z_dots_new[2], scalar_Z_dots_new[3], scalar_Z_dots_new[4], scalar_Z_dots_new[5],
                            row_Z_dots_new[0], row_Z_dots_new[1], row_Z_dots_new[2],
                            col_Z_dots_new[0], col_Z_dots_new[1], col_Z_dots_new[2]);


        }

    }



/*
       // Export initial X/V/A_vector to text files to be analysed in python.
    std:: fstream X2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D-BranesNis4/thermalised_X.txt", std:: ios:: out);
    X2_vector_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (matrix Matrix : X_vector_new)
    {
        X2_vector_Export << Matrix << std::endl;
    }


    std:: fstream V2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D-BranesNis4/thermalised_V.txt", std:: ios:: out);
    V2_vector_Export << std::fixed << std::setprecision(15);
    for (matrix Matrix : V_vector_new)
    {
        V2_vector_Export << Matrix << std::endl;
    }

    std:: fstream A2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D-BranesNis4/thermalised_A.txt", std:: ios:: out);
    A2_vector_Export << std::fixed << std::setprecision(15);
    // Print to text file
    for (matrix Matrix : A_vector_new)
    {
        A2_vector_Export << Matrix << std::endl;
    }


    X2_vector_Export.close();
    V2_vector_Export.close();
    A2_vector_Export.close();

*/



    std::cout << "Finished" << std::time(nullptr)-start;

    return 0;
}
