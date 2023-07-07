#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <random>
#include <complex>
#include <vector>
#include "eigen/Eigen/Dense"

// Define timestep
const double delta_t = 1e-4;
const double seconds_thermalised = 1;
const double g = 1;

// Repeat simulation for 1000 seconds.
const int simulation_repetitions = seconds_thermalised / delta_t;
// Number of D0-Branes
const int N = 2;
const int rows = N;
const int cols = N;

// Dimension of space
const int dim = 9;


//typedef double R_or_C;
//typedef Eigen:: Matrix<double, N, N> matrix;


//To go from real matrices to complex matrices delete the above typedefs and use these typedefs.
//And change the commented out code in generateHermitianMatrix()

typedef std::complex<double> R_or_C;
typedef Eigen:: Matrix<std::complex<double>, N, N> matrix;




matrix commutator(matrix A, matrix B)
{
    return A * B - B * A;
} 


// Cillian's Hamiltonian
double H(
    double g, 
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    // Compute kinetic energy T
    R_or_C T = 1.0/(2.0 * pow(g,2)) * (V1 * V1 + V2 * V2 + V3 * V3 + V4 * V4 + V5 * V5 + V6 * V6 + V7 * V7 + V8 * V8 + V9 * V9).trace();

    matrix X[9] = {X1,X2,X3,X4,X5,X6,X7,X8,X9}; 

    matrix commutator_sum = matrix::Zero(rows, cols);  
    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            if(i != j)
            {
                commutator_sum += commutator(X[i],X[j])*commutator(X[i],X[j]); //can likely be more efficient by less function calls

            }
        }
    }
    R_or_C U = - 1.0/(4.0*pow(g,2)) * commutator_sum.trace();



    return std:: abs(T + U);
}

// Cillian's Gauss' law
matrix gauss_law(
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    matrix result;
    for (int i = 0; i < 9; i ++)
    {
        result = commutator(X1,V1) + commutator(X2,V2) + commutator(X3,V3) + commutator(X4,V4) + commutator(X5,V5) +
        commutator(X6,V6) + commutator(X7,V7) + commutator(X8,V8) + commutator(X9,V9);
    }
    return result;
}


// Random Hermitian matrix

static std::random_device rd;

static std::mt19937 rng(std::time(nullptr)); 
std::normal_distribution<double> dist(0.0, 0.1);

matrix generateHermitianMatrix(int rows, int cols) 
{   
    matrix Matrix(rows, cols);

    for (int i = 0; i < rows; ++i) 
    {
        for (int j = 0; j <= i; ++j) 
        {
            //double element = dist(rng);


            // For Complex elements Uncomment this
            double real = dist(rng);
            double imag = dist(rng);
            R_or_C element(real, imag);
            
            if (i == j)
            {
                Matrix(i,j) = real;
            }

            if (i != j) 
            {
                Matrix(i,j) = element;
                Matrix(j,i) = std::conj(element);
            }
        }
    }
    // Make matrix traceless by replacing last entry with the negative of the sum of the other entries.
    // Calculate the sum of diagonal elements in X1
    R_or_C diagonalSum = 0.0;
    for (int i = 0; i < rows-1; ++i) 
    {
        diagonalSum += Matrix(i,i);
    }
    Matrix(rows - 1,cols - 1) = -diagonalSum;

    return Matrix;
}


// Acceleration of each coordinate matrix
matrix Acceleration(const int j, matrix* X_vector, int rows, int cols, const double g)
{
    matrix commutator_sum = matrix::Zero(rows, cols);

    matrix X = X_vector[j];
    for (int i = 0; i < dim; ++i)
    {
        if (i != j)
        {  
            matrix temp_commutator = commutator(X_vector[i], commutator(X, X_vector[i]));
            

            commutator_sum += temp_commutator;
            
        }   
    // Comment out the line below to go back to the lagrangian e.q.m. With no -1/(g^2)
    // commutator_sum = -1.0 / (g * g) * commutator_sum;
    }
    return commutator_sum;
}

// Acceleration of each coordinate matrix
matrix Acceleration2(const int i, matrix* X_vector, int rows, int cols, const double g)
{
    matrix commutator_sum = matrix::Zero(rows, cols);
    matrix temp_commutator = matrix::Zero(rows, cols);

    matrix X = X_vector[i];
    double a1 = X(0,0).real();
    double b1 = X(1,1).real();
    double c1 = X(2,2).real();
    double d1 = X(0,1).real();
    double e1 = X(0,1).imag();
    double f1 = X(0,2).real();
    double g1 = X(0,2).imag();
    double h1 = X(0,3).real();
    double j1 = X(0,3).imag();
    double k1 = X(1,2).real();
    double l1 = X(1,2).imag();
    double m1 = X(1,3).real();
    double n1 = X(1,3).imag();
    double o1 = X(2,3).real();
    double p1 = X(2,3).imag();
    for (int j = 0; j < dim; ++j)
    {
        if (i != j)
        {  
            matrix X_other = X_vector[j];
            double a2 = X_other(0,0).real();
            double b2 = X_other(1,1).real();
            double c2 = X_other(2,2).real();
            double d2 = X_other(0,1).real();
            double e2 = X_other(0,1).imag();
            double f2 = X_other(0,2).real();
            double g2 = X_other(0,2).imag();
            double h2 = X_other(0,3).real();
            double j2 = X_other(0,3).imag();
            double k2 = X_other(1,2).real();
            double l2 = X_other(1,2).imag();
            double m2 = X_other(1,3).real();
            double n2 = X_other(1,3).imag();
            double o2 = X_other(2,3).real();
            double p2 = X_other(2,3).imag();

            double zero_zero_entry = -2*a1*d2*d2 - 2*a1*e2*e2 - 2*a1*f2*f2 - 2*a1*g2*g2 - 4*a1*h2*h2 - 4*a1*j2*j2 + 2*a2*d1*d2 + 2*a2*e1*e2 + 2*a2*f1*f2 + 2*a2*g1*g2 + 4*a2*h1*h2 + 4*a2*j1*j2 + 2*b1*d2*d2 + 2*b1*e2*e2 - 2*b1*h2*h2 - 2*b1*j2*j2 - 2*b2*d1*d2 - 2*b2*e1*e2 + 2*b2*h1*h2 + 2*b2*j1*j2 + 2*c1*f2*f2 + 2*c1*g2*g2 - 2*c1*h2*h2 - 2*c1*j2*j2 - 2*c2*f1*f2 - 2*c2*g1*g2 + 2*c2*h1*h2 + 2*c2*j1*j2 - 2*d1*f2*k2 - 2*d1*g2*l2 - 2*d1*h2*m2 - 2*d1*j2*n2 - 2*d2*f1*k2 + 4*d2*f2*k1 - 2*d2*g1*l2 + 4*d2*g2*l1 - 2*d2*h1*m2 + 4*d2*h2*m1 - 2*d2*j1*n2 + 4*d2*j2*n1 + 2*e1*f2*l2 - 2*e1*g2*k2 + 2*e1*h2*n2 - 2*e1*j2*m2 + 2*e2*f1*l2 - 4*e2*f2*l1 - 2*e2*g1*k2 + 4*e2*g2*k1 + 2*e2*h1*n2 - 4*e2*h2*n1 - 2*e2*j1*m2 + 4*e2*j2*m1 - 2*f1*h2*o2 - 2*f1*j2*p2 - 2*f2*h1*o2 + 4*f2*h2*o1 - 2*f2*j1*p2 + 4*f2*j2*p1 + 2*g1*h2*p2 - 2*g1*j2*o2 + 2*g2*h1*p2 - 4*g2*h2*p1 - 2*g2*j1*o2 + 4*g2*j2*o1;
            R_or_C zero_one_entry( a1*a2*d2 - a1*b2*d2 - a1*f2*k2 - a1*g2*l2 - 3*a1*h2*m2 - 3*a1*j2*n2 - a2*a2*d1 - a2*b1*d2 + 2*a2*b2*d1 + 2*a2*f1*k2 - a2*f2*k1 + 2*a2*g1*l2 - a2*g2*l1 + 3*a2*h1*m2 + 3*a2*j1*n2 + b1*b2*d2 - b1*f2*k2 - b1*g2*l2 - 3*b1*h2*m2 - 3*b1*j2*n2 - b2*b2*d1 - b2*f1*k2 + 2*b2*f2*k1 - b2*g1*l2 + 2*b2*g2*l1 + 3*b2*h2*m1 + 3*b2*j2*n1 + 2*c1*f2*k2 + 2*c1*g2*l2 - 2*c1*h2*m2 - 2*c1*j2*n2 - c2*f1*k2 - c2*f2*k1 - c2*g1*l2 - c2*g2*l1 + c2*h1*m2 + c2*h2*m1 + c2*j1*n2 + c2*j2*n1 - 4*d1*e2*e2 - d1*f2*f2 - d1*g2*g2 - d1*h2*h2 - d1*j2*j2 - d1*k2*k2 - d1*l2*l2 - d1*m2*m2 - d1*n2*n2 + 4*d2*e1*e2 + d2*f1*f2 + d2*g1*g2 + d2*h1*h2 + d2*j1*j2 + d2*k1*k2 + d2*l1*l2 + d2*m1*m2 + d2*n1*n2 - 3*e2*f1*g2 + 3*e2*f2*g1 - 3*e2*h1*j2 + 3*e2*h2*j1 + 3*e2*k1*l2 - 3*e2*k2*l1 + 3*e2*m1*n2 - 3*e2*m2*n1 - f1*m2*o2 - f1*n2*p2 - f2*m1*o2 + 2*f2*m2*o1 - f2*n1*p2 + 2*f2*n2*p1 + g1*m2*p2 - g1*n2*o2 + g2*m1*p2 - 2*g2*m2*p1 - g2*n1*o2 + 2*g2*n2*o1 - h1*k2*o2 + h1*l2*p2 - h2*k1*o2 + 2*h2*k2*o1 + h2*l1*p2 - 2*h2*l2*p1 - j1*k2*p2 - j1*l2*o2 - j2*k1*p2 + 2*j2*k2*p1 - j2*l1*o2 + 2*j2*l2*o1 , a1*a2*e2 - a1*b2*e2 + a1*f2*l2 - a1*g2*k2 + 3*a1*h2*n2 - 3*a1*j2*m2 - a2*a2*e1 - a2*b1*e2 + 2*a2*b2*e1 - 2*a2*f1*l2 + a2*f2*l1 + 2*a2*g1*k2 - a2*g2*k1 - 3*a2*h1*n2 + 3*a2*j1*m2 + b1*b2*e2 + b1*f2*l2 - b1*g2*k2 + 3*b1*h2*n2 - 3*b1*j2*m2 - b2*b2*e1 + b2*f1*l2 - 2*b2*f2*l1 - b2*g1*k2 + 2*b2*g2*k1 - 3*b2*h2*n1 + 3*b2*j2*m1 - 2*c1*f2*l2 + 2*c1*g2*k2 + 2*c1*h2*n2 - 2*c1*j2*m2 + c2*f1*l2 + c2*f2*l1 - c2*g1*k2 - c2*g2*k1 - c2*h1*n2 - c2*h2*n1 + c2*j1*m2 + c2*j2*m1 + 4*d1*d2*e2 - 4*d2*d2*e1 + 3*d2*f1*g2 - 3*d2*f2*g1 + 3*d2*h1*j2 - 3*d2*h2*j1 - 3*d2*k1*l2 + 3*d2*k2*l1 - 3*d2*m1*n2 + 3*d2*m2*n1 - e1*f2*f2 - e1*g2*g2 - e1*h2*h2 - e1*j2*j2 - e1*k2*k2 - e1*l2*l2 - e1*m2*m2 - e1*n2*n2 + e2*f1*f2 + e2*g1*g2 + e2*h1*h2 + e2*j1*j2 + e2*k1*k2 + e2*l1*l2 + e2*m1*m2 + e2*n1*n2 - f1*m2*p2 + f1*n2*o2 - f2*m1*p2 + 2*f2*m2*p1 + f2*n1*o2 - 2*f2*n2*o1 - g1*m2*o2 - g1*n2*p2 - g2*m1*o2 + 2*g2*m2*o1 - g2*n1*p2 + 2*g2*n2*p1 + h1*k2*p2 + h1*l2*o2 + h2*k1*p2 - 2*h2*k2*p1 + h2*l1*o2 - 2*h2*l2*o1 - j1*k2*o2 + j1*l2*p2 - j2*k1*o2 + 2*j2*k2*o1 + j2*l1*p2 - 2*j2*l2*p1 );
            R_or_C zero_two_entry( a1*a2*f2 - a1*c2*f2 - a1*d2*k2 + a1*e2*l2 - 3*a1*h2*o2 - 3*a1*j2*p2 - a2*a2*f1 - a2*c1*f2 + 2*a2*c2*f1 + 2*a2*d1*k2 - a2*d2*k1 - 2*a2*e1*l2 + a2*e2*l1 + 3*a2*h1*o2 + 3*a2*j1*p2 + 2*b1*d2*k2 - 2*b1*e2*l2 - 2*b1*h2*o2 - 2*b1*j2*p2 - b2*d1*k2 - b2*d2*k1 + b2*e1*l2 + b2*e2*l1 + b2*h1*o2 + b2*h2*o1 + b2*j1*p2 + b2*j2*p1 + c1*c2*f2 - c1*d2*k2 + c1*e2*l2 - 3*c1*h2*o2 - 3*c1*j2*p2 - c2*c2*f1 - c2*d1*k2 + 2*c2*d2*k1 + c2*e1*l2 - 2*c2*e2*l1 + 3*c2*h2*o1 + 3*c2*j2*p1 + d1*d2*f2 - 3*d1*e2*g2 - d1*m2*o2 - d1*n2*p2 - d2*d2*f1 + 3*d2*e1*g2 + 2*d2*m1*o2 - d2*m2*o1 + 2*d2*n1*p2 - d2*n2*p1 + e1*e2*f2 - e1*m2*p2 + e1*n2*o2 - e2*e2*f1 + 2*e2*m1*p2 - e2*m2*p1 - 2*e2*n1*o2 + e2*n2*o1 - 4*f1*g2*g2 - f1*h2*h2 - f1*j2*j2 - f1*k2*k2 - f1*l2*l2 - f1*o2*o2 - f1*p2*p2 + 4*f2*g1*g2 + f2*h1*h2 + f2*j1*j2 + f2*k1*k2 + f2*l1*l2 + f2*o1*o2 + f2*p1*p2 - 3*g2*h1*j2 + 3*g2*h2*j1 - 3*g2*k1*l2 + 3*g2*k2*l1 + 3*g2*o1*p2 - 3*g2*o2*p1 - h1*k2*m2 - h1*l2*n2 - h2*k1*m2 + 2*h2*k2*m1 - h2*l1*n2 + 2*h2*l2*n1 - j1*k2*n2 + j1*l2*m2 - j2*k1*n2 + 2*j2*k2*n1 + j2*l1*m2 - 2*j2*l2*m1 , a1*a2*g2 - a1*c2*g2 - a1*d2*l2 - a1*e2*k2 + 3*a1*h2*p2 - 3*a1*j2*o2 - a2*a2*g1 - a2*c1*g2 + 2*a2*c2*g1 + 2*a2*d1*l2 - a2*d2*l1 + 2*a2*e1*k2 - a2*e2*k1 - 3*a2*h1*p2 + 3*a2*j1*o2 + 2*b1*d2*l2 + 2*b1*e2*k2 + 2*b1*h2*p2 - 2*b1*j2*o2 - b2*d1*l2 - b2*d2*l1 - b2*e1*k2 - b2*e2*k1 - b2*h1*p2 - b2*h2*p1 + b2*j1*o2 + b2*j2*o1 + c1*c2*g2 - c1*d2*l2 - c1*e2*k2 + 3*c1*h2*p2 - 3*c1*j2*o2 - c2*c2*g1 - c2*d1*l2 + 2*c2*d2*l1 - c2*e1*k2 + 2*c2*e2*k1 - 3*c2*h2*p1 + 3*c2*j2*o1 + d1*d2*g2 + 3*d1*e2*f2 + d1*m2*p2 - d1*n2*o2 - d2*d2*g1 - 3*d2*e1*f2 - 2*d2*m1*p2 + d2*m2*p1 + 2*d2*n1*o2 - d2*n2*o1 + e1*e2*g2 - e1*m2*o2 - e1*n2*p2 - e2*e2*g1 + 2*e2*m1*o2 - e2*m2*o1 + 2*e2*n1*p2 - e2*n2*p1 + 4*f1*f2*g2 - 4*f2*f2*g1 + 3*f2*h1*j2 - 3*f2*h2*j1 + 3*f2*k1*l2 - 3*f2*k2*l1 - 3*f2*o1*p2 + 3*f2*o2*p1 - g1*h2*h2 - g1*j2*j2 - g1*k2*k2 - g1*l2*l2 - g1*o2*o2 - g1*p2*p2 + g2*h1*h2 + g2*j1*j2 + g2*k1*k2 + g2*l1*l2 + g2*o1*o2 + g2*p1*p2 + h1*k2*n2 - h1*l2*m2 + h2*k1*n2 - 2*h2*k2*n1 - h2*l1*m2 + 2*h2*l2*m1 - j1*k2*m2 - j1*l2*n2 - j2*k1*m2 + 2*j2*k2*m1 - j2*l1*n2 + 2*j2*l2*n1 );
            R_or_C zero_three_entry( 4*a1*a2*h2 + 2*a1*b2*h2 + 2*a1*c2*h2 - 4*a2*a2*h1 + 2*a2*b1*h2 - 4*a2*b2*h1 + 2*a2*c1*h2 - 4*a2*c2*h1 + 3*a2*d1*m2 - 3*a2*d2*m1 - 3*a2*e1*n2 + 3*a2*e2*n1 + 3*a2*f1*o2 - 3*a2*f2*o1 - 3*a2*g1*p2 + 3*a2*g2*p1 + b1*b2*h2 + b1*c2*h2 + 3*b1*d2*m2 - 3*b1*e2*n2 + b1*f2*o2 - b1*g2*p2 - b2*b2*h1 + b2*c1*h2 - 2*b2*c2*h1 - 3*b2*d2*m1 + 3*b2*e2*n1 + b2*f1*o2 - 2*b2*f2*o1 - b2*g1*p2 + 2*b2*g2*p1 + c1*c2*h2 + c1*d2*m2 - c1*e2*n2 + 3*c1*f2*o2 - 3*c1*g2*p2 - c2*c2*h1 + c2*d1*m2 - 2*c2*d2*m1 - c2*e1*n2 + 2*c2*e2*n1 - 3*c2*f2*o1 + 3*c2*g2*p1 + d1*d2*h2 - 3*d1*e2*j2 - d1*k2*o2 + d1*l2*p2 - d2*d2*h1 + 3*d2*e1*j2 + 2*d2*k1*o2 - d2*k2*o1 - 2*d2*l1*p2 + d2*l2*p1 + e1*e2*h2 + e1*k2*p2 + e1*l2*o2 - e2*e2*h1 - 2*e2*k1*p2 + e2*k2*p1 - 2*e2*l1*o2 + e2*l2*o1 + f1*f2*h2 - 3*f1*g2*j2 - f1*k2*m2 - f1*l2*n2 - f2*f2*h1 + 3*f2*g1*j2 + 2*f2*k1*m2 - f2*k2*m1 + 2*f2*l1*n2 - f2*l2*n1 + g1*g2*h2 + g1*k2*n2 - g1*l2*m2 - g2*g2*h1 - 2*g2*k1*n2 + g2*k2*n1 + 2*g2*l1*m2 - g2*l2*m1 - 4*h1*j2*j2 - h1*m2*m2 - h1*n2*n2 - h1*o2*o2 - h1*p2*p2 + 4*h2*j1*j2 + h2*m1*m2 + h2*n1*n2 + h2*o1*o2 + h2*p1*p2 - 3*j2*m1*n2 + 3*j2*m2*n1 - 3*j2*o1*p2 + 3*j2*o2*p1 , 4*a1*a2*j2 + 2*a1*b2*j2 + 2*a1*c2*j2 - 4*a2*a2*j1 + 2*a2*b1*j2 - 4*a2*b2*j1 + 2*a2*c1*j2 - 4*a2*c2*j1 + 3*a2*d1*n2 - 3*a2*d2*n1 + 3*a2*e1*m2 - 3*a2*e2*m1 + 3*a2*f1*p2 - 3*a2*f2*p1 + 3*a2*g1*o2 - 3*a2*g2*o1 + b1*b2*j2 + b1*c2*j2 + 3*b1*d2*n2 + 3*b1*e2*m2 + b1*f2*p2 + b1*g2*o2 - b2*b2*j1 + b2*c1*j2 - 2*b2*c2*j1 - 3*b2*d2*n1 - 3*b2*e2*m1 + b2*f1*p2 - 2*b2*f2*p1 + b2*g1*o2 - 2*b2*g2*o1 + c1*c2*j2 + c1*d2*n2 + c1*e2*m2 + 3*c1*f2*p2 + 3*c1*g2*o2 - c2*c2*j1 + c2*d1*n2 - 2*c2*d2*n1 + c2*e1*m2 - 2*c2*e2*m1 - 3*c2*f2*p1 - 3*c2*g2*o1 + d1*d2*j2 + 3*d1*e2*h2 - d1*k2*p2 - d1*l2*o2 - d2*d2*j1 - 3*d2*e1*h2 + 2*d2*k1*p2 - d2*k2*p1 + 2*d2*l1*o2 - d2*l2*o1 + e1*e2*j2 - e1*k2*o2 + e1*l2*p2 - e2*e2*j1 + 2*e2*k1*o2 - e2*k2*o1 - 2*e2*l1*p2 + e2*l2*p1 + f1*f2*j2 + 3*f1*g2*h2 - f1*k2*n2 + f1*l2*m2 - f2*f2*j1 - 3*f2*g1*h2 + 2*f2*k1*n2 - f2*k2*n1 - 2*f2*l1*m2 + f2*l2*m1 + g1*g2*j2 - g1*k2*m2 - g1*l2*n2 - g2*g2*j1 + 2*g2*k1*m2 - g2*k2*m1 + 2*g2*l1*n2 - g2*l2*n1 + 4*h1*h2*j2 - 4*h2*h2*j1 + 3*h2*m1*n2 - 3*h2*m2*n1 + 3*h2*o1*p2 - 3*h2*o2*p1 - j1*m2*m2 - j1*n2*n2 - j1*o2*o2 - j1*p2*p2 + j2*m1*m2 + j2*n1*n2 + j2*o1*o2 + j2*p1*p2 );
            double one_one_entry = 2*a1*d2*d2 + 2*a1*e2*e2 - 2*a1*m2*m2 - 2*a1*n2*n2 - 2*a2*d1*d2 - 2*a2*e1*e2 + 2*a2*m1*m2 + 2*a2*n1*n2 - 2*b1*d2*d2 - 2*b1*e2*e2 - 2*b1*k2*k2 - 2*b1*l2*l2 - 4*b1*m2*m2 - 4*b1*n2*n2 + 2*b2*d1*d2 + 2*b2*e1*e2 + 2*b2*k1*k2 + 2*b2*l1*l2 + 4*b2*m1*m2 + 4*b2*n1*n2 + 2*c1*k2*k2 + 2*c1*l2*l2 - 2*c1*m2*m2 - 2*c1*n2*n2 - 2*c2*k1*k2 - 2*c2*l1*l2 + 2*c2*m1*m2 + 2*c2*n1*n2 - 2*d1*f2*k2 - 2*d1*g2*l2 - 2*d1*h2*m2 - 2*d1*j2*n2 + 4*d2*f1*k2 - 2*d2*f2*k1 + 4*d2*g1*l2 - 2*d2*g2*l1 + 4*d2*h1*m2 - 2*d2*h2*m1 + 4*d2*j1*n2 - 2*d2*j2*n1 + 2*e1*f2*l2 - 2*e1*g2*k2 + 2*e1*h2*n2 - 2*e1*j2*m2 - 4*e2*f1*l2 + 2*e2*f2*l1 + 4*e2*g1*k2 - 2*e2*g2*k1 - 4*e2*h1*n2 + 2*e2*h2*n1 + 4*e2*j1*m2 - 2*e2*j2*m1 - 2*k1*m2*o2 - 2*k1*n2*p2 - 2*k2*m1*o2 + 4*k2*m2*o1 - 2*k2*n1*p2 + 4*k2*n2*p1 + 2*l1*m2*p2 - 2*l1*n2*o2 + 2*l2*m1*p2 - 4*l2*m2*p1 - 2*l2*n1*o2 + 4*l2*n2*o1;
            R_or_C one_two_entry( 2*a1*d2*f2 + 2*a1*e2*g2 - 2*a1*m2*o2 - 2*a1*n2*p2 - a2*d1*f2 - a2*d2*f1 - a2*e1*g2 - a2*e2*g1 + a2*m1*o2 + a2*m2*o1 + a2*n1*p2 + a2*n2*p1 + b1*b2*k2 - b1*c2*k2 - b1*d2*f2 - b1*e2*g2 - 3*b1*m2*o2 - 3*b1*n2*p2 - b2*b2*k1 - b2*c1*k2 + 2*b2*c2*k1 + 2*b2*d1*f2 - b2*d2*f1 + 2*b2*e1*g2 - b2*e2*g1 + 3*b2*m1*o2 + 3*b2*n1*p2 + c1*c2*k2 - c1*d2*f2 - c1*e2*g2 - 3*c1*m2*o2 - 3*c1*n2*p2 - c2*c2*k1 - c2*d1*f2 + 2*c2*d2*f1 - c2*e1*g2 + 2*c2*e2*g1 + 3*c2*m2*o1 + 3*c2*n2*p1 + d1*d2*k2 + 3*d1*e2*l2 - d1*h2*o2 - d1*j2*p2 - d2*d2*k1 - 3*d2*e1*l2 + 2*d2*h1*o2 - d2*h2*o1 + 2*d2*j1*p2 - d2*j2*p1 + e1*e2*k2 + e1*h2*p2 - e1*j2*o2 - e2*e2*k1 - 2*e2*h1*p2 + e2*h2*p1 + 2*e2*j1*o2 - e2*j2*o1 + f1*f2*k2 - 3*f1*g2*l2 - f1*h2*m2 - f1*j2*n2 - f2*f2*k1 + 3*f2*g1*l2 + 2*f2*h1*m2 - f2*h2*m1 + 2*f2*j1*n2 - f2*j2*n1 + g1*g2*k2 + g1*h2*n2 - g1*j2*m2 - g2*g2*k1 - 2*g2*h1*n2 + g2*h2*n1 + 2*g2*j1*m2 - g2*j2*m1 - 4*k1*l2*l2 - k1*m2*m2 - k1*n2*n2 - k1*o2*o2 - k1*p2*p2 + 4*k2*l1*l2 + k2*m1*m2 + k2*n1*n2 + k2*o1*o2 + k2*p1*p2 - 3*l2*m1*n2 + 3*l2*m2*n1 + 3*l2*o1*p2 - 3*l2*o2*p1 , 2*a1*d2*g2 - 2*a1*e2*f2 + 2*a1*m2*p2 - 2*a1*n2*o2 - a2*d1*g2 - a2*d2*g1 + a2*e1*f2 + a2*e2*f1 - a2*m1*p2 - a2*m2*p1 + a2*n1*o2 + a2*n2*o1 + b1*b2*l2 - b1*c2*l2 - b1*d2*g2 + b1*e2*f2 + 3*b1*m2*p2 - 3*b1*n2*o2 - b2*b2*l1 - b2*c1*l2 + 2*b2*c2*l1 + 2*b2*d1*g2 - b2*d2*g1 - 2*b2*e1*f2 + b2*e2*f1 - 3*b2*m1*p2 + 3*b2*n1*o2 + c1*c2*l2 - c1*d2*g2 + c1*e2*f2 + 3*c1*m2*p2 - 3*c1*n2*o2 - c2*c2*l1 - c2*d1*g2 + 2*c2*d2*g1 + c2*e1*f2 - 2*c2*e2*f1 - 3*c2*m2*p1 + 3*c2*n2*o1 + d1*d2*l2 - 3*d1*e2*k2 + d1*h2*p2 - d1*j2*o2 - d2*d2*l1 + 3*d2*e1*k2 - 2*d2*h1*p2 + d2*h2*p1 + 2*d2*j1*o2 - d2*j2*o1 + e1*e2*l2 + e1*h2*o2 + e1*j2*p2 - e2*e2*l1 - 2*e2*h1*o2 + e2*h2*o1 - 2*e2*j1*p2 + e2*j2*p1 + f1*f2*l2 + 3*f1*g2*k2 - f1*h2*n2 + f1*j2*m2 - f2*f2*l1 - 3*f2*g1*k2 + 2*f2*h1*n2 - f2*h2*n1 - 2*f2*j1*m2 + f2*j2*m1 + g1*g2*l2 - g1*h2*m2 - g1*j2*n2 - g2*g2*l1 + 2*g2*h1*m2 - g2*h2*m1 + 2*g2*j1*n2 - g2*j2*n1 + 4*k1*k2*l2 - 4*k2*k2*l1 + 3*k2*m1*n2 - 3*k2*m2*n1 - 3*k2*o1*p2 + 3*k2*o2*p1 - l1*m2*m2 - l1*n2*n2 - l1*o2*o2 - l1*p2*p2 + l2*m1*m2 + l2*n1*n2 + l2*o1*o2 + l2*p1*p2 );
            R_or_C one_three_entry( a1*a2*m2 + 2*a1*b2*m2 + a1*c2*m2 + 3*a1*d2*h2 + 3*a1*e2*j2 + a1*k2*o2 - a1*l2*p2 - a2*a2*m1 + 2*a2*b1*m2 - 4*a2*b2*m1 + a2*c1*m2 - 2*a2*c2*m1 - 3*a2*d2*h1 - 3*a2*e2*j1 + a2*k1*o2 - 2*a2*k2*o1 - a2*l1*p2 + 2*a2*l2*p1 + 4*b1*b2*m2 + 2*b1*c2*m2 - 4*b2*b2*m1 + 2*b2*c1*m2 - 4*b2*c2*m1 + 3*b2*d1*h2 - 3*b2*d2*h1 + 3*b2*e1*j2 - 3*b2*e2*j1 + 3*b2*k1*o2 - 3*b2*k2*o1 - 3*b2*l1*p2 + 3*b2*l2*p1 + c1*c2*m2 + c1*d2*h2 + c1*e2*j2 + 3*c1*k2*o2 - 3*c1*l2*p2 - c2*c2*m1 + c2*d1*h2 - 2*c2*d2*h1 + c2*e1*j2 - 2*c2*e2*j1 - 3*c2*k2*o1 + 3*c2*l2*p1 + d1*d2*m2 + 3*d1*e2*n2 - d1*f2*o2 + d1*g2*p2 - d2*d2*m1 - 3*d2*e1*n2 + 2*d2*f1*o2 - d2*f2*o1 - 2*d2*g1*p2 + d2*g2*p1 + e1*e2*m2 - e1*f2*p2 - e1*g2*o2 - e2*e2*m1 + 2*e2*f1*p2 - e2*f2*p1 + 2*e2*g1*o2 - e2*g2*o1 + 2*f1*h2*k2 - 2*f1*j2*l2 - f2*h1*k2 - f2*h2*k1 + f2*j1*l2 + f2*j2*l1 + 2*g1*h2*l2 + 2*g1*j2*k2 - g2*h1*l2 - g2*h2*l1 - g2*j1*k2 - g2*j2*k1 + h1*h2*m2 - 3*h1*j2*n2 - h2*h2*m1 + 3*h2*j1*n2 + j1*j2*m2 - j2*j2*m1 + k1*k2*m2 - 3*k1*l2*n2 - k2*k2*m1 + 3*k2*l1*n2 + l1*l2*m2 - l2*l2*m1 - 4*m1*n2*n2 - m1*o2*o2 - m1*p2*p2 + 4*m2*n1*n2 + m2*o1*o2 + m2*p1*p2 - 3*n2*o1*p2 + 3*n2*o2*p1 , a1*a2*n2 + 2*a1*b2*n2 + a1*c2*n2 + 3*a1*d2*j2 - 3*a1*e2*h2 + a1*k2*p2 + a1*l2*o2 - a2*a2*n1 + 2*a2*b1*n2 - 4*a2*b2*n1 + a2*c1*n2 - 2*a2*c2*n1 - 3*a2*d2*j1 + 3*a2*e2*h1 + a2*k1*p2 - 2*a2*k2*p1 + a2*l1*o2 - 2*a2*l2*o1 + 4*b1*b2*n2 + 2*b1*c2*n2 - 4*b2*b2*n1 + 2*b2*c1*n2 - 4*b2*c2*n1 + 3*b2*d1*j2 - 3*b2*d2*j1 - 3*b2*e1*h2 + 3*b2*e2*h1 + 3*b2*k1*p2 - 3*b2*k2*p1 + 3*b2*l1*o2 - 3*b2*l2*o1 + c1*c2*n2 + c1*d2*j2 - c1*e2*h2 + 3*c1*k2*p2 + 3*c1*l2*o2 - c2*c2*n1 + c2*d1*j2 - 2*c2*d2*j1 - c2*e1*h2 + 2*c2*e2*h1 - 3*c2*k2*p1 - 3*c2*l2*o1 + d1*d2*n2 - 3*d1*e2*m2 - d1*f2*p2 - d1*g2*o2 - d2*d2*n1 + 3*d2*e1*m2 + 2*d2*f1*p2 - d2*f2*p1 + 2*d2*g1*o2 - d2*g2*o1 + e1*e2*n2 + e1*f2*o2 - e1*g2*p2 - e2*e2*n1 - 2*e2*f1*o2 + e2*f2*o1 + 2*e2*g1*p2 - e2*g2*p1 + 2*f1*h2*l2 + 2*f1*j2*k2 - f2*h1*l2 - f2*h2*l1 - f2*j1*k2 - f2*j2*k1 - 2*g1*h2*k2 + 2*g1*j2*l2 + g2*h1*k2 + g2*h2*k1 - g2*j1*l2 - g2*j2*l1 + h1*h2*n2 + 3*h1*j2*m2 - h2*h2*n1 - 3*h2*j1*m2 + j1*j2*n2 - j2*j2*n1 + k1*k2*n2 + 3*k1*l2*m2 - k2*k2*n1 - 3*k2*l1*m2 + l1*l2*n2 - l2*l2*n1 + 4*m1*m2*n2 - 4*m2*m2*n1 + 3*m2*o1*p2 - 3*m2*o2*p1 - n1*o2*o2 - n1*p2*p2 + n2*o1*o2 + n2*p1*p2 );
            double two_two_entry = 2*a1*f2*f2 + 2*a1*g2*g2 - 2*a1*o2*o2 - 2*a1*p2*p2 - 2*a2*f1*f2 - 2*a2*g1*g2 + 2*a2*o1*o2 + 2*a2*p1*p2 + 2*b1*k2*k2 + 2*b1*l2*l2 - 2*b1*o2*o2 - 2*b1*p2*p2 - 2*b2*k1*k2 - 2*b2*l1*l2 + 2*b2*o1*o2 + 2*b2*p1*p2 - 2*c1*f2*f2 - 2*c1*g2*g2 - 2*c1*k2*k2 - 2*c1*l2*l2 - 4*c1*o2*o2 - 4*c1*p2*p2 + 2*c2*f1*f2 + 2*c2*g1*g2 + 2*c2*k1*k2 + 2*c2*l1*l2 + 4*c2*o1*o2 + 4*c2*p1*p2 + 4*d1*f2*k2 + 4*d1*g2*l2 - 2*d2*f1*k2 - 2*d2*f2*k1 - 2*d2*g1*l2 - 2*d2*g2*l1 - 4*e1*f2*l2 + 4*e1*g2*k2 + 2*e2*f1*l2 + 2*e2*f2*l1 - 2*e2*g1*k2 - 2*e2*g2*k1 - 2*f1*h2*o2 - 2*f1*j2*p2 + 4*f2*h1*o2 - 2*f2*h2*o1 + 4*f2*j1*p2 - 2*f2*j2*p1 + 2*g1*h2*p2 - 2*g1*j2*o2 - 4*g2*h1*p2 + 2*g2*h2*p1 + 4*g2*j1*o2 - 2*g2*j2*o1 - 2*k1*m2*o2 - 2*k1*n2*p2 + 4*k2*m1*o2 - 2*k2*m2*o1 + 4*k2*n1*p2 - 2*k2*n2*p1 + 2*l1*m2*p2 - 2*l1*n2*o2 - 4*l2*m1*p2 + 2*l2*m2*p1 + 4*l2*n1*o2 - 2*l2*n2*o1;
            R_or_C two_three_entry( a1*a2*o2 + a1*b2*o2 + 2*a1*c2*o2 + 3*a1*f2*h2 + 3*a1*g2*j2 + a1*k2*m2 + a1*l2*n2 - a2*a2*o1 + a2*b1*o2 - 2*a2*b2*o1 + 2*a2*c1*o2 - 4*a2*c2*o1 - 3*a2*f2*h1 - 3*a2*g2*j1 + a2*k1*m2 - 2*a2*k2*m1 + a2*l1*n2 - 2*a2*l2*n1 + b1*b2*o2 + 2*b1*c2*o2 + b1*f2*h2 + b1*g2*j2 + 3*b1*k2*m2 + 3*b1*l2*n2 - b2*b2*o1 + 2*b2*c1*o2 - 4*b2*c2*o1 + b2*f1*h2 - 2*b2*f2*h1 + b2*g1*j2 - 2*b2*g2*j1 - 3*b2*k2*m1 - 3*b2*l2*n1 + 4*c1*c2*o2 - 4*c2*c2*o1 + 3*c2*f1*h2 - 3*c2*f2*h1 + 3*c2*g1*j2 - 3*c2*g2*j1 + 3*c2*k1*m2 - 3*c2*k2*m1 + 3*c2*l1*n2 - 3*c2*l2*n1 + 2*d1*f2*m2 + 2*d1*g2*n2 + 2*d1*h2*k2 + 2*d1*j2*l2 - d2*f1*m2 - d2*f2*m1 - d2*g1*n2 - d2*g2*n1 - d2*h1*k2 - d2*h2*k1 - d2*j1*l2 - d2*j2*l1 - 2*e1*f2*n2 + 2*e1*g2*m2 - 2*e1*h2*l2 + 2*e1*j2*k2 + e2*f1*n2 + e2*f2*n1 - e2*g1*m2 - e2*g2*m1 + e2*h1*l2 + e2*h2*l1 - e2*j1*k2 - e2*j2*k1 + f1*f2*o2 + 3*f1*g2*p2 - f2*f2*o1 - 3*f2*g1*p2 + g1*g2*o2 - g2*g2*o1 + h1*h2*o2 - 3*h1*j2*p2 - h2*h2*o1 + 3*h2*j1*p2 + j1*j2*o2 - j2*j2*o1 + k1*k2*o2 + 3*k1*l2*p2 - k2*k2*o1 - 3*k2*l1*p2 + l1*l2*o2 - l2*l2*o1 + m1*m2*o2 - 3*m1*n2*p2 - m2*m2*o1 + 3*m2*n1*p2 + n1*n2*o2 - n2*n2*o1 - 4*o1*p2*p2 + 4*o2*p1*p2 , a1*a2*p2 + a1*b2*p2 + 2*a1*c2*p2 + 3*a1*f2*j2 - 3*a1*g2*h2 + a1*k2*n2 - a1*l2*m2 - a2*a2*p1 + a2*b1*p2 - 2*a2*b2*p1 + 2*a2*c1*p2 - 4*a2*c2*p1 - 3*a2*f2*j1 + 3*a2*g2*h1 + a2*k1*n2 - 2*a2*k2*n1 - a2*l1*m2 + 2*a2*l2*m1 + b1*b2*p2 + 2*b1*c2*p2 + b1*f2*j2 - b1*g2*h2 + 3*b1*k2*n2 - 3*b1*l2*m2 - b2*b2*p1 + 2*b2*c1*p2 - 4*b2*c2*p1 + b2*f1*j2 - 2*b2*f2*j1 - b2*g1*h2 + 2*b2*g2*h1 - 3*b2*k2*n1 + 3*b2*l2*m1 + 4*c1*c2*p2 - 4*c2*c2*p1 + 3*c2*f1*j2 - 3*c2*f2*j1 - 3*c2*g1*h2 + 3*c2*g2*h1 + 3*c2*k1*n2 - 3*c2*k2*n1 - 3*c2*l1*m2 + 3*c2*l2*m1 + 2*d1*f2*n2 - 2*d1*g2*m2 - 2*d1*h2*l2 + 2*d1*j2*k2 - d2*f1*n2 - d2*f2*n1 + d2*g1*m2 + d2*g2*m1 + d2*h1*l2 + d2*h2*l1 - d2*j1*k2 - d2*j2*k1 + 2*e1*f2*m2 + 2*e1*g2*n2 - 2*e1*h2*k2 - 2*e1*j2*l2 - e2*f1*m2 - e2*f2*m1 - e2*g1*n2 - e2*g2*n1 + e2*h1*k2 + e2*h2*k1 + e2*j1*l2 + e2*j2*l1 + f1*f2*p2 - 3*f1*g2*o2 - f2*f2*p1 + 3*f2*g1*o2 + g1*g2*p2 - g2*g2*p1 + h1*h2*p2 + 3*h1*j2*o2 - h2*h2*p1 - 3*h2*j1*o2 + j1*j2*p2 - j2*j2*p1 + k1*k2*p2 - 3*k1*l2*o2 - k2*k2*p1 + 3*k2*l1*o2 + l1*l2*p2 - l2*l2*p1 + m1*m2*p2 + 3*m1*n2*o2 - m2*m2*p1 - 3*m2*n1*o2 + n1*n2*p2 - n2*n2*p1 + 4*o1*o2*p2 - 4*o2*o2*p1 );

            temp_commutator(0,0) = zero_zero_entry;
            temp_commutator(0,1) = zero_one_entry;
            temp_commutator(0,2) = zero_two_entry; 
            temp_commutator(0,3) = zero_three_entry;  
            temp_commutator(1,0) = std::conj(zero_one_entry);  
            temp_commutator(1,1) = one_one_entry;  
            temp_commutator(1,2) = one_two_entry;  
            temp_commutator(1,3) = one_three_entry;  
            temp_commutator(2,0) = std::conj(zero_two_entry);  
            temp_commutator(2,1) = std::conj(one_two_entry);  
            temp_commutator(2,2) = two_two_entry;  
            temp_commutator(2,3) = two_three_entry;  
            temp_commutator(3,0) = std::conj(zero_three_entry);  
            temp_commutator(3,1) = std::conj(one_three_entry);  
            temp_commutator(3,2) = std::conj(two_three_entry);  
            temp_commutator(3,3) = - zero_zero_entry - one_one_entry - two_two_entry;  

                        
                                    
            commutator_sum += temp_commutator;
        }   
    }
    return commutator_sum;
}

int main() 
{


    // Create  vectors to store the matrices
    matrix X_vector[dim];
    matrix V_vector[dim];
    matrix A_vector[dim];

    // Create a zero matrix in order to populate the V_vector with it.
    matrix zero_matrix = matrix::Zero(rows, cols);


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

/*
    // Generate and store X1, X2, X3, X4, X5, X6, X7, X8, and X9
    for (int i = 0; i < dim; ++i) 
    {
        X_vector[i] = generateHermitianMatrix(rows, cols);
    }
*/



    // Generate and store V1, V2, V3, V4, V5, V6, V7, V8, and V9
    for (int i = 0; i < dim; ++i) 
    {
        V_vector[i] = zero_matrix;
    }  

    // Generate and store A1, A2, A3, A4, A5, A6, A7, A8, and A9
    for (int i = 0; i < dim; ++i) 
    {
        A_vector[i] = Acceleration2(i, X_vector, rows, cols, g);
    }  


    // Create  vectors to store the new matrices
    matrix X_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix V_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix A_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};

    // Write simulation to thermalise system
    for (int j = 0; j < seconds_thermalised / delta_t; ++j)
    {

        // velocity Verlet 1 to get new positions from old positions, momentums and rate of change of momentums

        for (int i = 0; i < 9; ++i)
        {
            X_vector_new[i] = X_vector[i] + V_vector[i] * delta_t + 0.5 * A_vector[i] * delta_t * delta_t;
        }

        // Generate and store new A1, A2, A3, A4, A5, A6, A7, A8, and A9
        for (int i = 0; i < 9; ++i) 
        {
            A_vector_new[i] = Acceleration2( i, X_vector_new, rows, cols, g);
        }  
        
        // Use Velocity Verlet 2 to get new momentums
        for (int i = 0; i < 9; ++i)
        {
            V_vector_new[i] = V_vector[i] + 0.5 * (A_vector_new[i] + A_vector[i]) * delta_t;
        }

        // Copy elements from X_vector_new to X_vector
        std::memcpy(X_vector, X_vector_new, sizeof(X_vector_new));  

        // Copy elements from X_vector_new to X_vector
        std::memcpy(V_vector, V_vector_new, sizeof(V_vector_new)); 

        // Copy elements from X_vector_new to X_vector
        std::memcpy(A_vector, A_vector_new, sizeof(A_vector_new)); 

        if (j % 1000 == 0)
        {
            //for (matrix el : V_vector)
            //{
            //    std::cout <<"Ideal " << el << std::endl;
            //}
            //std::cout  << std::endl << gauss_law(X_vector_new[0], X_vector_new[1], X_vector_new[2], X_vector_new[3], X_vector_new[4], X_vector_new[5], X_vector_new[6], X_vector_new[7], X_vector_new[8],
            //                       V_vector_new[0], V_vector_new[1], V_vector_new[2], V_vector_new[3], V_vector_new[4], V_vector_new[5], V_vector_new[6], V_vector_new[7], V_vector_new[8]);

            std::cout << std::endl;
            std::cout << "H" << std::setprecision(15) << H(1.0, 
                            X_vector_new[0], X_vector_new[1], X_vector_new[2], X_vector_new[3], X_vector_new[4], X_vector_new[5], X_vector_new[6], X_vector_new[7], X_vector_new[8],
                          V_vector_new[0], V_vector_new[1], V_vector_new[2], V_vector_new[3], V_vector_new[4], V_vector_new[5], V_vector_new[6], V_vector_new[7], V_vector_new[8]);


        }

    }




       // Export initial X/V/A_vector to text files to be analysed in python.
    std:: fstream X2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-BranesNis2/thermalised_X.txt", std:: ios:: out);
    X2_vector_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (matrix Matrix : X_vector_new)
    {
        X2_vector_Export << Matrix << std::endl;
    }


    std:: fstream V2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-BranesNis2/thermalised_V.txt", std:: ios:: out);
    V2_vector_Export << std::fixed << std::setprecision(15);
    for (matrix Matrix : V_vector_new)
    {
        V2_vector_Export << Matrix << std::endl;
    }

    std:: fstream A2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-BranesNis2/thermalised_A.txt", std:: ios:: out);
    A2_vector_Export << std::fixed << std::setprecision(15);
    // Print to text file
    for (matrix Matrix : A_vector_new)
    {
        A2_vector_Export << Matrix << std::endl;
    }


    X2_vector_Export.close();
    V2_vector_Export.close();
    A2_vector_Export.close();
 


    return 0;
}