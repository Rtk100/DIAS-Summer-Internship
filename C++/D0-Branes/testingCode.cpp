#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <random>
#include <complex>
#include <vector>
#include "eigen/Eigen/Dense"


// Define timestep
const double delta_t = 10e-5;
// Repeat simulation for 1000 seconds.
const int simulation_repetitions = 1000 / delta_t;


typedef Eigen::MatrixXcd matrix;

// Random Hermitian matrix
matrix generateHermitianMatrix(int rows, int cols) 
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, 1.0);

    matrix Matrix(rows, cols);

    for (int i = 0; i < rows; ++i) 
    {
        for (int j = 0; j <= i; ++j) 
        {
            double real = dist(gen);
            double imag = dist(gen);
            std::complex<double> element(real, imag);

            Matrix(i,j) = element;
            if (i != j) 
            {
                Matrix(j,i) = std::conj(element);
            }
        }
    }
    // Make matrix traceless by replacing last entry with the negative of the sum of the other entries.
    // Calculate the sum of diagonal elements in X1
    std::complex<double> diagonalSum = 0.0;
    for (int i = 0; i < rows-1; ++i) 
    {
        diagonalSum += Matrix(i,i);
    }
    Matrix(rows - 1,cols - 1) = -diagonalSum;

    return Matrix;
}

// Commutator of two matrices
matrix commutator(const matrix& A, const matrix& B, int rows, int cols) 
{


    matrix result(rows, cols);

    for (int i = 0; i < rows; ++i) 
    {
        for (int j = 0; j < cols; ++j) 
        {
            result(i,j) = 0.0;
            for (int k = 0; k < cols; ++k) 
            {
                result(i,j) += (A(i,k) * B(k,j)) - (B(i,k) * A(k,j));
            }
        }
    }

    return result;
}

// Acceleration of each coordinate matrix
matrix Acceleration(matrix X, std::vector<matrix> X_vector, const double delta_t, matrix commutator_sum, int rows, int cols)
{
    
    for (int i = 0; i < 9; ++i)
    {
         matrix temp_commutator = commutator(X, X_vector[i], rows, cols);
         matrix one_summand = commutator(X_vector[i], temp_commutator, rows, cols);
        
        for (int i = 0; i < rows; ++i) 
        {
            for (int j = 0; j < cols; ++j) 
            {
                commutator_sum(i,j) += one_summand(i,j);
            }
        }
    }   
    return commutator_sum;
}

int main() 
{
    const int rows = 3;
    const int cols = 3;

    // Create  vectors to store the matrices
    std::vector<matrix> X_vector;
    std::vector<matrix> V_vector;
    std::vector<matrix> A_vector;


    // Generate and store X1, X2, X3, X4, X5, X6, X7, X8, and X9
    for (int i = 0; i < 9; ++i) 
    {
        X_vector.push_back(generateHermitianMatrix(rows, cols));
    }

    // Create a zero matrix in order to populate the V_vector with it.
    matrix zero_matrix = Eigen::MatrixXd::Zero(rows, cols);

    // Generate and store V1, V2, V3, V4, V5, V6, V7, V8, and V9
    for (int i = 0; i < 9; ++i) 
    {
        V_vector.push_back(zero_matrix);
    }  

    // Generate and store A1, A2, A3, A4, A5, A6, A7, A8, and A9
    for (int i = 0; i < 9; ++i) 
    {
        A_vector.push_back( Acceleration( X_vector[i], X_vector, delta_t, zero_matrix, rows, cols) );
    }  


    // Write simulation to thermalise system
    std::vector<matrix> X_vector_new;
    std::vector<matrix> V_vector_new;
    std::vector<matrix> A_vector_new;

    for (int i = 0; i < 1 / delta_t; ++i)
    {
        // Create  vectors to store the new matrices

                std::cout << X_vector[0];

        // velocity Verlet 1 to get new positions from old positions, momentums and rate of change of momentums

        for (int i = 0; i < 9; ++i)
        {
            X_vector_new[i] = X_vector[i] + V_vector[i] * delta_t + 0.5 * A_vector[i] * delta_t * delta_t;
        }
        // Generate and store new A1, A2, A3, A4, A5, A6, A7, A8, and A9
        for (int i = 0; i < 9; ++i) 
        {
            A_vector.push_back( Acceleration( X_vector_new[i], X_vector_new, delta_t, zero_matrix, rows, cols) );
        }  
        
        // Use Velocity Verlet 2 to get new momentums
        for (int i = 0; i < 9; ++i)
        {
            V_vector_new[i] = V_vector[i] + 0.5 * (A_vector_new[i] + A_vector[i]) * delta_t;
        }

        X_vector = X_vector_new;
        V_vector = V_vector_new;
        A_vector = A_vector_new;
    }


    return 0;
}