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
const int N = 3;

typedef Eigen::Matrix<double, N, N> matrix;

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
            double element = dist(gen);

            Matrix(i,j) = element;
            if (i != j) 
            {
                Matrix(j,i) = std::conj(element);
            }
        }
    }
    // Make matrix traceless by replacing last entry with the negative of the sum of the other entries.
    // Calculate the sum of diagonal elements in X1
    double diagonalSum = 0.0;
    for (int i = 0; i < rows-1; ++i) 
    {
        diagonalSum += Matrix(i,i);
    }
    Matrix(rows - 1,cols - 1) = -diagonalSum;

    return Matrix;
}

matrix commutator(matrix A, matrix B)
{
    return A * B - B * A;
} 

// Acceleration of each coordinate matrix
matrix Acceleration(const int i, std::vector<matrix> X_vector, int rows, int cols)
{
    matrix commutator_sum = matrix::Zero(rows, cols);
    std::cout << "zeros" << commutator_sum;
    matrix X = X_vector[i];
    for (int j = 0; j < 4; ++j)
    {
        if (i != j)
        {  
            matrix temp_commutator = commutator(X_vector[j], commutator(X, X_vector[j]));
            
            commutator_sum += temp_commutator;
            
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
    matrix zero_matrix = matrix::Zero(rows, cols);
   
    Eigen::Matrix3d test, test1, test2, test3;
    test << 1,2,3,
            4,5,6,
            7,8,9;



    matrix answer2 = Acceleration(0, X_vector, 3, 3);

    for (matrix matrix : X_vector)
    {
        std::cout << std::endl << matrix << std::endl << "New" ;
    }
    

   
    return 0;

}