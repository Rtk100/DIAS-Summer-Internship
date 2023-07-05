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
const int N = 2;

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
matrix Acceleration(const int i, matrix* X_vector, int rows, int cols)
{
    matrix commutator_sum = matrix::Zero(rows, cols);
    std::cout << "zeros" << commutator_sum;
    matrix temp_commutator = matrix::Zero(rows, cols);
    

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
// Acceleration of each coordinate matrix
matrix Acceleration2(const int i, matrix* X_vector, int rows, int cols)
{
    matrix commutator_sum = matrix::Zero(rows, cols);
    std::cout << "zeros" << commutator_sum;
    matrix temp_commutator = matrix::Zero(rows, cols);

    matrix X = X_vector[i];
    for (int j = 0; j < 4; ++j)
    {
        if (i != j)
        {  
            matrix X_other = X_vector[j];
            double g = X_other(0,0) * X(1,1) - X_other(1,1) * X(0,0);
            double b = X_other(0,1);
            double h = X(0,1) * b * g;

            temp_commutator(0,0) = h * 2 * b;
            temp_commutator(0,1) = h * (X_other(0,0) - X_other(1,1));
            temp_commutator(1,0) = temp_commutator(0,1); 
            temp_commutator(1,1) = - temp_commutator(0,0);  
            
            commutator_sum += temp_commutator;
            
        }   
    }
    return commutator_sum;

}

int main() 
{
    for (int i = 0; i < 9; ++i)
    {
        std::cout << i;
    }

    for (int i = 0; i < 9; i++)
    {
        std::cout << i;
    }
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

    matrix test_list2[2];

    test_list2[0](0,0) = 1.0;
    test_list2[0](0,1) = 2.0;
    test_list2[0](1,0) = 3.0;
    test_list2[0](1,1) = 4.0;
    test_list2[1](0,0) = 5.0;
    test_list2[1](0,1) = 6.0;
    test_list2[1](1,0) = 7.0;
    test_list2[1](1,1) = 8.0;

    std::cout << "\n" << "Correct A" <<  Acceleration2(0, test_list2, 2, 2);
    std::cout << "\n" << "Equals A?" << Acceleration(0, test_list2, 2, 2);

   
    return 0;

}