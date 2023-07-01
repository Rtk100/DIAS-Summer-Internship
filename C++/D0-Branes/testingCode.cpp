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
   
    Eigen::Matrix3d test, test1, test2, test3;
    test << 0, 0, 0, 0, 1, 0, 0, 0, 0;


    test1 << 1,0,0,0,0,0,0,0,0;

    test2 << 3,7,7,7,7,7,7,7,6;

    test3 << 3,9,9,9,9,9,99,9,6;

    matrix test_old[2] = {test, test1};
    matrix test_list2[2] = {test2, test3};

for (int i = 0; i < 9; i++)
    {
        std::cout << std::endl << i << std::endl << std::endl;


        std::cout << "Before; 1:" << test_old[0] << std::endl << test_old[1];
        std::cout << std::endl << "Before; 2:" << test_list2[0] << std::endl << test_list2[1];

        test_list2[0] = test_old[0] * 10.0;
        test_list2[1] = test_old[1] * 10.0;

        std::memcpy(test_old, test_list2, sizeof(test_list2));  
/*
    for (int i = 0; i < 9; i++)
    {
        test_list1[i] = test_list2[i]
    }
*/

        std::cout << std::endl << "After; 1:" << test_old[0] << std::endl << test_old[1];
        std::cout << std::endl << "After; 2:" << test_list2[0] << std::endl << test_list2[1];
    }
    

   
    return 0;

}