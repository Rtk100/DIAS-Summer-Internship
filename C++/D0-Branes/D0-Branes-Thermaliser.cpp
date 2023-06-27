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

    // Export initial X/V/A_vector to text files to be analysed in python.

    std:: fstream X_vector_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/X_vector_initial.txt", std:: ios:: out);
    X_vector_Export << std::fixed << std::setprecision(15);

    //Print to text file
    for (matrix Matrix : X_vector)
        {
            X_vector_Export << Matrix << ",";
        }

    X_vector_Export.close();


    std:: fstream V_vector_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/V_vector_initial.txt", std:: ios:: out);
    V_vector_Export << std::fixed << std::setprecision(15);

for (matrix Matrix : V_vector)
    {
    
        V_vector_Export << Matrix << ",";
    }

    V_vector_Export.close();

    std:: fstream A_vector_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/A_vector_initial.txt", std:: ios:: out);
    A_vector_Export << std::fixed << std::setprecision(15);
    // Print to text file
    for (matrix Matrix : A_vector)
    {
        A_vector_Export << Matrix << ",";
    }

    A_vector_Export.close();

    // Write simulation to evolve system





    /*
    
    


    std:: fstream ThermalisationExport("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/ThermalisedPart3.txt", std:: ios:: out);
    ThermalisationExport << std::fixed << std::setprecision(15);
    ThermalisationExport << x_1_new << x_2_new << x_3_new << x_4_new << x_5_new << x_6_new << x_7_new << x_8_new << x_9_new;
    ThermalisationExport << v_1_new << v_2_new << v_3_new << v_4_new << v_5_new << v_6_new << v_7_new << v_8_new << v_9_new;
    ThermalisationExport << a_1_new << a_2_new << a_3_new << a_4_new << a_5_new << a_6_new << a_7_new << a_8_new << a_9_new;
    */
    return 0;
}