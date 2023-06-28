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
const double delta_t = 10e-5;
// Repeat simulation for 1000 seconds.
const int simulation_repetitions = 1000 / delta_t;
// Number of D0-Branes
const int N = 3;


typedef double R_or_C;
typedef Eigen:: Matrix<double, N, N> matrix;
/*
To go from real matrices to complex matrices delete the above typedefs and use these typedefs.
And change the commented out code in generateHermitianMatrix()

typedef std::complex<double> R_or_C;
typedef Eigen:: Matrix<std::complex<double>, N, N> matrix;

*/

// Random Hermitian matrix
matrix generateHermitianMatrix(int rows, int cols) 
{
    std::random_device rd;
    std::mt19937 rng(std::time(nullptr));
    std::normal_distribution<double> dist(0.0, 1.0);

    matrix Matrix(rows, cols);

    for (int i = 0; i < rows; ++i) 
    {
        for (int j = 0; j <= i; ++j) 
        {
            double element = dist(rng);
            /* For Complex elements Uncomment this
            double real = dist(rng);
            double imag = dist(rng);
            R_or_C element(real, imag);
            */

            Matrix(i,j) = element;
            if (i != j) 
            {
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

// Commutator of two matrices

matrix commutator(matrix A, matrix B)
{
    return A * B - B * A;
} 

// Original commutator with loops.
/*
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
*/

// Acceleration of each coordinate matrix
matrix Acceleration(const int j, std::vector<matrix> X_vector, int rows, int cols)
{
    matrix commutator_sum = Eigen::MatrixXd::Zero(rows, cols);

    matrix X = X_vector[j];
    for (int i = 0; i < 9; ++i)
    {
        if (i != j)
        {  
            matrix temp_commutator = commutator(X_vector[i], commutator(X, X_vector[i]));
            

            commutator_sum = commutator_sum + temp_commutator;
            
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
        A_vector.push_back( Acceleration(i, X_vector, rows, cols) );
    }  

    // Export initial X/V/A_vector to text files to be analysed in python.

    std:: fstream X_vector_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/X_vector_initial.txt", std:: ios:: out);
    X_vector_Export << std::fixed << std::setprecision(15);

    //Print to text file
    for (matrix Matrix : X_vector)
        {
            X_vector_Export << Matrix << ";";
        }


    std:: fstream V_vector_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/V_vector_initial.txt", std:: ios:: out);
    V_vector_Export << std::fixed << std::setprecision(15);

    for (matrix Matrix : V_vector)
        {
            V_vector_Export << Matrix << ";";
        }

    std:: fstream A_vector_Export("C:/Users/robtk/OneDrive/Desktop/DIAS Internship/Raw data/Harmonic oscillator warm up/A_vector_initial.txt", std:: ios:: out);
    A_vector_Export << std::fixed << std::setprecision(15);
    // Print to text file
    for (matrix Matrix : A_vector)
    {
        A_vector_Export << Matrix << ";";
    }

    // Create  vectors to store the new matrices
    std::vector<matrix> X_vector_new;
    std::vector<matrix> V_vector_new;
    std::vector<matrix> A_vector_new;

    // Make the vectors the correct size
        for (int i = 0; i < 9; ++i) 
    {
        X_vector_new.push_back(zero_matrix);
        V_vector_new.push_back(zero_matrix);
        A_vector_new.push_back(zero_matrix);

    }  

    std::cout << "original" << X_vector[0];
    // Write simulation to thermalise system
    std::cout<< "Number of simulation repetitions:" << 1.0 / delta_t;
    for (int j = 0; j < 10000; ++j)
    {


        // velocity Verlet 1 to get new positions from old positions, momentums and rate of change of momentums

        for (int i = 0; i < 9; ++i)
        {
            X_vector_new[i] = X_vector[i] + V_vector[i] * delta_t + 0.5 * A_vector[i] * delta_t * delta_t;
        }


        std::cout <<std::endl << "X_vector new[0] =" << X_vector[0] + V_vector[0] * delta_t + 0.5 * A_vector[0] * delta_t * delta_t;


        // Generate and store new A1, A2, A3, A4, A5, A6, A7, A8, and A9
        for (int i = 0; i < 9; ++i) 
        {
            A_vector_new[i] = Acceleration( i, X_vector_new, rows, cols);
        }  
        
        // Use Velocity Verlet 2 to get new momentums
        for (int i = 0; i < 9; ++i)
        {
            V_vector_new[i] = V_vector[i] + 0.5 * (A_vector_new[i] + A_vector[i]) * delta_t;
        }

        std::cout << std::endl<< "Before"<< X_vector[0]<< std::endl;
        X_vector = X_vector_new;
        V_vector = V_vector_new;
        A_vector = A_vector_new;
        std::cout << std::endl<< "After"<< X_vector[0]<< std::endl;
        std::cout << std::endl<< "Ideal" << X_vector_new[0] << std::endl;

        if (j % 1000 == 0)
        {
            std::cout << A_vector_new[0];

        }

    }



    for ( matrix Matrix : X_vector_new)
    {
        X_vector_Export << Matrix << ";";
    }

    for ( matrix Matrix : V_vector_new)
    {
        X_vector_Export << Matrix << ";";
    }

    for ( matrix Matrix : A_vector_new)
    {
        A_vector_Export << Matrix << ";";
    }
    X_vector_Export.close();
    V_vector_Export.close();
    A_vector_Export.close();

    return 0;
}