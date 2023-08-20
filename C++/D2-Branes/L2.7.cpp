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

// Define timestep
const long double delta_t = 1e-3;
const long double seconds_thermalised = 10;

// Repeat simulation for 1000 seconds.
const int simulation_repetitions = seconds_thermalised / delta_t;
// Number of D0-Branes
const int N = 3;
const long double g = 1.0/sqrt(N);
const int rows = N;
const int cols = N;

// Dimension of space
const int dim = 9;


//To go from real matrices to complex matrices delete the above typedefs and use these typedefs.
//And change the commented out code in generateHermitianMatrix()
typedef double l_or_d;
typedef std::complex<l_or_d> R_or_C;
typedef Eigen:: Matrix<std::complex<l_or_d>, N, N> matrix;
typedef Eigen:: Matrix<std::complex<l_or_d>, 1, N> row_vector;
typedef Eigen:: Matrix<std::complex<l_or_d>, N, 1> col_vector;



matrix commutator(matrix A, matrix B)
{
    return A * B - B * A;
} 

// Initialising all Xs and Vs and As.
l_or_d X_scalars_vector[27] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
matrix X_matrices_vector[9] = {matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols)};
l_or_d V_scalars_vector[27] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
matrix V_matrices_vector[9] = {matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols)};
l_or_d A_scalars_vector[27] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
matrix A_matrices_vector[9] = {matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols)};


// Initialising all Zs and Z_dots.
// Z_scalars_vector = {Z12, Z13, Z21, Z23, Z31, Z32};
R_or_C C_zero(0.0, 0.0);
R_or_C Z_scalars_vector[6] = {C_zero, C_zero, C_zero, C_zero, C_zero, C_zero};
R_or_C Z_dot_scalars_vector[6] = {C_zero, C_zero, C_zero, C_zero, C_zero, C_zero};

// Z_row_vectors_vector = {Z14, Z24, Z34};
row_vector Z_row_vectors_vector[3] = {matrix::Zero(1, cols), matrix::Zero(1, cols), matrix::Zero(1, cols)};
row_vector Z_dot_row_vectors_vector[3] = {matrix::Zero(1, cols), matrix::Zero(1, cols), matrix::Zero(1, cols)};

// Z_row_vectors_vector = {Z41, Z42, Z43};
col_vector Z_col_vectors_vector[3] = {matrix::Zero(rows, 1), matrix::Zero(rows, 1), matrix::Zero(rows, 1)};
col_vector Z_dot_col_vectors_vector[3] = {matrix::Zero(rows, 1), matrix::Zero(rows, 1), matrix::Zero(rows, 1)};

// Initialising all phis and phi_dots


// Initialising all Hs 
matrix H(matrix A, )

// Initialising all Gs


// Initialising all Fs


// Initialising all Ys

// Total energy
long double H( long double g, l_or_d* X_scalars_vector, matrix* X_matrices_vector, 
               l_or_d* V_scalars_vector, matrix* V_matrices_vector, l_or_d* A_scalars_vector,
               matrix* A_matrices_vector)
{
    l_or_d X_sum = 0.0;
    // Compute kinetic energy K
    for(l_or_d X: X_scalars_vector)
    {
        X_sum += X*X;
    }

    for(matrix X: X_matrices_vector)
    {
        X_sum += (X*X).trace();
    }

    l_or_d Z_sum = 0.0;
    for()
    {
        Z_sum 
    }    

    l_or_d K = 0.5 * (  X_sum + Z_sum);



    return K - 0.5*g*g*(V_gauge, V_D, V_F);
}

// Cillian's Gauss' law
matrix gauss_law(
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    matrix result =  matrix::Zero(rows, cols);
    for (int i = 0; i < 9; i ++)
    {
        result = commutator(X1,V1) + commutator(X2,V2) + commutator(X3,V3) + commutator(X4,V4) + commutator(X5,V5) +
        commutator(X6,V6) + commutator(X7,V7) + commutator(X8,V8) + commutator(X9,V9);
    }
    return result;
}


// Random Hermitian matrix

// Acceleration of each coordinate matrix
matrix Acceleration(const int j, matrix* X_vector, int rows, int cols, const long double g)
{

    matrix temp_commutator = matrix::Zero(rows, cols);
    matrix commutator_sum = matrix::Zero(rows, cols);

    matrix X = X_vector[j];
    for (int i = 0; i < dim; ++i)
    {
        if (i != j)
        {  
            temp_commutator = commutator(X_vector[i], commutator(X, X_vector[i]));
            

            commutator_sum += temp_commutator;
            
        }   
    // Comment out the line below to go back to the lagrangian e.q.m. With no -1/(g^2)
    // commutator_sum = -1.0 / (g * g) * commutator_sum;
    }

    return g*g*commutator_sum;
}

int main() 
{
    std::cout << "g" << g << '\n';
    // Create  vectors to store the matrices
    matrix X_vector[dim] = {matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols)};
    matrix V_vector[dim]= {matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols)};
    matrix A_vector[dim]= {matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols), matrix::Zero(rows, cols)};

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
                X_vector[i](row, col) = X_vector[i](row, col)/g;
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


    // Create  vectors to store the new matrices
    matrix X_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix V_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix A_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};

    // Defining all phis and phi_dots




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
            A_vector_new[i] = Acceleration( i, X_vector_new, rows, cols, g);
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
            std::cout << "H" << std::setprecision(15) << H(g, 
                            X_vector_new[0], X_vector_new[1], X_vector_new[2], X_vector_new[3], X_vector_new[4], X_vector_new[5], X_vector_new[6], X_vector_new[7], X_vector_new[8],
                          V_vector_new[0], V_vector_new[1], V_vector_new[2], V_vector_new[3], V_vector_new[4], V_vector_new[5], V_vector_new[6], V_vector_new[7], V_vector_new[8]);


        }

    }




       // Export initial X/V/A_vector to text files to be analysed in python.
    std:: fstream X2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D2-Branes/thermalised_X.txt", std:: ios:: out);
    X2_vector_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (matrix Matrix : X_vector_new)
    {
        X2_vector_Export << Matrix << std::endl;
    }


    std:: fstream V2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D2-Branes/thermalised_V.txt", std:: ios:: out);
    V2_vector_Export << std::fixed << std::setprecision(15);
    for (matrix Matrix : V_vector_new)
    {
        V2_vector_Export << Matrix << std::endl;
    }

    std:: fstream A2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D2-Branes/thermalised_A.txt", std:: ios:: out);
    A2_vector_Export << std::fixed << std::setprecision(15);
    // Print to text file
    for (matrix Matrix : A_vector_new)
    {
        A2_vector_Export << Matrix << std::endl;
    }


    X2_vector_Export.close();
    V2_vector_Export.close();
    A2_vector_Export.close();


    std::cout << "Finished" << std::time(nullptr)-start;
    return 0;
}
