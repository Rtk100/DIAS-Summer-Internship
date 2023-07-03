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
const double delta_t = 1e-5;
const double seconds_perturbed = 0.1;
const double g = 1;

// Repeat simulation for 1000 seconds.
const int simulation_repetitions = seconds_perturbed / delta_t;
// Number of D0-Branes
const int N = 2;
const int rows = N;
const int cols = N;

// Dimension of space
const int dim = 9;

//typedef std::complex<double> R_or_C;
//typedef Eigen:: Matrix<std::complex<double>, N, N> matrix;

typedef double R_or_C;
typedef Eigen:: Matrix<double, N, N> matrix;

matrix commutator(matrix A, matrix B)
{
    return A * B - B * A;
} 

matrix anti_commutator(matrix A, matrix B)
{
    return A * B + B * A;
} 

static std::random_device rd;
static std::mt19937 rng(std::time(nullptr)); 
std::normal_distribution<double> dist(0.0, 1e-8);

// Cillian's Hamiltonian
double H(
    double g, 
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    // Compute kinetic energy T
    R_or_C T = 1.0/(2.0*pow(g,2)) * (V1*V1 + V2*V2 + V3*V3 + V4*V4 + V5*V5 + V6*V6 + V7*V7 + V8*V8 + V9*V9).trace();

    matrix X[9] = {X1,X2,X3,X4,X5,X6,X7,X8,X9}; 

    matrix commutator_sum = matrix::Zero(rows, cols);  
    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            if(i == j)
                continue;
            commutator_sum += commutator(X[i],X[j])*commutator(X[i],X[j]); //can likely be more efficient by less function calls
        }
    }
    R_or_C U = - 1.0/(4.0*pow(g,2)) * commutator_sum.trace();
    return std:: abs(T + U);
}

// Cillian's Gauss' law
matrix gauss_law(matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
                 matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    matrix result;
    for (int i = 0; i < 9; i++)
    {
        result = commutator(X1,V1) + commutator(X2,V2) + commutator(X3,V3) + commutator(X4,V4) + commutator(X5,V5) +
        commutator(X6,V6) + commutator(X7,V7) + commutator(X8,V8) + commutator(X9,V9);
    }
    return result;
}

// Acceleration of each coordinate matrix
matrix PerturbingAcceleration(const int j, matrix* X_vector, int rows, int cols, const double g)
{
    matrix commutator_sum = matrix::Zero(rows, cols);
    matrix sum_X = matrix::Zero(rows, cols);

    R_or_C c_1= 0.0; 
    R_or_C c_2= 0.0;

    matrix X = X_vector[j];
    for (int i = 0; i < dim; i++)
    {
        c_1 = dist(rng);
        c_2 = dist(rng);
        sum_X = matrix::Zero(rows, cols);

        if (i != j)
        {  
            matrix temp_commutator = commutator(X_vector[i], commutator(X, X_vector[i]));
            for( int k = 0; k < dim; k++)
            {
                if (k != i)
                {
                    
                    sum_X += X_vector[k] * X_vector[k];
                }

            }
            matrix perturbation = 2.0 * c_1 * X + 2.0 * c_2 * anti_commutator(X, sum_X );

            commutator_sum += (temp_commutator + perturbation);


        }   

    // Comment out the line below to go back to the lagrangian e.q.m. With no -1/(g^2)
    //commutator_sum = -1.0 / (g * g) * commutator_sum;
    }
    return commutator_sum;
}


int main() 
{



    // Create  vectors to store the matrices
    matrix X2_vector[dim];
    matrix V2_vector[dim];
    matrix A2_vector[dim];



    // Generate and store X1, X2, X3, X4, X5, X6, X7, X8, and X9
        // Create an array to store the matrices
    // Open the text file for reading
    std::ifstream inputX("Thermalised_X.txt");
    if (!inputX.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }

    // Read the values from the file and store them in the matrices
    for (int i = 0; i < dim; ++i) 
    {

        for (int row = 0; row < rows; ++row) 
        {
            for (int col = 0; col < cols; ++col) 
            {
                 inputX >> X2_vector[i](row, col);
            }
        }
    }

    // Close the input file
    inputX.close();

        // Create an array to store the matrices
    // Open the text file for reading
    std::ifstream inputV("Thermalised_V.txt");
    if (!inputV.is_open()) {
        std::cerr << "Failed to open the V file." << std::endl;
        return 1;
    }


    // Read the values from the file and store them in the matrices
    for (int i = 0; i < dim; ++i) 
    {

        for (int row = 0; row < rows; ++row) 
        {
            for (int col = 0; col < cols; ++col) 
            {
                 inputV >> V2_vector[i](row, col);
            }
        }
    }

    // Close the input file
    inputV.close();

        // Create an array to store the matrices
    // Open the text file for reading
    std::ifstream inputA("Thermalised_A.txt");
    if (!inputA.is_open()) {
        std::cerr << "Failed to open the file." << std::endl;
        return 1;
    }

    // Read the values from the file and store them in the matrices
    for (int i = 0; i < dim; ++i) 
    {

        for (int row = 0; row < rows; ++row) 
        {
            for (int col = 0; col < cols; ++col) 
            {
                 inputA >> A2_vector[i](row, col);
            }
        }
    }

    // Close the input file
    inputA.close();

// If you want to test if the coords are loaded in right (They are)


    std::cout << "Start" << std::endl << "X";

    // Print out starting values
    for (matrix Matrix : X2_vector)
    {
        std::cout << Matrix << std::endl;
    }

    std::cout <<std::endl << "V";

    // Print out starting values
    for (matrix Matrix : V2_vector)
    {
        std::cout << Matrix << std::endl;
    }

    std::cout << std::endl << "A";

    // Print out starting values
    for (matrix Matrix : A2_vector)
    {
        std::cout << Matrix << std::endl;
    }



    matrix zero_matrix = matrix::Zero(rows, cols);

    matrix X2_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix V2_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    matrix A2_vector_new[9] = {zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix, zero_matrix};
    
    // Perturb the system to get the two set of conditions, for 1 second before turning off the deformation.
    for (int j = 0; j < 1.0 / delta_t; ++j)
    {
        //for (int h = 0; h < dim; h++)
        //{
        //    std::cout << std::endl << X2_vector[h];
       // }
        // velocity Verlet 1 to get new positions from old positions, momentums and rate of change of momentums

        for (int i = 0; i < 9; ++i)
        {
            X2_vector_new[i] = X2_vector[i] + V2_vector[i] * delta_t + 0.5 * A2_vector[i] * delta_t * delta_t;
        }

        // Generate and store new A1, A2, A3, A4, A5, A6, A7, A8, and A9
        for (int i = 0; i < 9; ++i) 
        {
            A2_vector_new[i] = PerturbingAcceleration( i, X2_vector_new, rows, cols, g);
        }  
        
        // Use Velocity Verlet 2 to get new momentums
        for (int i = 0; i < 9; ++i)
        {
            V2_vector_new[i] = V2_vector[i] + 0.5 * (A2_vector_new[i] + A2_vector[i]) * delta_t;
        }
        
        // Copy elements from X_vector_new to X_vector
        std::memcpy(X2_vector, X2_vector_new, sizeof(X2_vector_new));  


        // Copy elements from X_vector_new to X_vector
        std::memcpy(V2_vector, V2_vector_new, sizeof(V2_vector_new));  

        // Copy elements from X_vector_new to X_vector
        std::memcpy(A2_vector, A2_vector_new, sizeof(A2_vector_new));  

        if (j % 500== 0)
        {
            std::cout << '\n';
            std::cout << "H" << std::setprecision(15) << H(g, 
                           X2_vector_new[0], X2_vector_new[1], X2_vector_new[2], X2_vector_new[3], X2_vector_new[4], X2_vector_new[5], X2_vector_new[6], X2_vector_new[7], X2_vector_new[8],
                            V2_vector_new[0], V2_vector_new[1], V2_vector_new[2], V2_vector_new[3], V2_vector_new[4], V2_vector_new[5], V2_vector_new[6], V2_vector_new[7], V2_vector_new[8]);
           std::cout << '\n' << std::setprecision(15) << "Gauss" << gauss_law(X2_vector_new[0], X2_vector_new[1], X2_vector_new[2], X2_vector_new[3], X2_vector_new[4], X2_vector_new[5], X2_vector_new[6], X2_vector_new[7], X2_vector_new[8],
                                                               V2_vector_new[0], V2_vector_new[1], V2_vector_new[2], V2_vector_new[3], V2_vector_new[4], V2_vector_new[5], V2_vector_new[6], V2_vector_new[7], V2_vector_new[8]);
        }
        

    }



    // Export initial X/V/A_vector to text files to be analysed in python.
    std:: fstream X2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-Branes/perturbed_X.txt", std:: ios:: out);
    X2_vector_Export << std::fixed << std::setprecision(15);
    //Print to text file
    for (matrix Matrix : X2_vector)
    {
        X2_vector_Export << Matrix << std::endl;
    }

    std:: fstream V2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-Branes/perturbed_V.txt", std:: ios:: out);
    V2_vector_Export << std::fixed << std::setprecision(15);
    for (matrix Matrix : V2_vector)
    {
        V2_vector_Export << Matrix << std::endl;
    }

    std:: fstream A2_vector_Export("C:/Users/robtk/DIAS-Summer-Internship/C++/D0-Branes/perturbed_A.txt", std:: ios:: out);
    A2_vector_Export << std::fixed << std::setprecision(15);
    // Print to text file
    for (matrix Matrix : A2_vector)
    {
        A2_vector_Export << Matrix << std::endl;
    }

    X2_vector_Export.close();
    V2_vector_Export.close();
    A2_vector_Export.close();


    return 0;
}