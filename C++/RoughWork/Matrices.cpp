#include <iostream>
#include "eigen/Eigen/Dense"

int main() 
{
  // Create a 3x3 matrix
  Eigen::Matrix3d matrix;

  // Set values in the matrix
  matrix << 1, 2, 3,
            4, 5, 6,
            7, 8, 90;

  // Print the matrix
  std::cout << "Matrix:\n" << matrix << std::endl;

  return 0;
}