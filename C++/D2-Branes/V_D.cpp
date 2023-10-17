#include <iostream>
#include <random>
#include <ctime>
#include <string>
#include <complex>
#include <fstream>
#include <sstream>
#include "eigen/Eigen/Dense"

/*
To run the code you will need to install the Eigen library and included as above. 
*/

//Initializing some constants, g, N (length of vectors), time step dt, and number of iterations
const long double g = 1;
const int N = 3;
const long double dt = 1e-4;
const int iterations = 10000;

/* 
Here I redefine the names for matrix, row and column data types because they will be used a lot 
Defining complex numbers as complex;
Defining complex matrices as matrix ;
Defining complex rows as row;
Defining complex column as col;
*/

typedef std::complex<long double> complex;
typedef std::complex<long double> Complex;

typedef Eigen:: Matrix<complex, N, N> matrix;
typedef Eigen:: Matrix<complex, 1, N> row;
typedef Eigen:: Matrix<complex, N, 1> col;

/*
Defining functions that return the modulus squared of the columnsl, rows, and 
just for a complex number because they will be used a lot
*/

long double colabsqr(col a)
{
    complex b = (a.adjoint()*a);
    return b.real();
}

long double rowabsqr(row a)
{
    complex b = (a*a.adjoint());
    return b.real();
}

long double modsqr(complex a)
{
    return pow(abs(a),2);
}

// Initialise a random engine to be used when creating random initial Z coordinates 
std:: mt19937 rng(10);
std:: normal_distribution<long double> gauss_dist(0, 1);

// Initialise the Z variables only and let U denote the \dot{Z} coordinates
/*
complex Z12 = complex(gauss_dist(rng),gauss_dist(rng));
complex Z13 = complex(gauss_dist(rng),gauss_dist(rng));
complex Z21 = complex(gauss_dist(rng),gauss_dist(rng));
complex Z23 = complex(gauss_dist(rng),gauss_dist(rng)); 
complex Z31 = complex(gauss_dist(rng),gauss_dist(rng)); 
complex Z32 = complex(gauss_dist(rng),gauss_dist(rng)); 

*/
complex Z12 = complex(0.207606,0.366979);
complex Z13 = complex(0.0365052,-0.344526);
complex Z21 = complex(0.0786534,-1.71609);
complex Z23 = complex(0.508537,1.50702); 
complex Z31 = complex(0.547891,0.496717); 
complex Z32 = complex(-0.937337,-1.09565); 



// initialise the Z_n which will be needed in the update function
complex Z12_n = complex(0,0);
complex Z13_n = complex(0,0);
complex Z21_n = complex(0,0);
complex Z23_n = complex(0,0); 
complex Z31_n = complex(0,0); 
complex Z32_n = complex(0,0); 

col Z41_n = col::Zero();
col Z42_n = col::Zero();
col Z43_n = col::Zero();
row Z14_n = row::Zero();
row Z24_n = row::Zero(); 
row Z34_n = row::Zero();

// Initialise U
complex U12 = complex(0,0);
complex U13 = complex(0,0);
complex U21 = complex(0,0);
complex U23 = complex(0,0);
complex U31 = complex(0,0);
complex U32 = complex(0,0);

col U41 = col::Zero();
col U42 = col::Zero();
col U43 = col::Zero();
row U14 = row::Zero();
row U24 = row::Zero();
row U34 = row::Zero();

// Initialise c's as random complex numbers
complex c1 = 2.5;
complex c2 = 3.3;
complex c3 = 1.2;
complex c4 = -c1 -c2 -c3;

// Define functions for the potential V_D and the Kinetic energy K
long double K(complex U12,complex U13,row U14,complex U21,complex U23,row U24,complex U31,complex U32,row U34, col U41, col U42,col U43)
{
    // For k = 1
    long double sum1 = (conj(U12)*U12 + conj(U13)*U13 + (U14.adjoint()*U14).trace()).real();
    
    // For k = 2
    long double sum2 = (conj(U21)*U21 + conj(U23)*U23 + (U24.adjoint()*U24).trace()).real();

    // For k = 3
    long double sum3 = (conj(U31)*U31 + conj(U32)*U32 + (U34.adjoint()*U34).trace()).real();

    // For k = 4
    long double sum4 = (U41.adjoint()*U41 + U42.adjoint()*U42 + U43.adjoint()*U43).trace().real();

    return 0.5 * (sum1 + sum2 + sum3 + sum4);
}

long double V_D(complex Z12,complex Z13,row Z14,complex Z21,complex Z23,row Z24,complex Z31,complex Z32,row Z34, col Z41, col Z42,col Z43)
{
    // For k = 1
    long double sum1 = pow((Z12*conj(Z12) + Z13*conj(Z13) + (Z14*Z14.adjoint()).trace() - conj(Z21)*Z21 - conj(Z31)*Z31 - (Z41.adjoint()*Z41).trace() - c1/(g*g)).real(),2);

    // For k = 2
    long double sum2 = pow((Z21*conj(Z21) + Z23*conj(Z23) + (Z24*Z24.adjoint()).trace() - conj(Z12)*Z12 - conj(Z32)*Z32 - (Z42.adjoint()*Z42).trace() - c2/(g*g)).real(),2);

    // For k = 3
    long double sum3 = pow((Z31*conj(Z31) + Z32*conj(Z32) + (Z34*Z34.adjoint()).trace() - conj(Z13)*Z13 - conj(Z23)*Z23 - (Z43.adjoint()*Z43).trace() - c3/(g*g)).real(),2);

    // For k = 4
    matrix sum4 = (
        Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34 - c4/(g*g)*matrix::Identity()) *
        (Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34 - c4/(g*g)*matrix::Identity()
        );

    return 0.5 * (sum1 + sum2 + sum3 + sum4.trace().real());
}

// Define the Force functions (defined as in dropbox -\ddot(Zij) = FZij)
complex FZ12(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z24, complex Z32, col Z42)
{
    return Z12*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    - Z12*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g));
}

complex FZ21(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z24, complex Z32, col Z42)
{
    return - Z21*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    + Z21*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g));
}

complex FZ13(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z34, complex Z32, col Z43)
{
    return Z13*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    - Z13*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g)); 
}

complex FZ31(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z34, complex Z32, col Z43)
{
    return - Z31*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    + Z31*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g));
}

complex FZ23(
    complex Z23, complex Z21, row Z24, complex Z12, complex Z32, col Z42,
    complex Z31, row Z34, complex Z13, col Z43)
{
    return Z23*(modsqr(Z23) + modsqr(Z21) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g))
    - Z23*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g));
}

complex FZ32(
    complex Z23, complex Z21, row Z24, complex Z12, complex Z32, col Z42,
    complex Z31, row Z34, complex Z13, col Z43)
{
    return - Z32*(modsqr(Z23) + modsqr(Z21) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g))
    + Z32*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g));
}

row FZ14(
    row Z14, complex Z12, complex Z13, complex Z21, complex Z31, col Z41, 
    col Z42, col Z43, row Z24, row Z34)
{
    return Z14*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) -c1/(g*g))
    - Z14*(Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    - c4/(g*g) * matrix::Identity()); 
}

col FZ41(
    row Z14, complex Z12, complex Z13, complex Z21, complex Z31, col Z41, 
    col Z42, col Z43, row Z24, row Z34)
{
    return - Z41*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) -c1/(g*g))
    + (Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    - c4/(g*g) * matrix::Identity())*Z41;
}

row FZ24(
    row Z24, complex Z21, complex Z23, complex Z12, complex Z32, col Z41,
    col Z42, col Z43, row Z14, row Z34)
{
    return Z24*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g))
    - Z24*(Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    - c4/(g*g) * matrix::Identity());
}

col FZ42(
    row Z24, complex Z21, complex Z23, complex Z12, complex Z32, col Z41,
    col Z42, col Z43, row Z14, row Z34)
{
    return - Z42*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g))
    + (Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    - c4/(g*g) * matrix::Identity())*Z42;
}

row FZ34(
    row Z34, complex Z31, complex Z32, complex Z13, complex Z23, col Z41,
    col Z42, col Z43, row Z14, row Z24)
{
    return Z34*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g))
    - Z34*(Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    - c4/(g*g) * matrix::Identity());
}

col FZ43(
    row Z34, complex Z31, complex Z32, complex Z13, complex Z23, col Z41,
    col Z42, col Z43, row Z14, row Z24)
{
    return - Z43*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(Z13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g))
    + (Z41*Z41.adjoint() + Z42*Z42.adjoint() + Z43*Z43.adjoint() - Z14.adjoint()*Z14 - Z24.adjoint()*Z24 - Z34.adjoint()*Z34
    - c4/(g*g) * matrix::Identity())*Z43;
}

// Define the Update function... must make sure that all of the paramaters of the force functions are in the 
// Correct order or else they will give unseen errors

/*
Update function description: The update function takes in all of the generalised coordinates and velocites and uses the Verlet algorithm
to "update" the coordinates to their value ate the next time step (dt later). Zij_n is necessary because to update the velocities we 
must take an average of the Force function now and at a time dt later. Thus we need the generalised coordinates at a time dt later (the Zij_n)
but we must still hold onto our generalised coordinates at the current time (Zij). Once the velocity has been updated we set 
Zij to be replaced by Zij_n and then repeat. Note I am using pointers in this function so the function does not return some value but 
instead it just goes to the memory address of the current generalised coordinates and replaces them with their values at a time dt later. 
If things are done this way then when the function is called you must pass the address of the variables as arguments rather than just the 
variables themselves. To access the address of a variable you use the ampersand symbol &. So when the update function is called in the 
code the arguments will be update(dt, &Z12, &Z12_n, &Z21, &Z21_n, ......., &U43) 
*/

void update(
    long double dt,
    complex* Z12, complex* Z12_n, complex* Z21, complex* Z21_n, complex* Z13, complex* Z13_n, complex* Z31, complex* Z31_n,
    complex* Z23, complex* Z23_n, complex* Z32, complex* Z32_n, col* Z41, col* Z41_n, col* Z42, col* Z42_n, col* Z43,
    col* Z43_n, row* Z14, row* Z14_n, row* Z24, row* Z24_n, row* Z34, row* Z34_n,
    complex* U12, complex* U13, complex* U21, complex* U23, complex* U31, complex* U32,
    row* U14, row* U24, row* U34, col* U41, col* U42, col* U43)
{
    *Z12_n = *Z12 + *U12*dt - complex(0.5, 0) * FZ12(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z24, *Z32, *Z42) * dt * dt;
    *Z13_n = *Z13 + *U13*dt - complex(0.5, 0) * FZ13(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z34, *Z32, *Z43) * dt * dt;
    *Z21_n = *Z21 + *U21*dt - complex(0.5, 0) * FZ21(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z24, *Z32, *Z42) * dt * dt;
    *Z23_n = *Z23 + *U23*dt - complex(0.5, 0) * FZ23(*Z23, *Z21, *Z24, *Z12, *Z32, *Z42, *Z31, *Z24, *Z13, *Z43) * dt * dt;
    *Z31_n = *Z31 + *U31*dt - complex(0.5, 0) * FZ31(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z34, *Z32, *Z43) * dt * dt;
    *Z32_n = *Z32 + *U32*dt - complex(0.5, 0) * FZ32(*Z23, *Z21, *Z24, *Z12, *Z32, *Z42, *Z31, *Z34, *Z31, *Z43) * dt * dt;
    *Z14_n = *Z14 + *U14*dt - complex(0.5, 0) * FZ14(*Z14, *Z12, *Z13, *Z21, *Z31, *Z41, *Z42, *Z43, *Z24, *Z34) * dt * dt;
    *Z24_n = *Z24 + *U24*dt - complex(0.5, 0) * FZ24(*Z24, *Z21, *Z23, *Z12, *Z32, *Z41, *Z42, *Z43, *Z14, *Z34) * dt * dt;
    *Z34_n = *Z34 + *U34*dt - complex(0.5, 0) * FZ34(*Z34, *Z31, *Z32, *Z13, *Z23, *Z41, *Z42, *Z43, *Z14, *Z24) * dt * dt;
    *Z41_n = *Z41 + *U41*dt - complex(0.5, 0) * FZ41(*Z14, *Z12, *Z13, *Z21, *Z31, *Z41, *Z42, *Z43, *Z24, *Z34) * dt * dt;
    *Z42_n = *Z42 + *U42*dt - complex(0.5, 0) * FZ42(*Z24, *Z21, *Z23, *Z12, *Z32, *Z41, *Z42, *Z43, *Z14, *Z34) * dt * dt;
    *Z43_n = *Z43 + *U43*dt - complex(0.5, 0) * FZ43(*Z34, *Z31, *Z32, *Z13, *Z23, *Z41, *Z42, *Z43, *Z14, *Z24) * dt * dt;

    *U12 = *U12 - complex(0.5, 0) * (FZ12(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z24, *Z32, *Z42) + FZ12(*Z12_n, *Z13_n, *Z14_n, *Z21_n, *Z31_n, *Z41_n, *Z23_n, *Z24_n, *Z32_n, *Z42_n)) * dt;
    *U13 = *U13 - complex(0.5, 0) * (FZ13(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z34, *Z32, *Z43) + FZ13(*Z12_n, *Z13_n, *Z14_n, *Z21_n, *Z31_n, *Z41_n, *Z23_n, *Z34_n, *Z32_n, *Z43_n)) * dt;
    *U21 = *U21 - complex(0.5, 0) * (FZ21(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z24, *Z32, *Z42) + FZ21(*Z12_n, *Z13_n, *Z14_n, *Z21_n, *Z31_n, *Z41_n, *Z23_n, *Z24_n, *Z32_n, *Z42_n)) * dt;
    *U23 = *U23 - complex(0.5, 0) * (FZ23(*Z23, *Z21, *Z24, *Z12, *Z32, *Z42, *Z31, *Z24, *Z13, *Z43) + FZ23(*Z23_n, *Z21_n, *Z24_n, *Z12_n, *Z32_n, *Z42_n, *Z31_n, *Z24_n, *Z13_n, *Z43_n)) * dt;
    *U31 = *U31 - complex(0.5, 0) * (FZ31(*Z12, *Z13, *Z14, *Z21, *Z31, *Z41, *Z23, *Z34, *Z32, *Z43) + FZ31(*Z23_n, *Z21_n, *Z24_n, *Z12_n, *Z32_n, *Z42_n, *Z31_n, *Z24_n, *Z13_n, *Z43_n)) * dt;
    *U32 = *U32 - complex(0.5, 0) * (FZ32(*Z23, *Z21, *Z24, *Z12, *Z32, *Z42, *Z31, *Z34, *Z31, *Z43) + FZ32(*Z23_n, *Z21_n, *Z24_n, *Z12_n, *Z32_n, *Z42_n, *Z31_n, *Z34_n, *Z31_n, *Z43_n)) * dt;
    *U14 = *U14 - complex(0.5, 0) * (FZ14(*Z14, *Z12, *Z13, *Z21, *Z31, *Z41, *Z42, *Z43, *Z24, *Z34) + FZ14(*Z14_n, *Z12_n, *Z13_n, *Z21_n, *Z31_n, *Z41_n, *Z42_n, *Z43_n, *Z24_n, *Z34_n)) * dt;
    *U24 = *U24 - complex(0.5, 0) * (FZ24(*Z24, *Z21, *Z23, *Z12, *Z32, *Z41, *Z42, *Z43, *Z14, *Z34) + FZ24(*Z24_n, *Z21_n, *Z23_n, *Z12_n, *Z32_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z34_n)) * dt;
    *U34 = *U34 - complex(0.5, 0) * (FZ34(*Z34, *Z31, *Z32, *Z13, *Z23, *Z41, *Z42, *Z43, *Z14, *Z24) + FZ34(*Z34_n, *Z31_n, *Z32_n, *Z13_n, *Z23_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z24_n)) * dt;
    *U41 = *U41 - complex(0.5, 0) * (FZ41(*Z14, *Z12, *Z13, *Z21, *Z31, *Z41, *Z42, *Z43, *Z24, *Z34) + FZ41(*Z14_n, *Z12_n, *Z13_n, *Z21_n, *Z31_n, *Z41_n, *Z42_n, *Z43_n, *Z24_n, *Z34_n)) * dt;
    *U42 = *U42 - complex(0.5, 0) * (FZ42(*Z24, *Z21, *Z23, *Z12, *Z32, *Z41, *Z42, *Z43, *Z14, *Z34) + FZ42(*Z24_n, *Z21_n, *Z23_n, *Z12_n, *Z32_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z34_n)) * dt;
    *U43 = *U43 - complex(0.5, 0) * (FZ43(*Z34, *Z31, *Z32, *Z13, *Z23, *Z41, *Z42, *Z43, *Z14, *Z24) + FZ43(*Z34_n, *Z31_n, *Z32_n, *Z13_n, *Z23_n, *Z41_n, *Z42_n, *Z43_n, *Z14_n, *Z24_n)) * dt;

    *Z12 = *Z12_n;
    *Z13 = *Z13_n;
    *Z21 = *Z21_n;
    *Z23 = *Z23_n;
    *Z31 = *Z31_n;
    *Z32 = *Z32_n;
    *Z14 = *Z14_n;
    *Z24 = *Z24_n;
    *Z34 = *Z34_n;
    *Z41 = *Z41_n;
    *Z42 = *Z42_n;
    *Z43 = *Z43_n;
}

// Initialise solution arrays for storing the recorded coordinates
complex Z12_sol[iterations/1000];
complex Z13_sol[iterations/1000];
complex Z21_sol[iterations/1000];
complex Z23_sol[iterations/1000];
complex Z31_sol[iterations/1000];
complex Z32_sol[iterations/1000];
row Z14_sol[iterations/1000];
row Z24_sol[iterations/1000];
row Z34_sol[iterations/1000];
col Z41_sol[iterations/1000];
col Z42_sol[iterations/1000];
col Z43_sol[iterations/1000];

int main()
{
    row Z14;
    Z14(0,0) = Complex(1,1);
    Z14(0,1) = Complex(2,1);
    Z14(0,2) = Complex(3,1);
    row Z24;
    Z24(0,0) = Complex(1,2);
    Z24(0,1) = Complex(2,2);
    Z24(0,2) = Complex(3,2);
    row Z34;
    Z34(0,0) = Complex(1,3);
    Z34(0,1) = Complex(2,3);
    Z34(0,2) = Complex(3,3);
    col Z41;
    Z41(0,0) = Complex(1,1);
    Z41(1,0) = Complex(1,2);
    Z41(2,0) = Complex(1,3);
    col Z42;
    Z42(0,0) = Complex(2,1);
    Z42(1,0) = Complex(2,2);
    Z42(2,0) = Complex(2,3);
    col Z43;
    Z43(0,0) = Complex(3,1);
    Z43(1,0) = Complex(3,2);
    Z43(2,0) = Complex(3,3);
    std::cout << Z12 << '\n';
    std::cout << Z13 << '\n';
    std::cout << Z21 << '\n';
    std::cout << Z23 << '\n';
    std::cout << Z31 << '\n';
    std::cout << Z32 << '\n';
    std::cout << Z14 << '\n';
    std::cout << Z24 << '\n';
    std::cout << Z34 << '\n';
    std::cout << Z41 << '\n';
    std::cout << Z42 << '\n';
    std::cout << Z43 << '\n';
    for (int i = 0; i < iterations; i++)
    {
        // Record the values of the coordinates if i%1000 == 0 (i.e. record once every thousand iterations because we don't need every single value)
        if (i % 1000 == 0)
        {
            Z12_sol[i/1000] = Z12;
            Z13_sol[i/1000] = Z13;
            Z21_sol[i/1000] = Z21;
            Z23_sol[i/1000] = Z23;
            Z31_sol[i/1000] = Z31;
            Z32_sol[i/1000] = Z32;
            Z14_sol[i/1000] = Z14;
            Z24_sol[i/1000] = Z24;
            Z34_sol[i/1000] = Z34;
            Z41_sol[i/1000] = Z41;
            Z42_sol[i/1000] = Z42;
            Z43_sol[i/1000] = Z43;
            
            // Print the energy to see if its conserved and also print the progress (% of program that has run) 
            std::cout << "E = " << K(U12,U13,U14,U21,U23,U24,U31,U32,U34,U41,U42,U43) + V_D(Z12,Z13,Z14,Z21,Z23,Z24,Z31,Z32,Z34,Z41,Z42,Z43) << std::endl;
            std::cout << "Progress: " << (float(i)/iterations)*100 << "%" << std::endl;
        }
        
        update(
            dt, &Z12, &Z12_n, &Z21, &Z21_n, &Z13, &Z13_n, &Z31, &Z31_n, &Z23, &Z23_n, 
            &Z32, &Z32_n, &Z41, &Z41_n, &Z42, &Z42_n, &Z43, &Z43_n, &Z14, &Z14_n, 
            &Z24, &Z24_n, &Z34, &Z34_n, &U12, &U13, &U21, &U23, &U31, &U32, &U14, &U24, &U34, &U41, &U42, &U43);


    }
    return 0;
}