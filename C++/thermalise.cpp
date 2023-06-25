#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <cmath>
#include <fstream>


// defining the Hamiltonian
double H(double lam1, double lam2, double x1, double x2, double y1, double y2,double px1, double px2, double py1, double py2)
{
    return 0.5*pow(lam1,2)*pow(x2,6) + 1.5*pow(lam1,2)*pow(x2,4)*pow(y2,2) + 1.5*pow(lam1,2)*pow(x2,2)*pow(y2,4) 
    + 0.5*pow(lam1,2)*pow(y2,6) + lam1*x1*pow(x2,3) - 3*lam1*x1*x2*pow(y2,2) + 3*lam1*pow(x2,2)*y1*y2 - lam1*y1*pow(y2,3)
    +0.5*pow(lam2,2)*pow(x1,6) + 1.5*pow(lam2,2)*pow(x1,4)*pow(y1,2) + 1.5*pow(lam2,2)*pow(x1,2)*pow(y1,4) + 0.5*pow(lam2,2)*pow(y1,6)
    +lam2*pow(x1,3)*x2 + 3*lam2*pow(x1,2)*y1*y2 - 3*lam2*x1*x2*pow(y1,2) -lam2*pow(y1,3)*y2 + 0.5*pow(px1,2) + 0.5*pow(px2,2) 
    +0.5*pow(py1,2) + 0.5*pow(py2,2) + 0.5*pow(x1,2) + 0.5*pow(x2,2) + 0.5*pow(y1,2) + 0.5*pow(y2,2); 
}

// defining derivatives of Hamiltonian 
double H_x1(double lam1, double lam2, double x1,double x2,double y1, double y2){

    return lam1*pow(x2,3) - 3*lam1*x2*pow(y2,2) + 3*pow(lam2,2)*pow(x1,5) + 6*pow(lam2,2)*pow(x1,3)*pow(y1,2)
    + 3*pow(lam2,2)*x1*pow(y1,4) + 3*lam2*pow(x1,2)*x2 +6*lam2*x1*y1*y2 -3*lam2*x2*pow(y1,2) + x1;
}

double H_x2(double lam1, double lam2, double x1,double x2,double y1, double y2){
    return 3*pow(lam1,2)*pow(x2,5) + 6*pow(lam1,2)*pow(x2,3)*pow(y2,2) + 3*pow(lam1,2)*x2*pow(y2,4) + 3*lam1*x1*pow(x2,2)
    - 3*lam1*x1*pow(y2,2) + 6*lam1*x2*y1*y2 + lam2*pow(x1,3) - 3*lam2*x1*pow(y1,2) + x2;
}

double H_y1(double lam1, double lam2, double x1,double x2,double y1, double y2){
    return 3*lam1*pow(x2,2)*y2 - lam1*pow(y2,3) + 3*pow(lam2,2)*pow(x1,4)*y1 + 6*pow(lam2,2)*pow(x1,2)*pow(y1,3) 
    + 3*pow(lam2,2)*pow(y1,5) + 3*lam2*pow(x1,2)*y2 - 6*lam2*x1*x2*y1 - 3*lam2*pow(y1,2)*y2 + y1;
}

double H_y2(double lam1, double lam2, double x1,double x2,double y1, double y2){
    return 3*pow(lam1,2)*pow(x2,4)*y2 + 6*pow(lam1,2)*pow(x2,2)*pow(y2,3) + 3*pow(lam1,2)*pow(y2,5) - 6*lam1*x1*x2*y2
    +3*lam1*pow(x2,2)*y1 - 3*lam1*y1*pow(y2,2) + 3*lam2*pow(x1,2)*y1 - lam2*pow(y1,3) + y2;
}

void update(
    double lam1,double lam2,double dt,double* x1,double* x2,double* y1,double* y2,double* px1,double* px2,double* py1,double* py2,
    double* H_x1_0, double* H_x2_0, double* H_y1_0, double* H_y2_0, double* x1_n, double* x2_n, double* y1_n, double* y2_n,
    double* H_x1_n, double* H_x2_n, double* H_y1_n, double* H_y2_n)
{
    

    *x1_n = *x1 + (*px1)*dt - 0.5*(*H_x1_0)*pow(dt,2);
    *x2_n = *x2 + (*px2)*dt - 0.5*(*H_x2_0)*pow(dt,2);
    *y1_n = *y1 + (*py1)*dt - 0.5*(*H_y1_0)*pow(dt,2);
    *y2_n = *y2 + (*py2)*dt - 0.5*(*H_y2_0)*pow(dt,2);

    *H_x1_n = H_x1(lam1,lam2,*x1_n,*x2_n,*y1_n,*y2_n);
    *H_x2_n = H_x2(lam1,lam2,*x1_n,*x2_n,*y1_n,*y2_n);
    *H_y1_n = H_y1(lam1,lam2,*x1_n,*x2_n,*y1_n,*y2_n);
    *H_y2_n = H_y2(lam1,lam2,*x1_n,*x2_n,*y1_n,*y2_n);

    *px1 = *px1 - 0.5*(*H_x1_0 + *H_x1_n)*dt;
    *px2 = *px2 - 0.5*(*H_x2_0 + *H_x2_n)*dt;
    *py1 = *py1 - 0.5*(*H_y1_0 + *H_y1_n)*dt;
    *py2 = *py2 - 0.5*(*H_y2_0 + *H_y2_n)*dt;

    *x1 = *x1_n;
    *x2 = *x2_n;
    *y1 = *y1_n;
    *y2 = *y2_n;

    *H_x1_0 = *H_x1_n;
    *H_x2_0 = *H_x2_n;
    *H_y1_0 = *H_y1_n;
    *H_y2_0 = *H_y2_n;
}

const int iterations = 1e7;

int main()
{
    // initializing lambdas, time step and iterations
    const double lam1 = 0.0001, lam2 = 0.0002; 
    const double dt = 1e-4;
    
    // provide initial conditions for coordinates
    double x1 = 100, x2 = -200, y1 = 170, y2 = 200, px1 = 90, px2 = 160, py1 = 200, py2 = -300;

    double H_x1_0 = H_x1(lam1,lam2,x1,x2,y1,y2);
    double H_x2_0 = H_x2(lam1,lam2,x1,x2,y1,y2);
    double H_y1_0 = H_y1(lam1,lam2,x1,x2,y1,y2);
    double H_y2_0 = H_y2(lam1,lam2,x1,x2,y1,y2);

    double x1_n = x1 + (px1)*dt - 0.5*(H_x1_0)*pow(dt,2);
    double x2_n = x2 + (px2)*dt - 0.5*(H_x2_0)*pow(dt,2);
    double y1_n = y1 + (py1)*dt - 0.5*(H_y1_0)*pow(dt,2);
    double y2_n = y2 + (py2)*dt - 0.5*(H_y2_0)*pow(dt,2);

    double H_x1_n = H_x1(lam1,lam2,x1_n,x2_n,y1_n,y2_n);
    double H_x2_n = H_x2(lam1,lam2,x1_n,x2_n,y1_n,y2_n);
    double H_y1_n = H_y1(lam1,lam2,x1_n,x2_n,y1_n,y2_n);
    double H_y2_n = H_y2(lam1,lam2,x1_n,x2_n,y1_n,y2_n);

    for (int i = 1; i < iterations; i++)
    {
        update(lam1,lam2,dt,&x1,&x2,&y1,&y2,&px1,&px2,&py1,&py2,&H_x1_0,&H_x2_0,&H_y1_0,&H_y2_0,&x1_n,&x2_n,&y1_n,&y2_n,&H_x1_n,&H_x2_n,&H_y1_n,&H_y2_n);
    }

    // Create and open a text file
    std:: ofstream result("C:/Users/cilli/DIAS.cpp/Run_files4.0/thermalised_coordinates.txt"/* specify where to save the thermalised configuration*/);

    // Write to the file
    result << x1 << "\n" << x2 << "\n" << y1 << "\n" << y2 << "\n" << 
    px1 << "\n" << px2 << "\n" << py1 << "\n" << py2 << std:: endl;

    // Close the file
    result.close();

    /* 
    This can be used to check if the hamiltonian is being conserved
    
    for (int i = 0; i < iterations / record; i++)
    {
        
        std:: cout << std:: setprecision(16);
        std:: cout << H(lam1,lam2,x1_sol[i],x2_sol[i],y1_sol[i],y2_sol[i],px1_sol[i],px2_sol[i],py1_sol[i],py2_sol[i]) << std:: endl;    
    }
    */

    // Write the current hamiltonian into a file to be used in perturbation of state
    std:: fstream Rhoensfile("C:/Users/cilli/DIAS.cpp/Run_files4.0/CurrentHamiltonian1.0.txt", std:: ios:: out);
    Rhoensfile << H(lam1,lam2,x1,x2,y1,y2,px1,px2,py1,py2);
    Rhoensfile.close();

    return 0;
}