#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <typeinfo>

// 
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

//Function to find the new px1 when the other coordinates have been perturbed
double get_px1_ds(double lam1, double lam2, double H, double x1, double x2, double y1, double y2, double px2, double py1, double py2)
{
    return sqrt(2*(H - 0.5*pow(lam1,2)*pow(x2,6) - 1.5*pow(lam1,2)*pow(x2,4)*pow(y2,2) - 1.5*pow(lam1,2)*pow(x2,2)*pow(y2,4)
    - 0.5*pow(lam1,2)*pow(y2,6) - lam1*x1*pow(x2,3) + 3*lam1*x1*x2*pow(y2,2) -3*lam1*pow(x2,2)*y1*y2 + lam1*y1*pow(y2,3)
    - 0.5*pow(lam2,2)*pow(x1,6) - 1.5*pow(lam2,2)*pow(x1,4)*pow(y1,2) - 1.5*pow(lam2,2)*pow(x1,2)*pow(y1,4) - 0.5*pow(lam2,2)*pow(y1,6)
    - lam2*pow(x1,3)*x2 - 3*lam2*pow(x1,2)*y1*y2 + 3*lam2*x1*x2*pow(y1,2) + lam2*pow(y1,3)*y2 - 0.5*pow(px2,2)
    - 0.5*pow(py1,2) - 0.5*pow(py2,2) - 0.5*pow(x1,2) - 0.5*pow(x2,2) - 0.5*pow(y1,2) - 0.5*pow(y2,2)));
}

//Function to update the coordinates at each time step
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

// Initialize iteration number and specify how often to record coordinates (save every 1000)
const int iterations = 200e4;
const int record = 1000;

// initialize solution arrays
double x1_sol[iterations / record], x2_sol[iterations / record], y1_sol[iterations / record], y2_sol[iterations / record],
px1_sol[iterations / record], px2_sol[iterations / record], py1_sol[iterations / record], py2_sol[iterations / record] = {};

double x1_sol_ds[iterations / record], x2_sol_ds[iterations / record], y1_sol_ds[iterations / record], y2_sol_ds[iterations / record],
px1_sol_ds[iterations / record], px2_sol_ds[iterations / record], py1_sol_ds[iterations / record], py2_sol_ds[iterations / record]; 

int main()
{    
    // initializing lambdas and time step
    const double lam1 = 0.0001, lam2 = 0.0002;
    const double dt = 1e-4;

    /*
    Sample number is the number of nearby points to simulate. Swapno said in the email to use 20-30 nearby points here I'm using 
    10 but can take as many as you like.
    */ 
    const int sample_number = 10;

    // Declaring variables for the initial conditions
    double x1, x2, y1, y2, px1, px2, py1, py2;
    double x1_ds, x2_ds, y1_ds, y2_ds, px1_ds, px2_ds, py1_ds, py2_ds;

    // Loading my initial conditions from the thermalised_coordinates file and assigning them to the variables declared above.
    std:: fstream initial_conditions("C:/Users/cilli/DIAS.cpp/Run_files4.0/thermalised_coordinates.txt", std:: ios:: in);
    initial_conditions >> x1 >> x2 >> y1 >> y2 >> px1 >> px2 >> py1 >> py2;
    initial_conditions.close();

    /*
    Loading in the current Hamiltonian for use in the get_px1_ds function defined outside of main. (the hamiltonian is
    computed and saved into the file CurrentHamiltonian.txt during the thermalise program)
    */ 
    double Hamiltonian;
    std:: fstream Rhoensfile("C:/Users/cilli/DIAS.cpp/Run_files4.0/CurrentHamiltonian1.0.txt", std:: ios :: in);
    Rhoensfile >> Hamiltonian;
    Rhoensfile.close();

    for (int i =0; i < sample_number; i++)
    {
        // Obtain perturbed coordinates 
        x1_ds = x1+0.1, x2_ds = x2+0.1, y1_ds = y1+0.1, y2_ds = y2+0.1, px2_ds = px2+0.1, py1_ds = py1+0.1, py2_ds = py2+0.1;
        px1_ds = get_px1_ds(lam1,lam2,Hamiltonian,x1_ds,x2_ds,y1_ds,y2_ds,px2_ds,py1_ds,py2_ds);

        // Setup variables for the update function 
        double H_x1_0 = H_x1(lam1,lam2,x1,x2,y1,y2), H_x1_0_ds = H_x1(lam1,lam2,x1_ds,x2_ds,y1_ds,y2_ds), 
        H_x2_0 = H_x2(lam1,lam2,x1,x2,y1,y2), H_x2_0_ds = H_x2(lam1,lam2,x1_ds,x2_ds,y1_ds,y2_ds), 
        H_y1_0 = H_y1(lam1,lam2,x1,x2,y1,y2), H_y1_0_ds = H_y1(lam1,lam2,x1_ds,x2_ds,y1_ds,y2_ds), 
        H_y2_0 = H_y2(lam1,lam2,x1,x2,y1,y2), H_y2_0_ds = H_y2(lam1,lam2,x1_ds,x2_ds,y1_ds,y2_ds),

        x1_n = x1 + (px1)*dt - 0.5*(H_x1_0)*pow(dt,2), x1_n_ds = x1_ds + (px1_ds)*dt - 0.5*(H_x1_0_ds)*pow(dt,2), 
        x2_n = x2 + (px2)*dt - 0.5*(H_x2_0)*pow(dt,2), x2_n_ds = x2_ds + (px2_ds)*dt - 0.5*(H_x2_0_ds)*pow(dt,2), 
        y1_n = y1 + (py1)*dt - 0.5*(H_y1_0)*pow(dt,2), y1_n_ds = y1_ds + (py1_ds)*dt - 0.5*(H_y1_0_ds)*pow(dt,2), 
        y2_n = y2 + (py2)*dt - 0.5*(H_y2_0)*pow(dt,2), y2_n_ds = y2_ds + (py2_ds)*dt - 0.5*(H_y2_0_ds)*pow(dt,2), 

        H_x1_n = H_x1(lam1,lam2,x1_n,x2_n,y1_n,y2_n), H_x1_n_ds = H_x1(lam1,lam2,x1_n_ds,x2_n_ds,y1_n_ds,y2_n_ds),
        H_x2_n = H_x2(lam1,lam2,x1_n,x2_n,y1_n,y2_n), H_x2_n_ds = H_x2(lam1,lam2,x1_n_ds,x2_n_ds,y1_n_ds,y2_n_ds),
        H_y1_n = H_y1(lam1,lam2,x1_n,x2_n,y1_n,y2_n), H_y1_n_ds = H_y1(lam1,lam2,x1_n_ds,x2_n_ds,y1_n_ds,y2_n_ds),
        H_y2_n = H_y2(lam1,lam2,x1_n,x2_n,y1_n,y2_n), H_y2_n_ds = H_y2(lam1,lam2,x1_n_ds,x2_n_ds,y1_n_ds,y2_n_ds);

        /*
        Checking that px1 and px1_ds are sufficiently close together (here sufficiently close being their distance less than 1)
        and also making sure we don't have Nan for px1_ds which you get if the sqrt in the get_px1_ds function is negative
        If both conditions not satisfied then update (i.e. simulate) for 1 second and check conditions again
        (The || in the condition statement is logical OR in c++)
        */
        while (sqrt(pow(px1 - px1_ds,2)) > 1 || std:: isnan(px1_ds))
        {
            for (int a = 0; a < 1000; a++)
            {
                update(lam1,lam2,dt,&x1,&x2,&y1,&y2,&px1,&px2,&py1,&py2,&H_x1_0,&H_x2_0,&H_y1_0,&H_y2_0,&x1_n,&x2_n,&y1_n,&y2_n,&H_x1_n,&H_x2_n,&H_y1_n,&H_y2_n);
                x1_ds = x1+0.1, x2_ds = x2+0.1, y1_ds = y1+0.1, y2_ds = y2+0.1, px2_ds = px2+0.1, py1_ds = py1+0.1, py2_ds = py2+0.1;
                px1_ds = get_px1_ds(lam1,lam2,Hamiltonian,x1_ds,x2_ds,y1_ds,y2_ds,px2_ds,py1_ds,py2_ds);    
            }
        }

        //Initialize first elements of solution arrays (Necessary to make update rule more efficient)
        x1_sol[0] = x1, x2_sol[0] = x2, y1_sol[0] = y1, y2_sol[0] = y2,
        px1_sol[0] = px1, px2_sol[0] = px2, py1_sol[0] = py1, py2_sol[0] = py2;

        x1_sol_ds[0] = x1_ds, x2_sol_ds[0] = x2_ds, y1_sol_ds[0] = y1_ds, y2_sol_ds[0] = y2_ds,
        px1_sol_ds[0] = px1_ds, px2_sol_ds[0] = px2_ds, py1_sol_ds[0] = py1_ds, py2_sol_ds[0] = py2_ds;

        for (int j = 1; j < iterations; j++)
        {
            update(lam1,lam2,dt,&x1,&x2,&y1,&y2,&px1,&px2,&py1,&py2,&H_x1_0,&H_x2_0,&H_y1_0,&H_y2_0,&x1_n,&x2_n,&y1_n,&y2_n,&H_x1_n,&H_x2_n,&H_y1_n,&H_y2_n);
            update(lam1,lam2,dt,&x1_ds,&x2_ds,&y1_ds,&y2_ds,&px1_ds,&px2_ds,&py1_ds,&py2_ds,&H_x1_0_ds,&H_x2_0_ds,&H_y1_0_ds,&H_y2_0_ds,&x1_n_ds,&x2_n_ds,&y1_n_ds,&y2_n_ds,&H_x1_n_ds,&H_x2_n_ds,&H_y1_n_ds,&H_y2_n_ds);
            
            if (j % record == 0)
            {
                x1_sol[j / record] = x1, x2_sol[j / record] = x2, y1_sol[j/ record] = y1, y2_sol[j / record] = y2;
                px1_sol[j / record] = px1, px2_sol[j / record] = px2, py1_sol[j / record] = py1, py2_sol[j / record] = py2;

                x1_sol_ds[j / record] = x1_ds, x2_sol_ds[j / record] = x2_ds, y1_sol_ds[j / record] = y1_ds, y2_sol_ds[j / record] = y2_ds;
                px1_sol_ds[j / record] = px1_ds, px2_sol_ds[j / record] = px2_ds, py1_sol_ds[j / record] = py1_ds, py2_sol_ds[j / record] = py2_ds;
            }
        }

        // Create paths and output files (using concatenation of strings allows me to change the file path for each simulation)
        std:: string path_out = "C:/Users/cilli/DIAS.cpp/Run_files4.0/run_" + std::to_string(i+1) + ".txt";
        std:: string path_out_ds = "C:/Users/cilli/DIAS.cpp/Run_files4.0/run_ds" + std::to_string(i+1) + ".txt";
        std:: fstream file_out(path_out,std::ios::out);
        std:: fstream file_out_ds(path_out_ds,std::ios::out);

        /*
        Writing all solution arrays into files to export to python for graphing / analysis.
        You have to write each array element to the file individually, as far as I'm aware, there's no way of saving
        an entire array to a file in c++ but this is almost definitely not the best and most efficient way to do it. 
        */
        for (int k=0; k < iterations / record; k++)
        {
            file_out << x1_sol[k] << ",";
            file_out_ds << x1_sol_ds[k] << ",";
        }
        for (int k=0; k < iterations / record; k++)
        {
            file_out << x2_sol[k] << ",";
            file_out_ds << x2_sol_ds[k] << ",";
        }
        for (int k=0; k < iterations / record; k++)
        {
            file_out << y1_sol[k] << ",";
            file_out_ds << y1_sol_ds[k] << ",";
        }
        for (int k=0; k < iterations / record; k++)
        {
            file_out << y2_sol[k] << ",";
            file_out_ds << y2_sol_ds[k] << ",";
        }
        for (int k=0; k < iterations / record; k++)
        {
            file_out << px1_sol[k] << ",";
            file_out_ds << px1_sol_ds[k] << ",";
        }
        for (int k=0; k < iterations / record; k++)
        {
            file_out << px2_sol[k] << ",";
            file_out_ds << px2_sol_ds[k] << ",";
        }
        for (int k=0; k < iterations / record; k++)
        {
            file_out << py1_sol[k] << ",";
            file_out_ds << py1_sol_ds[k] << ",";
        }
        for (int k=0; k < iterations / record - 1; k++)
        {
            file_out << py2_sol[k] << ",";
            file_out_ds << py2_sol_ds[k] << ",";
        }

        file_out << py2_sol[iterations / record - 1];
        file_out_ds << py2_sol_ds[iterations / record -1];
        file_out.close();
        file_out_ds.close();
    }    
    return 0;
}