#include <iostream>

using namespace std;

int main() {

    // Load in initial coordinates
    std:: fstream initial_coords("initial_coords.txt", std:: ios:: in);

    // Initialise the position and momentum variables
    double x1, p1;

    initial_coords >> x1 >> p1;
    initial_coords.close();

    // Print out other coordinates into txt file.
    cout << x1;

    return 0;
}