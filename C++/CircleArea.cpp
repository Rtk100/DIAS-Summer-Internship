#include <iostream>
#include <cmath>

using namespace std;

int main() {

    cout << "Enter a radius:";

    const double pi = 3.14159;
    double radius = 0;

    cin >> radius;

    cout << "area of circle is: " << pi * pow(radius, 2);

    return 0;
}