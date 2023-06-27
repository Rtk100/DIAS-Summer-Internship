#include <iostream>
#include <cmath>

const double delta_t = 1e-6;

int main() {
    std:: cout << -10 * pow(delta_t, 2);

    std:: cout << "Enter a radius:";

    const double pi = 3.14159;
    double radius = 0;

    std:: cin >> radius;

    std:: cout << "area of circle is: " << pi * pow(radius, 2);

    return 0;
}