#include <iostream>

double acceleration(int j, double h)
{   
    double F;
    F += h;
    F += j;
    return F;
}
int main() {

    for (int j = 0;j < 3 ; ++j)
    {
        std::cout << "j"<<j<<'\n'<<"F" << acceleration(j, -1.0);

    }

    return 0;
}


