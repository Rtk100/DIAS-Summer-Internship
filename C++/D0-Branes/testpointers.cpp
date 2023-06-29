#include <iostream>

void printArray(int* arr, int size) 
{   
    for (int i = 0; i < size; i++) 
    {
        int X = arr[i];
        std::cout << X << " ";
    }
    std::cout << std::endl;
}

int main() {
    int myArray[] = {1, 2, 99, 4, 5};
    int size = sizeof(myArray) / sizeof(myArray[0]);

    printArray(myArray, size);

    return 0;
}