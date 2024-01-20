#include <iostream>
#include <complex>

using namespace std;

int main()
{
    complex<double> matrix[10];
    double real = 1.0;

    for (int i = 0; i < 10; i++) {
        matrix[i] = complex<double>(0.0, 1.0);
        matrix[i].imag(matrix[i].imag() + real);
    }

    for (int i = 0; i < 10; i++) {
        cout << matrix[i] << endl;
    }

    return 0;
}