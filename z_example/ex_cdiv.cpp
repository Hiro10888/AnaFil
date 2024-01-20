#include <complex>
using namespace std;

complex<double> cdiv(complex<double> a, complex<double> b) {
    double db;
    complex<double> x;
    db = b.real() * b.real() + b.imag() * b.imag();
    if (db == 0.0) {
        printf(" *****  divided by 0\n");
        Error = -1;
    }
    x.real((a.real() * b.real() + a.imag() * b.imag()) / db);
    x.imag((a.imag() * b.real() - a.real() * b.imag()) / db);
    return x;
}