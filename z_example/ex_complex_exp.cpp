#define _USE_MATH_DEFINES
#define VMAX 100

#include "matplotlibcpp.h"
#include <bits/stdc++.h>
// https://qiita.com/e869120/items/518297c6816adb67f9a5
//#include <iostream>
//#include <complex>
//#include <fstream>

using namespace std;
namespace plt = matplotlibcpp;

int main() 
{
    vector<double> t(VMAX), re(VMAX), im(VMAX);
    ofstream ofs("output.dat");

    const double pi = acos(-1);
    const complex<double> i(0, 1);
    int N = 0;

    for (double x = -pi; x <= pi; x += 2*pi/t.size()) {     
        t[N] = x;
        re[N] = exp(i * x).real();
        im[N] = exp(i * x).imag();
        ofs << t[N] << " " << re[N]  << " " << im[N] << std::endl;
        N += 1;
    }

    plt::plot(t, re);
    plt::plot(t, im);          
    plt::show();
}