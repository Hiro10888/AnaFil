//#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "matplotlibcpp.h"

#include <complex>

using namespace std;
namespace plt = matplotlibcpp;

#define DEBUG 0
int main() 
{
    // Prepare data.

    double f_start = 0.1, f_end = 10;
    int n = 100;    
    vector<double> lf(n), f(n), p(n), z(n), w(n);
    for(int i=0; i < n; ++i) {
        lf.at(i) += (double)i*log10(f_end/f_start)/n;
        f.at(i) = pow(10, lf.at(i))/f_end;
        p.at(i) = f.at(i)*f.at(i);
        z.at(i) = 0.0;
        w.at(i) = -f.at(i)*f.at(i); 

#if DEBUG == 1
	cout << "lf" << i << " = " << lf.at(i) << endl;
	cout << "f" << i << " = " << f.at(i) << endl;
	cout << "p" << i << " = " << z.at(i) << endl;
	cout << "w" << i << " = " << w.at(i) << endl;    
#endif	

    }

    plt::figure_size(1200, 780);

    // Plot a Magnitude
    // subplot2grid(nrows, ncols, row, col);
    plt::subplot2grid(2, 1, 0, 0);
    plt::named_semilogx("Magnitude",f,p);
    plt::title("ex_vector.cpp");
    plt::ylabel("Magnitude  [dB]");
    plt::ylim(0, 100);
    plt::xlim(0.1, 10.0);
    plt::legend();
    plt::grid("true");

    // Plot a Phase
    plt::subplot2grid(2, 1, 1, 0);
    plt::named_semilogx("Phase",f,w,"r");
    plt::xlabel("Frequency  [Hz]");
    plt::ylabel("Phase  [deg]");
    plt::ylim(-200, 200);
    plt::xlim(0.1, 10.0);
    plt::legend();
    plt::grid("true");

#if DEBUG == 0
    plt::show();
#endif	

#if DEBUG == 1
    // save figure
    const char* filename = "ex_vector.png";
    std::cout << "Saving result to " << filename << std::endl;;
    plt::save(filename);
#endif	
}
