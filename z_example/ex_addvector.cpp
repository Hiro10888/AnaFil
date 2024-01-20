#include <iostream>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;

int main() {
    int n = 100;
    vector<double> p(n), q(n), r(n);
    for(int i = 0; i < n; ++ i){
        p[i] = i;
        q[i] = i;
        r[i] = p[i] + q[i];
    }
    plt::plot(p, q);
    plt::plot(p, r);
    plt::show();
    return 0;
}


/*
int main() {
    int n = 100;
    vector<double> x(n), y(n), z(n);
    for(int i = 0; i < n; ++ i){
        x[i] = i;
        y[i] = i;
        z[i] = x[i] + y[i];
    }
    plt::plot(x, y);
    plt::plot(x, z);
    plt::show();
    return 0;
}
*/