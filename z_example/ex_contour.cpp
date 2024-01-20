#include <iostream>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using namespace std;

int main() {
    int n = 100;
    vector<double> x(n), y(n);
    for(int i = 0; i < n; ++ i){
        x[i] = i;
        y[i] = i;
    }
    plt::contour(x, y);
    plt::show();
    return 0;
}