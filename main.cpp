#include <iostream>
#include <ctime>
#include "calculator.h"

int main() {

    Calculator cal;
    cal.setGlobalParams(127, 1e-6);
    std::vector< std::vector< double > > layer_params;
    std::vector<double> lay1, lay2;
    lay1.push_back(300*1e-6);
    lay1.push_back(300*1e-6);
    lay1.push_back(100*1e-6);
    lay1.push_back(1);
    lay2.push_back(600*1e-6);
    lay2.push_back(600*1e-6);
    lay2.push_back(100*1e-6);
    lay2.push_back(0);
    layer_params.push_back(lay1);
    layer_params.push_back(lay2);
    cal.setLocalParams(100*1e-3, 0., 100*1e-3, 10*1e-3, layer_params);
    std::cout << cal.print() << std::endl;
    cal.solve(cal.params_);
    cal.saveResult();

}
