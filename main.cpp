/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 */

#include <iostream>
#include <ctime>
#include <vector>
#include <string>

#include <FDCalculator/calculator.h>
#include <FDCalculator/biosensor_information.h>

using FDCalculator::bio_params;
using FDCalculator::solve_explicit;
using FDCalculator::solve_explicit_slow;
using FDCalculator::setGlobalParams;

void callback_crunched(void *ptr, double time, std::string info) {
    if (!info.empty()) {
        std::cout << info.c_str();
    }
    std::cout << "simulated "<< time << " (s) \n" << std::endl;
}

int main() {
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

    std::vector<double> t;
    bio_params * p = new bio_params;
    setGlobalParams(p);
    setLocalParams(p, 100*1e-3, 0., 100*1e-3, 10*1e-3, layer_params);
    print(p);

    solve_explicit(p, NULL, &callback_crunched, &t, &t, &t, &t);
    solve_explicit_slow(p, NULL, &callback_crunched, &t, &t, &t, &t);
}
