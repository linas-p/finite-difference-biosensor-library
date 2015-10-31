/* Copyright 2015
 * Linas Petkevicius
 * Vilnius University
 * GNU General Public License
 * */

#include <gtest/gtest.h>

#include <iostream>
#include <ctime>
#include <vector>
#include <string>

#include <FDCalculator/calculator.h>
#include <FDCalculator/biosensor_information.h>

using FDCalculator::bio_params;
using FDCalculator::solve_explicit;
using FDCalculator::setGlobalParams;

void callback_crunched(void *ptr, double time, std::string info) {
    if (!info.empty()) {
        std::cout << info.c_str();
    }
    std::cout << "simulated "<< time << " (s) \n" << std::endl;
}

TEST(fake, iliustration) {
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
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
