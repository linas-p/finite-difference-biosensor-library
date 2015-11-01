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
using FDCalculator::solve_explicit_slow;
using FDCalculator::setGlobalParams;
using FDCalculator::set_initial;

void callback_crunched(void *ptr, double time, std::string info) {
    if (!info.empty()) {
        std::cout << info.c_str();
    }
    std::cout << "simulated "<< time << " (s) \n" << std::endl;
}

TEST(set_initial, empty) {
    std::vector<double> a;
    EXPECT_EQ(set_initial(&a, 1, 1, 1), false);
    a.push_back(1.);
    EXPECT_EQ(set_initial(&a, 1, 1, 1), false);
    a.push_back(2.);
    EXPECT_EQ(set_initial(&a, 1, 1, 1), true);
    a.push_back(3.);
    EXPECT_EQ(set_initial(&a, 0, 2, 3), true);
    EXPECT_EQ(a[0], 2);
    EXPECT_EQ(a[a.size()-1], 3);

}

TEST(set_initial, normal_cases){
	std::vector<double> a= {1,2,3};
    EXPECT_EQ(set_initial(&a, 0, 3, 1), true);
    EXPECT_EQ(a[0], 3);
    EXPECT_EQ(a[1], 0);
    EXPECT_EQ(a[2], 1);

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
    solve_explicit_slow(p, NULL, &callback_crunched, &t, &t, &t, &t);
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
