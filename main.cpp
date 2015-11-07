/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 */

#include <FDCalculator/calculator.h>
#include <FDCalculator/biosensor_information.h>

#include <ctime>
#include <vector>
#include <string>
#include <iostream>



using FDCalculator::bio_params;
using FDCalculator::proc_params;
using FDCalculator::resp_method;
using FDCalculator::biser_params;
using FDCalculator::solve_explicit;
using FDCalculator::solve_explicit_slow;
using FDCalculator::solve_explicit_slow2;
using FDCalculator::setEqParams;
using FDCalculator::setGlobalParams;

void callback_crunched(void *ptr, double time, std::string info) {
    if (!info.empty()) {
        std::cout << info.c_str();
    }
    std::cout << "simulated "<< time << " (s) \n" << std::endl;
}
void fill(biser_params * b) {
    b->out_file_name = "test.dat";
    b->dt = 8.33333e-10;
    b->cond_t = 1;
    b->n = 100;
    b->layer_count = 2;
    b->ne = 1;
    b->resp_t_meth = resp_method::DEFAULT_TIME;
    b->equations = new proc_params[b->layer_count];
    std::vector<double> D, e;
    D.push_back(300*1e-6);
    D.push_back(600*1e-6);
    e.push_back(1.);
    e.push_back(0.);
    setEqParams(&b->equations[0], b->n, b->layer_count, 10*1e-3, 100*1e-3, \
                100*1e-6, 100*1e-3, b->dt, D, e);
    setEqParams(&b->equations[1], b->n, b->layer_count, 10*1e-3, 100*1e-3, \
                100*1e-6, 100*1e-3, b->dt, D, e);
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
    std::cout << print(p) << std::endl;

    solve_explicit(p, NULL, &callback_crunched, &t, &t, &t, &t);
    solve_explicit_slow(p, NULL, &callback_crunched, &t, &t, &t, &t);

    biser_params * b = new biser_params;
    fill(b);
    std::cout << print(b) << std::endl;
    solve_explicit_slow2(b, NULL, &callback_crunched, &t, &t, &t, &t);
}
