/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 */

#ifndef INCLUDE_FDCALCULATOR_CALCULATOR_H_
#define INCLUDE_FDCALCULATOR_CALCULATOR_H_



#include <ctime>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>

#include <FDCalculator/biosensor_information.h>
#include <FDCalculator/utils.h>
#include <FDCalculator/constants.h>

namespace FDCalculator {
void solve_explicit_slow(struct bio_params *bio_info, void *ptr,  \
                    void (*callback_crunched)(void *, double, std::string),  \
                    std::vector<double> * I,  \
                    std::vector<double> * S,  \
                    std::vector<double> * T,  \
                    std::vector<double> * P);
void solve_explicit(struct bio_params *bio_info, void *ptr,  \
                    void (*callback_crunched)(void *, double, std::string),  \
                    std::vector<double> * I,  \
                    std::vector<double> * S,  \
                    std::vector<double> * T,  \
                    std::vector<double> * P);
std::string print(struct bio_params *params);

void saveResult(struct bio_params *params);

void setGlobalParams(bio_params *, unsigned int N = 100, double dt = 1e-6, \
                     resp_method type = DEFAULT_TIME, double time = 200, double ne = 1,  \
                     const std::string file = "output.dat");

void setLocalParams(bio_params *, double, double, double, double, std::vector<std::vector<double> >);
}  //  namespace FDCalculator

#endif  //  INCLUDE_FDCALCULATOR_CALCULATOR_H_
