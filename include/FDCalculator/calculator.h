/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 */

#ifndef INCLUDE_FDCALCULATOR_CALCULATOR_H_
#define INCLUDE_FDCALCULATOR_CALCULATOR_H_

#include <FDCalculator/biosensor_information.h>
#include <FDCalculator/utils.h>
#include <FDCalculator/constants.h>

#include <ctime>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>


namespace FDCalculator {

bool solve_explicit_slow(const struct bio_params *bio_info, void *ptr,  \
                         void (*callback_crunched)(void *, double, std::string),  \
                         std::vector<double> * I,  \
                         std::vector<double> * S,  \
                         std::vector<double> * T,  \
                         std::vector<double> * P);
bool solve_explicit_slow2(const struct biser_params *bio_info, void *ptr,  \
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
bool explicit_biseric(const struct bio_params *bio_info, void *ptr,  \
                      void (*callback_crunched)(void *, double, std::string),  \
                      std::vector<double> * I,  \
                      std::vector<double> * S,  \
                      std::vector<double> * T,  \
                      std::vector<double> * P);
std::string print(const struct bio_params *params);
std::string print(const struct biser_params *params);

bool set_initial(std::vector<double> *, double, double, double);



void saveResult(struct bio_params *params);

void setGlobalParams(struct bio_params *, unsigned int N = 100, \
                     double dt = 1e-6, resp_method type = DEFAULT_TIME, \
                     double time = 200, double ne = 1,  \
                     const std::string file = "output.dat");

bool setLocalParams(bio_params *, double, double, double, double, \
                    std::vector<std::vector<double> >);

bool setLocalParams(biser_params *, std::vector<std::vector<double> >);

bool setEqParams(proc_params * p,
                 const int n,
                 const int layers,
                 const double vmax,
                 const double km,
                 const double delta,
                 const double b0,
                 const double dt,
                 const std::vector<double> diff,
                 const std::vector<double> enz);

}  //  namespace FDCalculator

#endif  //  INCLUDE_FDCALCULATOR_CALCULATOR_H_
