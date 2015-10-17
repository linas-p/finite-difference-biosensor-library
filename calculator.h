#ifndef CALCULATOR_H
#define CALCULATOR_H
#include "biosensor_information.h"
#include "utils.h"
#include "constants.h"
#include <cmath>
#include <ctime>
#include <cfloat>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

class Calculator {

private:


    //Srovė
    double i_, last_i_;

    //Rezultatų failas
    std::ofstream output_file_;

public:

    std::vector<double> i_rez_, t_rez_;
	std::vector<double> s_rez_, p_rez_;

    //Sistemos parametrai
    bio_params *params_;

    Calculator();
    ~Calculator();

    //Parametru iniciavimas
    void setGlobalParams(unsigned int N = 100, double dt = 1e-6, resp_method type = DEFAULT_TIME, double time = 200, double ne = 1, const std::string file = "output.dat");
    void setLocalParams(double, double, double, double, std::vector<std::vector<double> >);

    //Pradiniu salygu spausdinimas
    std::string print(struct bio_params *params);

    static void solve(struct bio_params *bio_info, void *ptr,
                          void (*callback_crunched)(void *, double, std::string),
                          std::vector<double> & I,
                          std::vector<double> & S,
                          std::vector<double> & T,
                          std::vector<double> & P);
    void saveResult();
};


#endif
