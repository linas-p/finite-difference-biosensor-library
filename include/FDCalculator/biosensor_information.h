/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 * */

#ifndef INCLUDE_FDCALCULATOR_BIOSENSOR_INFORMATION_H_
#define INCLUDE_FDCALCULATOR_BIOSENSOR_INFORMATION_H_
#include <string>

namespace FDCalculator {

enum resp_method {
    DEFAULT_TIME,  //  0 - iki pusiausvyros
    MIN_TIME,      //  1 - iki pusiausvyros su nurodytu minimaliu laiku
    FIXED_TIME     //  2 - fiksuotas laikas
};

struct layer_params {
    //  Difuzijos koeficientai (m^2/s)
    double Ds;
    double Dp;
    //  Žingsnis pagal x (m)
    double dx;
    //  Laukas nurodo ar tai fermento sluoksnis
    bool enz_layer;
};

struct bio_params {
    //  Biojutiklio sluoksnių masyvas
    struct layer_params *layers;
    //  Išvedimo failas
    std::string out_file_name;

    //  Žingsnis pagal laiką (s), dažniau parenkmas pagal stabilumo sąlygą.
    double dt;
    //  atsako laikas (s) (fiksuotas, minimalus)
    double cond_t;
    //  Substrato koncentracija tirpale (mol/m^3)
    double s0;
    //  Produkto koncentracija tirpale (mol/m^3)
    double p0;
    //  Pusiausvyros konstantos (mol/m^3)
    double km;
    //  Maksimalus pasiekiamas greitis
    double vmax;

    //  Į kiek dalių dalinami sluoksniai.
    unsigned int n;
    //  Biojutiklio sluoksnių skaičius
    unsigned int layer_count;
    //  Elektronų, dalyvaujančių krūvio pernešime, skaičius
    int ne;
    //  Metodas, kuriuo bus nustatomas atsako laikas.
    enum resp_method resp_t_meth;
};
}  // namespace FDCalculator

#endif  //  INCLUDE_FDCALCULATOR_BIOSENSOR_INFORMATION_H_
