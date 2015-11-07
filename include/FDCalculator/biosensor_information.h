/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 * */

#ifndef INCLUDE_FDCALCULATOR_BIOSENSOR_INFORMATION_H_
#define INCLUDE_FDCALCULATOR_BIOSENSOR_INFORMATION_H_
#include <string>
#include <vector>

namespace FDCalculator {

const int MAX_LAYERS = 4;

enum resp_method {
    DEFAULT_TIME,  //  0 - iki pusiausvyros
    MIN_TIME,      //  1 - iki pusiausvyros su nurodytu minimaliu laiku
    FIXED_TIME     //  2 - fiksuotas laikas
};

struct proc_params {  // For each equation
    std::vector<double> current;
    std::vector<double> last;
    std::vector<double> kinetics;

    //  Difuzijos koeficientai (m^2/s)
    double D[MAX_LAYERS];
    //  Laukas nurodo ar tai fermento sluoksnis
    double enz_layer[MAX_LAYERS];

    //  Maksimalus pasiekiamas greitis
    double vmax;
    //  Pusiausvyros konstantos (mol/m^3)
    double km;
    //  Žingsnis pagal x (m)
    double dx;
    // Boundary value
    double b0;
    double vmaxdt;
    std::vector<double>::iterator i, j;
    void init(const int layers, const int grid_sz) {
        last.resize(grid_sz);
        current.resize(grid_sz);
        kinetics.resize(grid_sz);
    }
    proc_params() {}
    ~proc_params() {
    }

    void next() {
        last = current;
    }
    void update_kinetics() {
        for (i = last.begin(), j = kinetics.begin(); i != last.end() && \
                j !=kinetics.end(); ++i, ++j) {
            *j = vmaxdt*(*i)*(km+*i);
        }
    }
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

struct biser_params {
    //  Biojutiklio sluoksnių masyvas
    struct proc_params * equations;
    //  Išvedimo failas
    std::string out_file_name;

    //  Žingsnis pagal laiką (s), dažniau parenkmas pagal stabilumo sąlygą.
    double dt;
    //  atsako laikas (s) (fiksuotas, minimalus)
    double cond_t;
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
