#ifndef BIOSENSORINFORMATION_H
#define BIOSENSORINFORMATION_H
#include <string>

enum resp_method
{
    DEFAULT_TIME, //0 - iki pusiausvyros
    MIN_TIME,     //1 - iki pusiausvyros su nurodytu minimaliu laiku
    FIXED_TIME    //2 - fiksuotas laikas
};

struct layer_params
{
    //Laukas nurodo ar tai fermento sluoksnis
    bool enz_layer;

    //Difuzijos koeficientai (m^2/s)
    double Ds;
    double Dp;

    //Sluoksnio storis (m)
    double d;
};

struct bio_params
{
    //Pusiausvyros konstantos (mol/m^3)
    double km;

    double v_max_;
    //Žingsnis pagal laiką (s), dažniau parenkmas pagal stabilumo sąlygą
    double dt;
    //Į kiek dalių dalinami sluoksniai
    int n;
    //Metodas, kuriuo bus nustatomas atsako laikas:
    enum resp_method resp_t_meth;
    //Minimalus atsako laikas (s)
    double min_t;
    //Fiksuotas atsako laikas (s)
    double resp_t;
    //Išvedimo failas
    std::string out_file_name;
    //Elektronų, dalyvaujančių krūvio pernešime, skaičius
    int ne;
    //Substrato koncentracija tirpale (mol/m^3)
    double s0;
    //Produkto koncentracija tirpale (mol/m^3)
    double p0;
    //Biojutiklio sluoksnių skaičius
    int layer_count;
    //Biojutiklio sluoksnių masyvas
    struct layer_params *layers;
};

#endif // BIOSENSORINFORMATION_H
