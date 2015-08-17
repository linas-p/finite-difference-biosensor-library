#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <iterator>
#include "explicit_calculator.h"
#include "utils.h"
#include "constants.h"

using namespace std;


void calculate_explicitly(struct bio_params *bio_info, void *ptr,
                          void (*callback_crunched)(void *, int, std::string),
                          std::vector<double> & I,
                          std::vector<double> & S,
                          std::vector<double> & T,
                          std::vector<double> & P
                         )
{
    int a;

    //Srovės tankis
    double i, last_i = 0;
    // Agreguoti laikas ir srove rezultui
    std::vector<double> i_list, t_list;

    //Kintamasis rodo kaip pakito srovės tankis nuo praėjusios iteracijos
    double di;

    //Medžiagų koncentracijų masyvai
    double *current_s, *last_s;
    double *current_p, *last_p;
    // Agreguoti produkto ir substrato sprendiniai rezultatui
    std::vector<double> s_list, p_list;

    //Žingsnių pagal erdvę masyvas
    double *space_steps;

    //Tinklo taškų skaičius per visus biojutiklio sluoksnius
    int point_count;

    //Iteracija pagal laiką
    long long int t = 0;

    //Simuliacijos laikas sekundėmis
    double execution_time;

    //Kintamasis nurodo ar jau pasiektas atsako laikas
    bool response_time_reached;

    //Kintamasis nurodo ties kuriuo sluoksniu esame
    int layer;

    //Kintamasis nurodo, kad esame ties sluoksnio kraštu (ties sluoksnių sandūra)
    bool is_boundary;

    //Kinetikos dedamoji
    double kinetics_part;

    //Rezultatų saugojimui skirtas failas
    std::ofstream output_file;


    //Sukuriami lokalūs kintamieji dėl optimizavimo
    double km                    = bio_info->km;
    double dt                    = bio_info->dt;
    int n                        = bio_info->n;
    enum resp_method resp_t_meth = bio_info->resp_t_meth;
    double min_t                 = bio_info->min_t;
    double resp_t                = bio_info->resp_t;
    std::string out_file_name          = bio_info->out_file_name;
    int ne                       = bio_info->ne;
    double s0                    = bio_info->s0;
    double p0                    = bio_info->p0;
    int layer_count              = bio_info->layer_count;
    bool enz_layer;
    double Ds, Ds0, Ds1;
    double Dp, Dp0, Dp1;
    double dx, dx0, dx1;
    double v_max = bio_info->v_max_;


    //Sukuriamas rezultatų saugojimui skirtas failas
    output_file.open(out_file_name.c_str());
    output_file.close();

    //Apskaičiuojamas tinklo taškų skaičius per visus biojutiklio sluoksnius
    point_count = layer_count * n + 1;

    //Medžiagų koncentracijų masyvams išskiriama atmintis
    last_s    = new double[point_count];//
    current_s = new double[point_count];//
    last_p    = new double[point_count];//
    current_p = new double[point_count];//

    //Priskiriamos pradinės ir kai kurios kraštinės sąlygos
    fill_array(last_s, point_count - 1, 0);
    fill_array(last_p, point_count - 1, 0);
    last_s[point_count - 1] = s0;
    last_p[point_count - 1] = p0;
    current_s[point_count - 1] = s0;
    current_p[point_count - 1] = p0;
    current_p[0] = 0;

	//Pradiniai taskai
    concatenate_vals( last_s, s_list, point_count);
    concatenate_vals( last_p, p_list, point_count);
    t_list.push_back(0.);
    i_list.push_back(0.);
    //Kiekvienam sluoksniui apskaičiuojami žingsniai pagal erdvę
    space_steps = new double[layer_count];//

    for (a = 0; a < layer_count; a++) {
        space_steps[a] = bio_info->layers[a].d / n;
        dt = std::min(dt, space_steps[a]*space_steps[a]/(std::max(bio_info->layers[a].Dp,bio_info->layers[a].Ds)*2));
    }
    // Stabilumo sąlyga parenkamas "geras" dt
    if(dt != bio_info->dt) {

        if (callback_crunched != NULL) {
            std::stringstream m;
            m << "New delta t set: " << dt << " form old: " << bio_info->dt << "\n";
            callback_crunched(ptr, 0, m.str());
        }
    }
    do {
        //Iteruojama per biojutiklio sluoksnius, skaičiuojamos medžiagų koncentracijos
        layer = 0;
        //Surenkami pirmojo sluoksnio parametrai
        enz_layer = bio_info->layers[layer].enz_layer;
        Ds = bio_info->layers[layer].Ds;
        Dp = bio_info->layers[layer].Dp;
        dx = space_steps[layer];

        for (a = 1; a < point_count - 1; a++) {
            //Nustatome ar tai nėra sluoksnių sandūra
            is_boundary = !(a % n);

            //Reikšmės sluoksnių sandūrose bus skaičiuojamos vėliau pagal derinimo sąlygas
            if (is_boundary) {
                //Nustatome kuriame sluoksnyje esame
                layer++;
                //Surenkami kito sluoksnio parametrai
                enz_layer = bio_info->layers[layer].enz_layer;
                Ds = bio_info->layers[layer].Ds;
                Dp = bio_info->layers[layer].Dp;
                dx = space_steps[layer];
            } else {
                //Įskaičiuojama difuzijos įtaka
                current_s[a] = dt * Ds * \
                               (last_s[a + 1] - 2 * last_s[a] + last_s[a - 1]) / (dx * dx) + \
                               last_s[a];

                current_p[a] = dt * Dp * \
                               (last_p[a + 1] - 2 * last_p[a] + last_p[a - 1]) / (dx * dx) + \
                               last_p[a];

                //Jeigu sluoksnis yra fermentinis, tuomet prisideda ir kinetikos dalis
                if (enz_layer) {
                    kinetics_part = dt * v_max * last_s[a] / \
                                    (km + last_s[a] );
                    current_s[a] -= kinetics_part;
                    current_p[a] += kinetics_part;
                }
            }
        }

        //Sluoksnių sandūroms pritaikomos derinimo sąlygos
        for (layer = 0; layer < layer_count - 1; layer++) {
            //Apskaičiuojame kuriame taške yra layer ir layer + 1 sluoksnių sandūra
            a = n * (layer + 1);

            Ds0 = bio_info->layers[layer].Ds;
            Dp0 = bio_info->layers[layer].Dp;
            dx0 = space_steps[layer];

            Ds1 = bio_info->layers[layer + 1].Ds;
            Dp1 = bio_info->layers[layer + 1].Dp;
            dx1 = space_steps[layer + 1];

            current_s[a] = (Ds1 * dx0 * current_s[a + 1] + Ds0 * dx1 * current_s[a - 1]) / \
                           (Ds1 * dx0 + Ds0 * dx1);
            current_p[a] = (Dp1 * dx0 * current_p[a + 1] + Dp0 * dx1 * current_p[a - 1]) / \
                           (Dp1 * dx0 + Dp0 * dx1);
        }

        //Kraštinė substrato nepratekėjimo sąlyga
        current_s[0] = current_s[1];

        //Skaičiuojamas srovės tankis
        i = ne * F * bio_info->layers[0].Dp * \
            (current_p[1] - current_p[0]) / space_steps[0];
        di = fabs(i - last_i);
        last_i = i;

        //Masyvai sukeičiami vietomis
        swap_arrays(&current_s, &last_s);
        swap_arrays(&current_p, &last_p);

        //Apskaičiuojamas laikas
        t++;
        execution_time = t * dt;

        //Spausdinami rezultatai
        if ((t % INTERVAL) == 0) {
            output_file.open(out_file_name.c_str(), std::ofstream::out | std::ofstream::app);
            output_file <<  i << "," << execution_time;
            output_file << "\n";

            output_file.close();
            i_list.push_back(i);
            t_list.push_back( execution_time );
            concatenate_vals( last_s, s_list, point_count);
            concatenate_vals( last_p, p_list, point_count);

            if (callback_crunched != NULL) {
                std::stringstream m;
                if(i > 1e-30) {
                    m << "Stationarity: " << (execution_time / i) * fabs(di / dt) << "\n";
                }

                callback_crunched(ptr, (int) execution_time, m.str());
            }
        }

        //Nustatoma ar tęsti simuliaciją
        switch (resp_t_meth) {
        case MIN_TIME:
            if (execution_time < min_t) {
                response_time_reached = false;
                break;
            }
            //Jeigu jau pasiekė minimalų laiką, tuomet tikrinama pagal DEFAULT_TIME sąlygas
        case DEFAULT_TIME:
            if (i > 1e-30)
                response_time_reached = ((execution_time / i) * fabs(di / dt) <= EPSILON);
            else
                response_time_reached = false;
            break;
        case FIXED_TIME:
            response_time_reached = (execution_time >= resp_t);
            break;
        }
    } while (!response_time_reached || t < 1e5);

    //Atspausdinamas paskutinis taškas
    output_file.open(out_file_name.c_str(), std::ofstream::out | std::ofstream::app);
    output_file <<  i << "," <<  execution_time << "\n";
    output_file.close();
    if (callback_crunched != NULL)
        callback_crunched(ptr, (int) execution_time, "");


    I = i_list;
    S = s_list;
    P = p_list;
    T = t_list;

    //Atlaisvinama atmintis
    delete [] current_s;
    delete [] last_s;
    delete [] current_p;
    delete [] last_p;
    delete [] space_steps;
}
