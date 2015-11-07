/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 */

#include <FDCalculator/calculator.h>

#include <string>
#include <algorithm>
#include <vector>

namespace FDCalculator {




bool solve_explicit_slow(const struct bio_params *bio_info, void *ptr,  \
                         void (*callback_crunched)(void *, double, std::string), \
                         std::vector<double> * I,  \
                         std::vector<double> * S,  \
                         std::vector<double> * T,  \
                         std::vector<double> * P) {
    //  Kintamasis nurodo ties kuriuo sluoksniu esame
    int layer, a;

    //  Srovės tankis
    double i, last_i = 0;

    //  Agreguoti laikas ir srove rezultui
    std::vector<double> i_list, t_list;
    i_list.reserve(500000);
    t_list.reserve(500000);

    //  Kintamasis rodo kaip pakito srovės tankis nuo praėjusios iteracijos
    double di;

    //  Agreguoti produkto ir substrato sprendiniai rezultatui
    std::vector<double> s_list, p_list;

    //  Tinklo taškų skaičius per visus biojutiklio sluoksnius
    int point_count;

    //  Iteracija pagal laiką
    unsigned int t = 0;

    //  Simuliacijos laikas sekundėmis
    double execution_time;

    //  Kintamasis nurodo ar jau pasiektas atsako laikas
    bool response_time_reached;

    //  Kintamasis nurodo, kad esame ties sluoksnio kraštu (sluoksnių sandūra)
    bool is_boundary;

    //  Rezultatų saugojimui skirtas failas
    std::ofstream output_file;

    //  Sukuriami lokalūs kintamieji dėl optimizavimo
    double km                    = bio_info->km;
    double dt                    = bio_info->dt;
    int n                        = bio_info->n;
    enum resp_method resp_t_meth = bio_info->resp_t_meth;
    double min_t                 = bio_info->cond_t;
    double resp_t                = bio_info->cond_t;
    std::string out_file_name          = bio_info->out_file_name;
    double s0                    = bio_info->s0;
    double p0                    = bio_info->p0;
    int layer_count              = bio_info->layer_count;
    double enz_layer;
    double Ds0, Ds1;
    double Dp0, Dp1;
    double dx0, dx1;
    double dx_sq, dtDs, dtDp;

    double i_w = bio_info->ne * F * bio_info->layers[0].Dp / \
                 bio_info->layers[0].dx;


    //  Sukuriamas rezultatų saugojimui skirtas failas
    output_file.open(out_file_name.c_str());
    //  output_file.close();

    //  Apskaičiuojamas tinklo taškų skaičius visiem biojutiklio sluoksniam
    point_count = layer_count * n + 1;

    //  Medžiagų koncentracijų masyvai
    std::vector<double> current_s(point_count, 0.), last_s(point_count, 0.);
    std::vector<double> current_p(point_count, 0.), last_p(point_count, 0.);
    std::vector<double> last_kinetics(point_count, 0.);

    //  Priskiriamos pradinės ir kai kurios kraštinės sąlygos
    if (unlikely(!(set_initial(&last_s, 0, 0, s0)  \
                   && set_initial(&last_p, 0, 0, p0)  \
                   && set_initial(&current_p, 0, 0, s0) \
                   && set_initial(&current_s, 0, 0, p0)))) {
        return false;
    }


    //  Pradiniai taskai
    S->insert(S->end(), last_s.begin(), last_s.end());
    P->insert(P->end(), last_p.begin(), last_p.end());
    t_list.push_back(0.);
    i_list.push_back(0.);

    double vmaxt = bio_info->vmax*dt;

    std::clock_t start;
    start = std::clock();
    do {
        //  Iteruojama per biojutiklio sluoksnius,
        //  skaičiuojamos medžiagų koncentracijos
        layer = 0;
        //  Surenkami pirmojo sluoksnio parametrai
        enz_layer = bio_info->layers[layer].enz_layer;
        dtDs = bio_info->layers[layer].Ds*dt;
        dtDp = bio_info->layers[layer].Dp*dt;
        dx_sq = pow(bio_info->layers[layer].dx, 2);


        for (a = 1; a < point_count - 1; a++) {
            last_kinetics[a] = vmaxt * last_s[a] /  \
                               (km + last_s[a]);
        }

        for (a = 1; a < point_count - 1; a++) {
            //  Nustatome ar tai nėra sluoksnių sandūra
            is_boundary = !(a % n);

            //  Reikšmės sluoksnių sandūrose bus skaičiuojamos vėliau
            if (unlikely(is_boundary)) {
                //  Nustatome kuriame sluoksnyje esame
                layer++;
                //  Surenkami kito sluoksnio parametrai
                enz_layer = bio_info->layers[layer].enz_layer;
                dtDs = bio_info->layers[layer].Ds*dt;
                dtDp = bio_info->layers[layer].Dp*dt;
                dx_sq = pow(bio_info->layers[layer].dx, 2);
            } else {
                //  Įskaičiuojama difuzijos įtaka
                current_s[a] =  dtDs *  \
                                (last_s[a + 1] - 2 * last_s[a] + last_s[a - 1]) / (dx_sq) +  \
                                last_s[a] - enz_layer*last_kinetics[a];
                //  no branching with kinetics

                current_p[a] =   dtDp *  \
                                 (last_p[a + 1] - 2 * last_p[a] + last_p[a - 1]) / (dx_sq) +  \
                                 last_p[a] + enz_layer*last_kinetics[a];
            }
        }

        //  Sluoksnių sandūroms pritaikomos derinimo sąlygos
        for (layer = 0; layer < layer_count - 1; layer++) {
            //  Apskaičiuojame kuriame taške yra layer ir layer + 1 sandūro
            a = n * (layer + 1);

            Ds0 = bio_info->layers[layer].Ds;
            Dp0 = bio_info->layers[layer].Dp;
            dx0 = bio_info->layers[layer].dx;

            Ds1 = bio_info->layers[layer + 1].Ds;
            Dp1 = bio_info->layers[layer + 1].Dp;
            dx1 = bio_info->layers[layer + 1].dx;

            current_s[a] = (Ds1 * dx0 * current_s[a + 1] +  \
                            Ds0 * dx1 * current_s[a - 1]) /  \
                           (Ds1 * dx0 + Ds0 * dx1);
            current_p[a] = (Dp1 * dx0 * current_p[a + 1] +  \
                            Dp0 * dx1 * current_p[a - 1]) /  \
                           (Dp1 * dx0 + Dp0 * dx1);
        }

        //  Kraštinė substrato nepratekėjimo sąlyga
        current_s[0] = current_s[1];
        current_s[point_count-1] = s0;
        current_p[0] = 0.;
        current_p[point_count-1] = 0.;
        //  Skaičiuojamas srovės tankis
        i = i_w * (current_p[1] - current_p[0]);
        di = fabs(i - last_i);
        last_i = i;

        //  Masyvai sukeičiami vietomis
        last_s = current_s;
        last_p = current_p;

        //  Apskaičiuojamas laikas
        t++;
        execution_time = t * dt;

        //  Spausdinami rezultatai
        if ((t % INTERVAL) == 0) {
            //  output_file.open(out_file_name.c_str(),
            //  std::ofstream::out | std::ofstream::app);
            output_file <<  i << "," << execution_time << "\n";
            //  output_file.close();

            i_list.push_back(i);
            t_list.push_back(execution_time);
            S->insert(S->end(), last_s.begin(), last_s.end());
            P->insert(P->end(), last_p.begin(), last_p.end());

            /*if (callback_crunched != NULL) {
                std::stringstream m;
                if(i > 1e-30) {
                    m << "Stationarity: " << (execution_time / i) * fabs(di / dt) << "\n";
                }

                callback_crunched(ptr, execution_time, m.str());
            }*/
        }

        //  Nustatoma ar tęsti simuliaciją
        //  switch (resp_t_meth) {
        switch (__builtin_expect(resp_t_meth, DEFAULT_TIME)) {
        case MIN_TIME:
            if (execution_time < min_t) {
                response_time_reached = false;
                break;
            }
            //  Jeigu jau pasiekė minimalų laiką,
            //  tuomet tikrinama pagal DEFAULT_TIME sąlygas
        case DEFAULT_TIME:
            if (i > 1e-30)
                response_time_reached = ((execution_time / i) *  \
                                         fabs(di / dt) <= EPSILON);
            else
                response_time_reached = false;
            break;
        case FIXED_TIME:
            response_time_reached = (execution_time >= resp_t);
            break;
        }
    } while (likely(!response_time_reached));

    std::cout << "t " << t << " " << execution_time << std::endl;
    double tt = (std::clock()-start)/ static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "loop time " << tt << std::endl;

    //  Atspausdinamas paskutinis taškas
    //  output_file.open(out_file_name.c_str(),
    //  std::ofstream::out | std::ofstream::app);
    output_file <<  i << "," <<  execution_time << "\n";
    output_file.close();
    if (callback_crunched != NULL)
        callback_crunched(ptr, execution_time, "");

    *I = i_list;
    *T = t_list;

    return true;
}


bool solve_explicit_slow2(const struct biser_params *bio_info, void *ptr,  \
                          void (*callback_crunched)(void *, double, std::string),  \
                          std::vector<double> * I,  \
                          std::vector<double> * S,  \
                          std::vector<double> * T,  \
                          std::vector<double> * P) {
    //  Kintamasis nurodo ties kuriuo sluoksniu esame
    int layer, a;

    //  Srovės tankis
    double i, last_i = 0;

    //  Agreguoti laikas ir srove rezultui
    std::vector<double> i_list, t_list;
    i_list.reserve(500000); t_list.reserve(500000);

    //  Kintamasis rodo kaip pakito srovės tankis nuo praėjusios iteracijos
    double di;

    //  Agreguoti produkto ir substrato sprendiniai rezultatui
    std::vector<double> s_list, p_list;

    //  Tinklo taškų skaičius per visus biojutiklio sluoksnius
    int point_count;

    //  Iteracija pagal laiką
    unsigned int t = 0;

    //  Simuliacijos laikas sekundėmis
    double execution_time;

    //  Kintamasis nurodo ar jau pasiektas atsako laikas
    bool response_time_reached;

    //  Kintamasis nurodo, kad esame ties sluoksnio kraštu (sluoksnių sandūra)
    bool is_boundary;

    //  Rezultatų saugojimui skirtas failas
    std::ofstream output_file;

    //  Sukuriami lokalūs kintamieji dėl optimizavimo
//    double km                    = bio_info->equations[0].km;
    double dt                    = bio_info->dt;
    int n                        = bio_info->n;
    enum resp_method resp_t_meth = bio_info->resp_t_meth;
    double min_t                 = bio_info->cond_t;
    double resp_t                = bio_info->cond_t;
    std::string out_file_name          = bio_info->out_file_name;
    double s0                    = bio_info->equations[0].b0;
    double p0                    = bio_info->equations[1].b0;
    int layer_count              = bio_info->layer_count;
    double enz_layer;
    double Ds0, Ds1;
    double Dp0, Dp1;
    double dx0, dx1;
    double dx_sq, dtDs, dtDp;

    double i_w = bio_info->ne * F * bio_info->equations[1].D[0] / \
                 bio_info->equations[1].dx;


    //  Sukuriamas rezultatų saugojimui skirtas failas
    output_file.open(out_file_name.c_str());
    //  output_file.close();

    //  Apskaičiuojamas tinklo taškų skaičius visiem biojutiklio sluoksniam
    point_count = layer_count * n + 1;

    //  Medžiagų koncentracijų masyvai

    //  Priskiriamos pradinės ir kai kurios kraštinės sąlygos
    if (unlikely(!(set_initial(&bio_info->equations[0].last, 0, 0, s0)  \
                   && set_initial(&bio_info->equations[1].last, 0, 0, p0)  \
                   && set_initial(&bio_info->equations[0].current, 0, 0, s0) \
                   && set_initial(&bio_info->equations[1].current, 0, 0, p0)))) {
        return false;
    }


    //  Pradiniai taskai
    S->insert(S->end(), bio_info->equations[0].last.begin(), \
              bio_info->equations[0].last.end());
    P->insert(P->end(), bio_info->equations[1].last.begin(), \
              bio_info->equations[1].last.end());
    t_list.push_back(0.);
    i_list.push_back(0.);

    std::clock_t start;
    start = std::clock();
    do {
        //  Iteruojama per sluoksnius, skaičiuojamos medžiagų koncentracijos
        layer = 0;
        //  Surenkami pirmojo sluoksnio parametrai
        enz_layer = bio_info->equations[0].enz_layer[layer];
        dtDs = bio_info->equations[0].D[layer]*dt;
        dtDp = bio_info->equations[1].D[layer]*dt;
        dx_sq = pow(bio_info->equations[0].dx, 2);

        bio_info->equations[0].update_kinetics();
        /*for (a = 1; a < point_count - 1; a++) {
            bio_info->equations[0].kinetics[a] = \
            vmaxt * bio_info->equations[0].last[a] /  \
            (km + bio_info->equations[0].last[a]);
        }*/

        for (a = 1; a < point_count - 1; a++) {
            //  Nustatome ar tai nėra sluoksnių sandūra
            is_boundary = !(a % n);

            //  Reikšmės sluoksnių sandūrose bus skaičiuojamos vėliau
            if (unlikely(is_boundary)) {
                //  Nustatome kuriame sluoksnyje esame
                layer++;
                //  Surenkami kito sluoksnio parametrai
                enz_layer = bio_info->equations[0].enz_layer[layer];
                dtDs = bio_info->equations[0].D[layer]*dt;
                dtDp = bio_info->equations[1].D[layer]*dt;
                dx_sq = pow(bio_info->equations[0].dx, 2);
            } else {
                //  Įskaičiuojama difuzijos įtaka
                bio_info->equations[0].current[a] =  dtDs *  \
                                                     (bio_info->equations[0].last[a + 1] - 2 * bio_info->equations[0].last[a] + bio_info->equations[0].last[a - 1]) / (dx_sq) +  \
                                                     bio_info->equations[0].last[a] - enz_layer*bio_info->equations[0].kinetics[a];
                //  no branching with kinetics

                bio_info->equations[1].current[a] =   dtDp *  \
                                                      (bio_info->equations[1].last[a + 1] - 2 * bio_info->equations[1].last[a] + bio_info->equations[1].last[a - 1]) / (dx_sq) +  \
                                                      bio_info->equations[1].last[a] + enz_layer*bio_info->equations[0].kinetics[a];
            }
        }

        //  Sluoksnių sandūroms pritaikomos derinimo sąlygos
        for (layer = 0; layer < layer_count - 1; layer++) {
            //  Apskaičiuojame kuriame taške yra layer ir layer + 1 sandūro
            a = n * (layer + 1);

            Ds0 = bio_info->equations[0].D[layer];
            Dp0 = bio_info->equations[1].D[layer];
            dx0 = bio_info->equations[0].dx;

            Ds1 = bio_info->equations[0].D[layer + 1];
            Dp1 = bio_info->equations[1].D[layer + 1];
            dx1 = bio_info->equations[0].dx;

            bio_info->equations[0].current[a] = (Ds1 * dx0 * bio_info->equations[0].current[a + 1] +  \
                                                 Ds0 * dx1 * bio_info->equations[0].current[a - 1]) /  \
                                                (Ds1 * dx0 + Ds0 * dx1);
            bio_info->equations[1].current[a] = (Dp1 * dx0 * bio_info->equations[1].current[a + 1] +  \
                                                 Dp0 * dx1 * bio_info->equations[1].current[a - 1]) /  \
                                                (Dp1 * dx0 + Dp0 * dx1);
        }

        //  Kraštinė substrato nepratekėjimo sąlyga
        bio_info->equations[0].current[0] = bio_info->equations[0].current[1];
        bio_info->equations[0].current[point_count-1] = s0;
        bio_info->equations[1].current[0] = 0.;
        bio_info->equations[1].current[point_count-1] = 0.;
        //  Skaičiuojamas srovės tankis
        i = i_w * (bio_info->equations[1].current[1] - \
                   bio_info->equations[1].current[0]);
        di = fabs(i - last_i);
        last_i = i;

        //  Masyvai sukeičiami vietomis
        bio_info->equations[0].next();
        bio_info->equations[1].next();


        //  Apskaičiuojamas laikas
        t++;
        execution_time = t * dt;

        //  Spausdinami rezultatai
        if ((t % INTERVAL) == 0) {
            //  output_file.open(out_file_name.c_str(),
            //  std::ofstream::out | std::ofstream::app);
            output_file <<  i << "," << execution_time << "\n";
            //  output_file.close();

            i_list.push_back(i);
            t_list.push_back(execution_time);
            S->insert(S->end(), bio_info->equations[0].last.begin(), \
                      bio_info->equations[0].last.end());
            P->insert(P->end(), bio_info->equations[1].last.begin(), \
                      bio_info->equations[1].last.end());

            /*if (callback_crunched != NULL) {
                std::stringstream m;
                if(i > 1e-30) {
                    m << "Stationarity: " << \
                    (execution_time / i) * fabs(di / dt) << "\n";
                }

                callback_crunched(ptr, execution_time, m.str());
            }*/
        }

        //  Nustatoma ar tęsti simuliaciją
        //  switch (resp_t_meth) {
        /*switch (__builtin_expect(resp_t_meth, DEFAULT_TIME)) {
        case MIN_TIME:
            if (execution_time < min_t) {
                response_time_reached = false;
                break;
            }
            //  Jeigu jau pasiekė minimalų laiką,
            //  tuomet tikrinama pagal DEFAULT_TIME sąlygas
        case DEFAULT_TIME:
            if (i > 1e-30)
                response_time_reached = ((execution_time / i) *  \
                                         fabs(di / dt) <= EPSILON);
            else
                response_time_reached = false;
            break;
        case FIXED_TIME:
            response_time_reached = (execution_time >= resp_t);
            break;
        }*/

        if (t == 336213) {
            response_time_reached = true;
        } else {
            response_time_reached = false;
        }
    } while (likely(!response_time_reached));

    std::cout << "t " << t << " " << execution_time << std::endl;
    double tt = (std::clock()-start)/ static_cast<double>(CLOCKS_PER_SEC);
    std::cout << "loop time " << tt << std::endl;

    //  Atspausdinamas paskutinis taškas
    //  output_file.open(out_file_name.c_str(),
    //  std::ofstream::out | std::ofstream::app);
    output_file <<  i << "," <<  execution_time << "\n";
    output_file.close();
    if (callback_crunched != NULL)
        callback_crunched(ptr, execution_time, "");

    *I = i_list;
    *T = t_list;

    return true;
}



bool set_initial(std::vector<double> * grid, double region, double boundary_1, \
                 double boundary_2) {
    int n = grid->size() - 1;
    if (n < 1) {
        return false;
    }

    (*grid)[0] = boundary_1;
    std::fill(grid->begin()+1, grid->end()-1, region);
    (*grid)[n] = boundary_2;
    return true;
}

}  //  namespace FDCalculator
