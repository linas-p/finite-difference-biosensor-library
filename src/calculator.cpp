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



void saveResult(struct bio_params *params) {
    /*output_file_.open(params->out_file_name.c_str(), std::ofstream::out | std::ofstream::app);
    for(unsigned int i = 0; i < i_rez_.size(); i++) {
        output_file_ << t_rez_[i] << "," << i_rez_[i] << "\n";
    }
    output_file_.close();*/
}
std::string print(struct bio_params *params) {
    std::stringstream m;
    m << "Global conditions:\n  n = " <<  params->n << " dt = "  \
      << params->dt << " type " << params->resp_t_meth   \
      << " ne = " << params->ne << " output goes to: "  \
      << params->out_file_name << "\n";

    m << "Initial conditions: \n  s0= " << params->s0  \
      << " p0= " << params->p0 << " km= " << params->km  \
      << " vmax= " << params->vmax  \
      << " number of layers= " << params->layer_count << " \n";

    for (unsigned int i = 0; i < params->layer_count; i++) {
        m << "  layer " << i << " params: Ds= " << params->layers[i].Ds <<  \
          " Dp= " << params->layers[i].Dp << " d= "  \
          << params->layers[i].dx*params->n << " with accuracy: "  \
          << params->layers[i].dx << " is enzyme? "  \
          << params->layers[i].enz_layer << "\n";
    }
    return m.str();
}

void solve_explicit(struct bio_params *bio_info, void *ptr,  \
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

    //  Kintamasis rodo kaip pakito srovės tankis nuo praėjusios iteracijos
    double di;

    //  Medžiagų koncentracijų masyvai
    double *current_s, *last_s;
    double *current_p, *last_p;

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

    //  Kinetikos dedamoji
    double kinetics_part;

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
    bool enz_layer;
    double Ds, Ds0, Ds1;
    double Dp, Dp0, Dp1;
    double dx, dx0, dx1;

    double i_w = bio_info->ne * F * bio_info->layers[0].Dp / \
                 bio_info->layers[0].dx;


    //  Sukuriamas rezultatų saugojimui skirtas failas
    output_file.open(out_file_name.c_str());
    //  output_file.close();

    //  Apskaičiuojamas tinklo taškų skaičius visiem biojutiklio sluoksniam
    point_count = layer_count * n + 1;

    //  Medžiagų koncentracijų masyvams išskiriama atmintis
    last_s    = new double[point_count];
    current_s = new double[point_count];
    last_p    = new double[point_count];
    current_p = new double[point_count];

    //  Priskiriamos pradinės ir kai kurios kraštinės sąlygos
    fill_array(last_s, point_count - 1, 0);
    fill_array(last_p, point_count - 1, 0);
    last_s[point_count - 1] = s0;
    last_p[point_count - 1] = p0;
    current_s[point_count - 1] = s0;
    current_p[point_count - 1] = p0;
    current_p[0] = 0;

    //  Pradiniai taskai
    //  concatenate_vals( last_s, s_list, point_count);
    //  concatenate_vals( last_p, p_list, point_count);
    t_list.push_back(0.);
    i_list.push_back(0.);
    //  Kiekvienam sluoksniui apskaičiuojami žingsniai pagal erdvę
    //  space_steps = new double[layer_count];//

    for (a = 0; a < layer_count; a++) {
        //  space_steps[a] = bio_info->layers[a].dx;
        dt = std::min(dt, pow(bio_info->layers[a].dx, 2)/  \
                      (std::max(bio_info->layers[a].Dp, bio_info->layers[a].Ds)*2));
    }
    //  Stabilumo sąlyga parenkamas "geras" dt
    if (dt != bio_info->dt) {
        if (callback_crunched != NULL) {
            std::stringstream m;
            m << "New delta t set: " << dt << " form old: " << bio_info->dt << "\n";
            //  m << print(bio_info);
            callback_crunched(ptr, 0, m.str());
        }
    }
    double vmaxt = bio_info->vmax*dt;


    std::clock_t start;
    start = std::clock();
    do {
        //  Iteruojama per biojutiklio sluoksnius, skaičiuojamos medžiagų koncentracijos
        layer = 0;
        //  Surenkami pirmojo sluoksnio parametrai
        enz_layer = bio_info->layers[layer].enz_layer;
        Ds = bio_info->layers[layer].Ds;
        Dp = bio_info->layers[layer].Dp;
        dx = bio_info->layers[layer].dx;

        for (a = 1; a < point_count - 1; a++) {
            //  Nustatome ar tai nėra sluoksnių sandūra
            is_boundary = !(a % n);

            //  Reikšmės sluoksnių sandūrose bus skaičiuojamos vėliau
            if (is_boundary) {
                //  Nustatome kuriame sluoksnyje esame
                layer++;
                //  Surenkami kito sluoksnio parametrai
                enz_layer = bio_info->layers[layer].enz_layer;
                Ds = bio_info->layers[layer].Ds;
                Dp = bio_info->layers[layer].Dp;
                dx = bio_info->layers[layer].dx;
            } else {
                //  Įskaičiuojama difuzijos įtaka
                current_s[a] = dt * Ds *  \
                               (last_s[a + 1] - 2 * last_s[a] + last_s[a - 1]) / (dx * dx) +  \
                               last_s[a];

                current_p[a] = dt * Dp *  \
                               (last_p[a + 1] - 2 * last_p[a] + last_p[a - 1]) / (dx * dx) +  \
                               last_p[a];

                //  Jeigu sluoksnis yra fermentinis, tuomet prisideda kinetika
                if (enz_layer) {
                    kinetics_part = vmaxt * last_s[a] /  \
                                    (km + last_s[a]);
                    current_s[a] -= kinetics_part;
                    current_p[a] += kinetics_part;
                }
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
        swap_arrays(&current_s, &last_s);
        swap_arrays(&current_p, &last_p);

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
            //  concatenate_vals( last_s, s_list, point_count);
            //  concatenate_vals( last_p, p_list, point_count);

            /*if (callback_crunched != NULL) {
                std::stringstream m;
                if(i > 1e-30) {
                    m << "Stationarity: " << (execution_time / i) * fabs(di / dt) << "\n";
                }

                callback_crunched(ptr, execution_time, m.str());
            }*/
        }

        //  Nustatoma ar tęsti simuliaciją
        switch (resp_t_meth) {
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
    } while (!response_time_reached);

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


    /*I = i_list;
    S = s_list;
    P = p_list;
    T = t_list;*/

    //  Atlaisvinama atmintis
    delete [] current_s;
    delete [] last_s;
    delete [] current_p;
    delete [] last_p;
}



void setGlobalParams(bio_params * p, unsigned int N, double dt,  \
                     resp_method type, double time , double ne,  \
                     const std::string file) {
    p->n = N;
    p->dt = dt;
    p->resp_t_meth = type;
    p->cond_t = time;
    p->ne = ne;
    p->out_file_name = file;
}

void setLocalParams(bio_params * p, double s0, double p0, double Km,  \
                    double Vmax, std::vector<std::vector<double> > inner) {
    p->s0 = s0;
    p->p0 = p0;
    p->km = Km;
    p->vmax = Vmax;
    p->layer_count = inner.size();
    p->layers = new layer_params[ inner.size() ];
    double dx(DBL_MAX), D(DBL_MIN), dt(p->dt);

    for (unsigned int i = 0; i < inner.size(); i++) {
        p->layers[i].Ds = inner[i][0];
        p->layers[i].Dp = inner[i][1];
        p->layers[i].dx = inner[i][2]/p->n;
        p->layers[i].enz_layer = inner[i][3];
        dx = std::min(dx, p->layers[i].dx);
        D = std::max(D, std::max(p->layers[i].Dp, p->layers[i].Ds));
    }
    dt = pow(dx, 2)/(2*D);
    if (dt < p->dt) {
        std::cout << "For converge conditions dt were changed from " \
                  << p->dt << " to " << dt << "!" << std::endl;
        p->dt = dt;
    }
}


}  //  namespace FDCalculator
