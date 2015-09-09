#include "calculator.h"

Calculator::Calculator(){

    params_ = new bio_params;

}

Calculator::~Calculator() {
    delete [] params_->layers;
    delete params_;
}

void Calculator::setGlobalParams(unsigned int N, double dt, resp_method type, double time , double ne, const std::string file) {
    params_->n = N;
    params_->dt = dt;
    params_->resp_t_meth = type;
    params_->cond_t = time;
    params_->ne = ne;
    params_->out_file_name = file;

}
void Calculator::saveResult() {
    output_file_.open(params_->out_file_name.c_str(), std::ofstream::out | std::ofstream::app);
    for(unsigned int i = 0; i < i_rez_.size(); i++) {
        output_file_ << t_rez_[i] << "," << i_rez_[i] << "\n";
    }
    output_file_.close();
}
void Calculator::setLocalParams(double s0, double p0, double Km, double Vmax, std::vector<std::vector<double> > inner) {
    params_->s0 = s0;
    params_->p0 = p0;
    params_->km = Km;
    params_->vmax = Vmax;
    params_->layer_count = inner.size();
    params_->layers = new layer_params[ inner.size() ];
    double dx(DBL_MAX), D(DBL_MIN), dt(params_->dt);

    for(unsigned long i = 0; i < inner.size(); i++) {
        params_->layers[i].Ds = inner[i][0];
        params_->layers[i].Dp = inner[i][1];
        params_->layers[i].dx = inner[i][2]/params_->n;
        params_->layers[i].enz_layer = inner[i][3];
        dx = std::min(dx, params_->layers[i].dx);
        D = std::max(D, std::max(params_->layers[i].Dp, params_->layers[i].Ds));
    }
    dt = pow(dx,2)/(2*D);
    if( dt < params_->dt ) {
        std::cout << "For converge conditions dt were changed from " << params_->dt << " to " << dt << "!" << std::endl;
        params_->dt = dt;

    }
}

std::string Calculator::print() {
    std::stringstream m;
    m << "Global conditions:\n  n = " <<  params_->n << " dt = " << params_->dt << " type " << params_->resp_t_meth << " ne = " << params_->ne << " output goes to: " << params_->out_file_name << "\n";

    m << "Initial conditions: \n  s0= " << params_->s0 << " p0= " << params_->p0 << " km= " << params_->km << " vmax= " << params_->vmax << " number of layers= " << params_->layer_count << " \n";

    for( unsigned int i = 0; i < params_->layer_count; i++) {
        m << "  layer " << i << " params: Ds= " << params_->layers[i].Ds << " Dp= " << params_->layers[i].Dp << " d= " << params_->layers[i].dx*params_->n << " with accuracy: " << params_->layers[i].dx << " is enzyme? " << params_->layers[i].enz_layer << "\n";
    }
    return m.str();
}

void Calculator::solve(struct bio_params *params) {

    unsigned int n(params_->n);
    unsigned int layer_count(params_->layer_count);
    unsigned int grid_points(layer_count * n + 1);

    double *last_s    = new double[grid_points];//$
    double *current_s = new double[grid_points];//$
    double *last_p    = new double[grid_points];//$
    double *current_p = new double[grid_points];//$
	double i(0.), last_i(0.);

    fill_array( last_s, grid_points, 0 );
    fill_array( last_p, grid_points, 0 );
    current_s[grid_points-1] = last_s[grid_points-1] = params_->s0;
    current_p[grid_points-1] = last_p[grid_points-1] = params_->p0;
    current_s[0] = last_s[0] = 0;
    current_p[0] = last_p[0] = 0;

    //Konstantos naudojamos cikle
    const double vmax(params_->vmax);
    const double km(params_->km);
    const double i_const(params_->ne * F * params_->layers[0].Dp / params_->layers[0].dx);
    double kinetics_part(0.), s0(params_->s0);

    layer_params current_layer, next_layer;
    double dx, dx0;
    double enzyme, layer_ws, layer_wp;
    unsigned int a, layer;

    //Iteracija pagal laiką
    const double dt(params_->dt);
    double execution_time(0.);
    unsigned int t(0);
    bool response_time_reached(false);

	i_rez_.push_back(i);
    t_rez_.push_back(execution_time);
   	concatenate_vals( last_s, s_rez_, grid_points);
   	concatenate_vals( last_p, p_rez_, grid_points);


    std::clock_t start;
    start = std::clock();
    do {

        for (layer = 0; layer < layer_count; layer++) {

            current_layer = params_->layers[layer];
            dx = pow(current_layer.dx, 2);
            enzyme = ((double)current_layer.enz_layer) * dt * vmax;
            layer_ws = current_layer.Ds * dt / dx;
            layer_wp = current_layer.Dp * dt / dx;

            for (a = n * layer + 1; a < n * (layer + 1); a++) {

                current_s[a] = last_s[a] + layer_ws * (last_s[a + 1] - 2 * last_s[a] + last_s[a - 1]);
                kinetics_part = enzyme*last_s[a]/(km+last_s[a]);
                current_s[a] -= kinetics_part;

                current_p[a] = last_p[a] + layer_wp * (last_p[a + 1] - 2 * last_p[a] + last_p[a - 1]);
                current_p[a] += kinetics_part;
            }
		}


        for (layer = 0; layer < layer_count - 1; layer++) {
            //Apskaičiuojame kuriame taške yra layer ir layer + 1 sluoksnių sandūra
            a = n * (layer + 1);

            current_layer = params_->layers[layer];
            next_layer = params_->layers[layer + 1];

            dx = next_layer.Ds * current_layer.dx;
            dx0 = current_layer.Ds * next_layer.dx;
            current_s[a] = ((dx * current_s[a + 1])+ (dx0 * current_s[a - 1]))/((dx + dx0));

            dx = next_layer.Dp * current_layer.dx;
            dx0 = current_layer.Dp * next_layer.dx;
            current_p[a] = ((dx * current_p[a + 1] ) + (dx0 * current_p[a - 1]))/((dx + dx0));
        }

        //Kraštinė substrato nepratekėjimo sąlyga
        current_s[0] = current_s[1];
        current_s[grid_points-1] = s0;
        current_p[0] = 0.;
        current_p[grid_points-1] = 0.;

        //Apskaičiuojamas laikas
        execution_time += dt;
        t++;

        //Skaičiuojamas srovės tankis
        i = i_const * (current_p[1] - current_p[0]);
        response_time_reached = (i > 1e-30) && ((execution_time / i) * (fabs(i - last_i) / dt) <= EPSILON);
        last_i = i;

        //Masyvai sukeičiami vietomis
        swap_arrays(&current_s, &last_s);
        swap_arrays(&current_p, &last_p);


        if ((t % (INTERVAL)) == 0) {
            i_rez_.push_back(i);
            t_rez_.push_back(execution_time);
        	concatenate_vals( last_s, s_rez_, grid_points);
        	concatenate_vals( last_p, p_rez_, grid_points);
        }

    } while(!response_time_reached );

    std::cout << "t " << t << " " << execution_time << std::endl;
    std::cout << "time " << (std::clock()-start)/ (double) CLOCKS_PER_SEC << std::endl;

	i_rez_.push_back(i);
    t_rez_.push_back(execution_time);
   	concatenate_vals( last_s, s_rez_, grid_points);
   	concatenate_vals( last_p, p_rez_, grid_points);

    delete [] last_s;
    delete [] last_p;
    delete [] current_p;
    delete [] current_s;
}

