/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 * */
#include <FDCalculator/utils.h>
#include <vector>

namespace FDCalculator {
void swap_arrays(double **array1, double **array2) {
    double *temp;
    temp = *array1;
    *array1 = *array2;
    *array2 = temp;
}

void fill_array(double *array, int length, double value) {
    int a;

    for (a = 0; a < length; a++)
        array[a] = value;
}

void concatenate_vals(double *x, std::vector<double> * out, int length) {
    std::vector<double> tmp(x, (x + length));
    out->insert(out->end(), tmp.begin(), tmp.end());
}
}  // namespace FDCalculator
