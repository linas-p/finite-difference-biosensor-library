/*
 *  Copyright (c) Linas Petkevicius 2015
 *  Vilnius University
 *  GNU General Public license
 * */

#ifndef INCLUDE_FDCALCULATOR_UTILS_H_
#define INCLUDE_FDCALCULATOR_UTILS_H_
#define unlikely(x)    __builtin_expect((bool)(x), 0)

#include <stdio.h>
#include <vector>

namespace FDCalculator {
void swap_arrays(double **array1, double **array2);

void fill_array(double *array, int length, double value);

void concatenate_vals(double *, std::vector<double> *, int);
}  // namespace FDCalculator

#endif  // INCLUDE_FDCALCULATOR_UTILS_H_
