#ifndef EXPLICITCALCULATOR_H
#define EXPLICITCALCULATOR_H

#include "biosensor_information.h"
#include <string>
#include <vector>

void calculate_explicitly(struct bio_params *bio_info, void *ptr,
                          void (*callback_crunched)(void *, int, std::string),
                          std::vector<double> & ,
                          std::vector<double> & ,
                          std::vector<double> & ,
                          std::vector<double> &
                         )

#endif
