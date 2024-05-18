//
// Created by Davide Costantini on 18/05/24.
//

#ifndef METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_BACKWARDSUBSTITUTION_H
#define METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_BACKWARDSUBSTITUTION_H

#include "types.h"

namespace BackwardSubstitution
{
    Eigen::VectorXf solve(const Matrix& A, const Eigen::VectorXf& b);
}

#endif //METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_BACKWARDSUBSTITUTION_H
