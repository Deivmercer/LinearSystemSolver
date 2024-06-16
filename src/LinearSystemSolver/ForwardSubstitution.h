//
// Created by Davide Costantini on 18/05/24.
//

#ifndef METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_FORWARDSUBSTITUTION_H
#define METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_FORWARDSUBSTITUTION_H

#include "types.h"

namespace ForwardSubstitution
{
    Eigen::VectorXf solve(const Matrix& L, const Eigen::VectorXf& b);
}

#endif //METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_FORWARDSUBSTITUTION_H
