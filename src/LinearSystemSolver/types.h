//
// Created by Davide Costantini on 18/05/24.
//

#ifndef METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_TYPES_H
#define METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_TYPES_H

#include "../libs/Eigen/SparseCore"

typedef Eigen::SparseMatrix<float, Eigen::RowMajor> Matrix;

typedef Eigen::Block<Matrix, 1, Eigen::Dynamic, Eigen::RowMajor> MatrixRow;

#endif //METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_TYPES_H
