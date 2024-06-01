//
// Created by Davide Costantini on 18/05/24.
//

#ifndef METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_UTILS_H
#define METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_UTILS_H

#include "types.h"

namespace Utils
{
    Matrix swapRows(const Matrix& matrix, int i, int j);

    bool checkSize(Matrix& matrix);

    bool checkDiagonalZero(Matrix& matrix);

    Eigen::VectorXf invertDiagonal(const Matrix& matrix);

    bool thresholdReached(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, const float tolerance);

    float euclideanNorm(const Eigen::VectorXf& x);

    float euclideanNorm(const Eigen::VectorXf& x, const Eigen::VectorXf& y);
}

#endif //METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_UTILS_H
