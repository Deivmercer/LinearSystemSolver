//
// Created by Davide Costantini on 17/05/24.
//

#ifndef GAUSS_SEIDEL_GAUSSSEIDEL_H
#define GAUSS_SEIDEL_GAUSSSEIDEL_H

#include "types.h"
#include "IterativeSolver.h"

class GaussSeidel : public IterativeSolver
{
    Matrix getLowerMatrix(const Matrix& matrix);

    Matrix getUpperMatrix(const Matrix& matrix);

    void getNextXk(const Matrix &A, const Eigen::VectorXf &b, const Eigen::VectorXf &xk) override;

    public:
        bool checkConvergence(const Matrix& A) override;
};

#endif //GAUSS_SEIDEL_GAUSSSEIDEL_H
