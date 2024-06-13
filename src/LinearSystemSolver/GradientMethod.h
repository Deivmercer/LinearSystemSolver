//
// Created by Davide Costantini on 19/05/24.
//

#ifndef GRADIENTMETHOD_H
#define GRADIENTMETHOD_H

#include "IterativeSolver.h"

class GradientMethod : public IterativeSolver
{
    void getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk) override;

    public:
        bool checkConvergence(const Matrix& A) override;
};

#endif //GRADIENTMETHOD_H
