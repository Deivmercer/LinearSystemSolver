//
// Created by Davide Costantini on 19/05/24.
//

#include "GradientMethod.h"

#include "Utils.h"

void GradientMethod::getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk)
{
    residual = b - (A * xk);
    Eigen::VectorXf y = A * residual;
    float a = residual.transpose() * residual;
    float bk = residual.transpose() * y;
    float alpha = a / bk;
    x = xk + (alpha * residual);
}

bool GradientMethod::checkConvergence(const Matrix& A)
{
    return Utils::isSymmetric(A) && Utils::isPositiveDefinite(A);
}
