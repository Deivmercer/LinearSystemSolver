#include "Jacobi.h"

void Jacobi::getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk)
{
    //x^(k+1) = x^k + P^(-1) * (Ax^(k) - b)
    //residual = b − Ax^(k)
    residual = b - (A * x);
    Eigen::VectorXf rk = invP.asDiagonal() * residual;
    x = x + rk;
}

