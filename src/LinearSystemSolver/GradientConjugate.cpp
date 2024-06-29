#include "GradientConjugate.h"

#include <iostream>
/*
 * 1: r^(k) = b − Ax^(k);  parto da r0 definito nel .h, da dopo il corrente sarà il successivo
 * 2: y^(k) = A*d^(k); yk
 * 3: z^(k) = A*r^(k); zk
 * 4: αk = (d^(k)· r^(k)) / (d^(k)· y^(k)); ak -> coefficiente, esso richiede un prodotto scalar
 * 5: x^(k+1) = x^(k) + αk * d^(k);
 * 6: r^(k+1) = b − Ax^(k+1);
 * 7: w^(k) = Ar^(k+1);
 * 8: βk = (d^(k)· w^(k))/(d^(k)· y^(k));
 * 9: d^(k+1) = r^(k+1) − βkd^(k);
 */
void GradientConjugate::getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk)
{
    Eigen::VectorXf yk = A * direction;
    Eigen::VectorXf zk = A * residual;
    float ak = (direction.dot(residual)) / (direction.dot(yk));
    x = x + (ak * direction); //xk+1
    residual = b - (A * x); //rk+1
    Eigen::VectorXf wk = A * residual;
    float bk = (direction.dot(wk)) / (direction.dot(yk));
    direction = residual - (bk * direction); //dk+1

}

bool GradientConjugate::checkConvergence(const Matrix& A)
{
    return Utils::isSymmetric(A) && Utils::isPositiveDefinite(A);
}