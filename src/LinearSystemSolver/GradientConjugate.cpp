#include "GradientConjugate.h"

/*
 * 1: r^(k) = b − Ax^(k); rk -> residuo
 * 2: y^(k) = A*d^(k); yk
 * 3: z^(k) = A*r^(k); zk
 * 4: αk = (d^(k)· r^(k)) / (d^(k)· y^(k)); ak -> coefficiente, esso richiede un prodotto scalar
 * 5: x^(k+1) = x^(k) + αkd^(k);
 * 6: r^(k+1) = b − Ax^(k+1);
 * 7: w^(k) = Ar^(k+1);
 * 8: βk = (d^(k)· w^(k))/(d^(k)· y^(k));
 * 9: d^(k+1) = r^(k+1) − βkd^(k);
 */
void GradientConjugate::getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk)
{
    residual = b - (A * xk);
    direction = residual;

    Eigen::VectorXf yk = A * direction;
    Eigen::VectorXf zk = A * residual;
    float ak = (direction.dot(residual)) / (direction.dot(yk));
    x = x + (ak * direction);
    residual = b - (A * xk);
    Eigen::VectorXf wk = A * residual;
    float bk = (direction.dot(wk)) / (direction.dot(yk));
    direction = residual - (bk * direction);
}

bool GradientConjugate::checkConvergence(const Matrix& A)
{
    return Utils::isSymmetric(A) && Utils::isPositiveDefinite(A);
}
