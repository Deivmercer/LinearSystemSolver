//
// Created by Davide Costantini on 17/05/24.
//

#include "GaussSeidel.h"
#include "BackwardSubstitution.h"

Matrix GaussSeidel::getLowerMatrix(const Matrix& matrix)
{
    return matrix.triangularView<Eigen::Lower>();
}

Matrix GaussSeidel::getUpperMatrix(const Matrix& matrix)
{
    return matrix.triangularView<Eigen::Upper | Eigen::ZeroDiag>() * -1;
}

void GaussSeidel::getNextXk(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk)
{
    // TODO: controllare teorema 4.3
    // Extracting the upper and lower triangular components

    Matrix P = GaussSeidel::getLowerMatrix(A);

    residual = b - (A * xk);
    Eigen::VectorXf y = BackwardSubstitution::solve(P, residual);
    x = xk + y;
}

bool GaussSeidel::checkConvergence(const Matrix& A)
{
    // TODO: stub
    return false;
}
