//
// Created by Davide Costantini on 13/06/24.
//

#include "CholeskyDecomposition.h"

Matrix CholeskyDecomposition::solve(const Matrix& A)
{
    Matrix Astar = A;
    Matrix R(A.rows(), A.cols());

    for (int k = 0; k < A.rows(); ++k)
    {
        R.coeffRef(k, k) = sqrt(Astar.coeff(k, k));
        if (R.coeff(k ,k) <= 0)
            throw std::invalid_argument("Could not compute Cholesky Decomposition with the given matrix.");
        for (int j = k + 1; j < A.rows(); ++j)
            R.coeffRef(k, j) = Astar.coeff(k, j) / R.coeff(k, k);
        float rho = 1 / Astar.coeff(k, k);
        for (int j = k + 1; j < A.rows(); ++j)
            for (int i = k + 1; i < A.rows(); ++i)
                Astar.coeffRef(i, j) = Astar.coeff(i, j) - rho * Astar.coeff(i, j) * Astar.coeff(k, j);
    }

    return R;
}
