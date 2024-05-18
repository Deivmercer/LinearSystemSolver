//
// Created by Davide Costantini on 18/05/24.
//

#include "Utils.h"

namespace Utils
{
    Matrix swapRows(const Matrix& matrix, int i, int j)
    {
        Eigen::SparseMatrix<float> permutation(matrix.rows(), matrix.cols());
        permutation.setIdentity();
        permutation.coeffRef(i, i) = 0;
        permutation.coeffRef(j, j) = 0;
        permutation.coeffRef(i, j) = 1;
        permutation.coeffRef(j, i) = 1;

        return permutation * matrix;
    }
}