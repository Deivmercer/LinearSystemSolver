//
// Created by Davide Costantini on 18/05/24.
//

#include "Utils.h"

#include <iostream>

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

    bool checkSize(Matrix& matrix)
    {
        int rowSize = matrix.rows();
        int colSize = matrix.cols();

        return (rowSize == colSize ? true : false);
    }

    bool checkDiagonalZero(Matrix& matrix)
    {
        int diag = matrix.rows();

        // check diagonal
        for (int index = 0; index < diag; index++)
        {
            if (matrix.coeff(index, index) == 0)
            {
                throw std::invalid_argument("Found zero at:" << std::endl << "row = " << index << std::endl << " - column = " << index);
            }
        }
        return false;
    }

    Eigen::VectorXf invertDiagonal(Matrix& matrix)
    {
        Eigen::VectorXf diag;

        // check diagonal
        for (int index = 0; index < matrix.rows(); index++)
        {
            if (matrix.coeff(index, index) == 0)
            {
                throw std::invalid_argument("Found zero at:" << std::endl << "row = " << index << std::endl << " - column = " << index);
            }
            diag[index] = 1 / matrix.coeff(index, index);
        }
        return diag;
    }
}
