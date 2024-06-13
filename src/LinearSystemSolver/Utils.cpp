//
// Created by Davide Costantini on 18/05/24.
//

#include "Utils.h"

#include <iostream>
#include <Eigen/Core>

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
                std::stringstream stream;
                stream << "Found zero at position " << index << std::endl;
                throw std::invalid_argument(stream.str());
            }
        }
        return false;
    }

    Eigen::VectorXf invertDiagonal(const Matrix& matrix)
    {
        Eigen::VectorXf diag;
        diag.resize(matrix.rows());

        // check diagonal
        for (int index = 0; index < matrix.rows(); index++)
        {
            if (matrix.coeff(index, index) == 0)
            {
                std::stringstream stream;
                stream << "Found zero at position " << index << std::endl;
                throw std::invalid_argument(stream.str());
            }
            diag[index] = 1 / matrix.coeff(index, index);
        }

        return diag;
    }

    bool thresholdReached(const Matrix& A, const Eigen::VectorXf& b, const Eigen::VectorXf& xk, const float tolerance)
    {
        // ||b - Axk|| / ||b|| < tol
        Eigen::VectorXf axk = A * xk;

        return (euclideanNorm(b, axk) / euclideanNorm(b)) < tolerance;
    }

    float euclideanNorm(const Eigen::VectorXf& x) {

        float sum = 0;
        for (int i = 0; i < x.size(); ++i)
            sum += x.coeff(i) * x.coeff(i);

        return sqrt(sum);
    }

    float euclideanNorm(const Eigen::VectorXf& x, const Eigen::VectorXf& y) {

        assert(x.size() == y.size());

        float sum = 0;
        for (int i = 0; i < x.size(); ++i)
            sum += (x.coeff(i) - y.coeff(i)) * (x.coeff(i) - y.coeff(i));

        return sqrt(sum);
    }
}
