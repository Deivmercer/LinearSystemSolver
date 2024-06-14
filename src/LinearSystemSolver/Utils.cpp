//
// Created by Davide Costantini on 18/05/24.
//

#include "Utils.h"

#include <iostream>
#include <Eigen/Core>
#include "CholeskyDecomposition.h"

Matrix Utils::swapRows(const Matrix &matrix, int i, int j) {
    Eigen::SparseMatrix<float> permutation(matrix.rows(), matrix.cols());
    permutation.setIdentity();
    permutation.coeffRef(i, i) = 0;
    permutation.coeffRef(j, j) = 0;
    permutation.coeffRef(i, j) = 1;
    permutation.coeffRef(j, i) = 1;

    return permutation * matrix;
}

bool Utils::checkSize(Matrix &matrix) {
    int rowSize = matrix.rows();
    int colSize = matrix.cols();

    return (rowSize == colSize ? true : false);
}

bool Utils::checkDiagonalZero(Matrix &matrix) {
    // check diagonal
    for (int index = 0; index < matrix.rows(); index++) {
        if (matrix.coeff(index, index) == 0) {
            std::stringstream stream;
            stream << "Found zero at position " << index << std::endl;
            throw std::invalid_argument(stream.str());
        }
    }
    return false;
}

Eigen::VectorXf Utils::invertDiagonal(const Matrix &matrix) {
    Eigen::VectorXf diag;
    diag.resize(matrix.rows());

    // check diagonal
    for (int index = 0; index < matrix.rows(); index++) {
        if (matrix.coeff(index, index) == 0) {
            std::stringstream stream;
            stream << "Found zero at position " << index << std::endl;
            throw std::invalid_argument(stream.str());
        }
        diag[index] = 1 / matrix.coeff(index, index);
    }

    return diag;
}

bool Utils::thresholdReached(const Matrix &A, const Eigen::VectorXf &b, const Eigen::VectorXf &xk,
                             const float tolerance) {
    // ||b - Axk|| / ||b|| < tol
    Eigen::VectorXf axk = A * xk;

    return (euclideanNorm(b, axk) / euclideanNorm(b)) < tolerance;
}

float Utils::euclideanNorm(const Eigen::VectorXf &x) {
    float sum = 0;
    for (int i = 0; i < x.size(); ++i)
        sum += x.coeff(i) * x.coeff(i);

    return sqrt(sum);
}

float Utils::euclideanNorm(const Eigen::VectorXf &x, const Eigen::VectorXf &y) {
    assert(x.size() == y.size());

    float sum = 0;
    for (int i = 0; i < x.size(); ++i)
        sum += (x.coeff(i) - y.coeff(i)) * (x.coeff(i) - y.coeff(i));

    return sqrt(sum);
}

bool Utils::isSymmetric(const Matrix &A) {
    for (int i = 0; i < A.rows(); ++i)
        for (int j = 0; j < A.cols(); ++j)
            if (A.coeff(i, j) != A.coeff(j, i))
                return false;

    return true;
}

bool Utils::isPositiveDefinite(const Matrix &A) {
    try {
        CholeskyDecomposition::solve(A);
    } catch (std::invalid_argument &e) {
        return false;
    }
    return true;
}


bool Utils::isStrictlyDiagonallyDominant(const Matrix &A) {
    Eigen::VectorXf diag;
    diag.resize(A.rows());
    float sumDiag;

    for (int i = 0; i < diag.size(); i++) {
        sumDiag = 0;
        //SOMMO VALORI DIAGONALE
        for (int j = 0; j < diag.size(); ++j) {
            if (j != i) {
                sumDiag += std::abs(A.coeff(i, j));
            }
        }
        //CONTROLLO REGOLA
        if (std::abs(A.coeff(i, i)) <= sumDiag) {
            std::cout << "NON CONVERGE ALL'ITERATA: " << i << std::endl;
            std::cout << "A.coeff: " << A.coeff(i, i);
            std::cout << " <= sumDiag: " << sumDiag << std::endl;
            return false;
        }
    }
    return true;
}
