#include <iostream>
#include "types.h"
#include "GaussSeidel.h"
#include "libs/unsupported/Eigen/SparseExtra"
#include "Utils.h"

int main()
{
    // Parsing the MatrixMarket file
    Matrix matrix;
    Eigen::loadMarket(matrix, "../sample-mtx/spa1.mtx");
    if(!Utils::checkSize(matrix))
    {
        throw std::invalid_argument("Matrix not valid");
    }
    // Initializing the initial guess
    Eigen::VectorXf x0;
    x0.resize(matrix.cols());
    x0.setOnes();

    // Initializing the b vector
    Eigen::VectorXf b;
    b.resize(matrix.cols());
    b.setOnes();

    GaussSeidel::Result r = GaussSeidel::solve(matrix, b, x0, 0.001f);

    return 0;
}
