#include <iostream>
#include "main.h"
#include "LinearSystemSolver/GaussSeidel.h"
#include "LinearSystemSolver/GradientConjugate.h"
#include "LinearSystemSolver/GradientMethod.h"
#include "LinearSystemSolver/Jacobi.h"
#include "LinearSystemSolver/libs/unsupported/Eigen/SparseExtra"
#include "LinearSystemSolver/Utils.h"

int main()
{
    DemoMatrix demo = DEMO_MATRICES[0];

    // Parsing the MatrixMarket file
    Matrix matrix;
    Eigen::loadMarket(matrix, demo.path());
    if(!Utils::checkSize(matrix))
    {
        throw std::invalid_argument("Matrix not valid");
    }

    std::cout << "Matrix: " << demo.name << std::endl;
    std::cout << "Tolerance: " << demo.tolerance << std::endl << std::endl;

    // Initializing the initial guess
    Eigen::VectorXf x0;
    x0.resize(matrix.cols());
    x0.setZero();

    // Initializing the b vector
    Eigen::VectorXf x;
    x.resize(matrix.cols());
    x.setOnes();
    Eigen::VectorXf b;
    b = matrix * x;


    GaussSeidel gs;
    gs.solve(matrix, b, x0, demo.tolerance);
    std::cout << "[GaussSeidel::solve] Iterations: " << gs.getIterations() << std::endl;
    std::cout << "[GaussSeidel::solve] x: " << gs.getX().transpose() << std::endl;
    std::cout << "[GaussSeidel::solve] residual: " << gs.getResidual().transpose() << std::endl << std::endl;
    float gserror = Utils::euclideanNorm(gs.getX(), x) / Utils::euclideanNorm(x);
    std::cout << "[GaussSeidel::solve] Relative error: " << gserror << std::endl << std::endl;

    GradientMethod gm;
    gm.solve(matrix, b, x0, demo.tolerance);
    std::cout << "[GaussSeidel::solve] Iterations: " << gm.getIterations() << std::endl;
    std::cout << "[GradientMethod::solve] x: " << gm.getX().transpose() << std::endl;
    std::cout << "[GradientMethod::solve] residual: " << gm.getResidual().transpose() << std::endl << std::endl;
    float gmerror = Utils::euclideanNorm(gm.getX(), x) / Utils::euclideanNorm(x);
    std::cout << "[GradientMethod::solve] Relative error: " << gmerror << std::endl << std::endl;

    Jacobi jm;
    jm.solve(matrix, b, x0, demo.tolerance);
    std::cout << "[Jacobi::solve] Iterations: " << jm.getIterations() << std::endl;
    std::cout << "[Jacobi::solve] x: " << jm.getX().transpose() << std::endl;
    std::cout << "[Jacobi::solve] residual: " << jm.getResidual().transpose() << std::endl << std::endl;
    float jmerror = Utils::euclideanNorm(jm.getX(), x) / Utils::euclideanNorm(x);
    std::cout << "[Jacobi::solve] Relative error: " << jmerror << std::endl << std::endl;

    GradientConjugate gc;
    gc.solve(matrix, b, x0, demo.tolerance);
    std::cout << "[GradientConjugate::solve] Iterations: " << gc.getIterations() << std::endl;
    std::cout << "[GradientConjugate::solve] x: " << gc.getX().transpose() << std::endl;
    std::cout << "[GradientConjugate::solve] direction: " << gc.getDirection().transpose() << std::endl;
    std::cout << "[GradientConjugate::solve] residual: " << gc.getResidual().transpose() << std::endl << std::endl;
    float gcerror = Utils::euclideanNorm(gc.getX(), x) / Utils::euclideanNorm(x);
    std::cout << "[GradientConjugate::solve] Relative error: " << gcerror << std::endl << std::endl;

    return 0;
}
