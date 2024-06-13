cmake_minimum_required(VERSION 3.28)
project(Metodi_del_Calcolo_Scientifico_1_bis)

set(CMAKE_CXX_STANDARD 17)

add_executable(Metodi_del_Calcolo_Scientifico_1_bis src/main.cpp
        src/main.h
        src/LinearSystemSolver/GaussSeidel.cpp
        src/LinearSystemSolver/GaussSeidel.h
        src/LinearSystemSolver/types.h
        src/LinearSystemSolver/Utils.cpp
        src/LinearSystemSolver/Utils.h
        src/LinearSystemSolver/BackwardSubstitution.cpp
        src/LinearSystemSolver/BackwardSubstitution.h
        src/LinearSystemSolver/Jacobi.h
        src/LinearSystemSolver/Jacobi.cpp
        src/LinearSystemSolver/GradientMethod.cpp
        src/LinearSystemSolver/GradientMethod.h
        src/LinearSystemSolver/GradientConjugate.h
        src/LinearSystemSolver/GradientConjugate.cpp
        src/LinearSystemSolver/IterativeSolver.cpp
        src/LinearSystemSolver/IterativeSolver.h)

target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC
        src/LinearSystemSolver/libs
)