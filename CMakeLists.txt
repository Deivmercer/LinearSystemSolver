cmake_minimum_required(VERSION 3.28)
project(Metodi_del_Calcolo_Scientifico_1_bis)

set(CMAKE_CXX_STANDARD 17)

add_executable(Metodi_del_Calcolo_Scientifico_1_bis src/main.cpp
        src/GaussSeidel.cpp
        src/GaussSeidel.h
        src/types.h
        src/Utils.cpp
        src/Utils.h
        src/BackwardSubstitution.cpp
        src/BackwardSubstitution.h
        src/Jacobi.h
        src/Jacobi.cpp
        src/GradientMethod.cpp
        src/GradientMethod.h)
