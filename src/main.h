//
// Created by Davide Costantini on 01/06/24.
//

#ifndef METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_MAIN_H
#define METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_MAIN_H

#include <string>

const std::string PATH = "../sample-mtx/";

struct DemoMatrix
{
    std::string name;
    float tolerance;

    std::string path() const
    {
        return PATH + name;
    }
};

const DemoMatrix DEMO_MATRICES[] = {
        {"spa1.mtx", 0.001},
        {"spa2.mtx", 0.00001},
        {"vem1.mtx", 0.0000001},
        {"vem2.mtx", 0.000000001}
};

#endif //METODI_DEL_CALCOLO_SCIENTIFICO_1_BIS_MAIN_H
