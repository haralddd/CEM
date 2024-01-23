//
// Created by haralg on 15.12.23.
//
#include <cstdlib>
#include "args-parser/all.hpp"

int main(int argc, char **argv) {

    parse_args(argc, argv);

    using namespace std;
    double wavelength = 632.8e-9; // wavelength of He-Ne laser
    double surf_length = 50 * wavelength;
    double surf_rms_height = wavelength / 20;
    double surf_cor_length = wavelength / 4;


    exit(EXIT_SUCCESS);
}