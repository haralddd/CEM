#ifndef _ARGUMENT_UTILS_H_
#define _ARGUMENT_UTILS_H_

enum SurfType { FLAT, BUMP, GAUSSIAN };
enum Polarization { P, S };


typedef struct options_struct {

    SurfType surf_type;
    Polarization polarization;
    int Nq;
    double permittivity;
    double permeability;
    double wavelength;
} OPTIONS;


OPTIONS *parse_args ( int argc, char **argv );

void help ( char const *exec, char const opt, char const *optarg );

#endif
