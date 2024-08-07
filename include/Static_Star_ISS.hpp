//
// Created by M. J. Steil on 2017.02.08.
//

#ifndef HTMAG_STATIC_ISS_STAR_HPP
#define HTMAG_STATIC_ISS_STAR_HPP


#include "../../eos/eos_irf.hpp"
#include <gsl/gsl_math.h>

class Static_Star_ISS{
public:

    Static_Star_ISS();
    Static_Star_ISS(double Min /*[Ms]*/, double Rin /*[km]*/);
    Static_Star_ISS(double rhoCin /*[km**-2]*/);

    void setMR(double Min /*[Ms]*/, double Rin /*[km]*/);
    void setZ(double Zin /*[1]*/);
    void sethCrhoC(double hC /*[1]*/, double rhoC /*[km**-2]*/);

    void setAux();

    eos_irf *star_eos;

    double rhoC, PC, hC, M, MB, R, Z,twoMR,twoMR3;

    // Thermodynamic properties
    double h(double r);
    double dh(double r);

    double P(double r);
    double dP(double r);

    double rho(double r);
    double drhodh(double r);

    double nbar(double r);
    double dnbardh(double r);

    // Metric potentials
    double expLambda(double r);
    double dLambda(double r);
    double expNu(double r);
    double dNu(double r);
    double m(double r);
    double dm(double r);

    // Limits
    void show_limits();
    void show_discontinuity();
};


#endif //NEUTRONSTAR_STATIC_ISS_STAR_HPP
