//
// Mag_Star_O2 metric perturbations
// Created by M. J. Steil on 2017.02.27.
//

#include "../../include/Mag_Star_O2.hpp"

// Analytic exterior solutions

double Mag_Star_O2::omegaQ00_ext(double r, double deltaJ20_00in) {
    double M = star->M;
    double M2 = M*M;
    double M4 = M2*M2;

    double y = 1.-2.*M/r;
    double y2 = y*y;

    double mu0 = mag_star->mu0;

    return 2.*deltaJ20_00in/gsl_pow_3(r) + (mu0*(-0.25 + y*(-0.375 + (0.75 - 0.125*y)*y - 0.75*log(y))))/M4;
}

double Mag_Star_O2::domegaQ00_ext(double r, double deltaJ20_00in) {
    double M = star->M;
    double M2 = M*M;
    double M4 = M2*M2;

    double y = 1.-2.*M/r;
    double y2 = y*y;

    double mu0 = mag_star->mu0;

    return -6.*deltaJ20_00in/gsl_pow_4(r) + (mu0*(-0.5625 + y*(1.875 + y*(-2.25 + (1.125 - 0.1875*y)*y - 0.375*log(y)) + 0.75*log(y)) - 0.375*log(y)))/M4/M;
}

// Complete metric potentials

double Mag_Star_O2::omegaQ00(double r) {
    if(mag_star_O2_l1_status>=2){
        if(r<=star->R&&r!=0.){
            return omegaQ00_of_r.f(r);
        }else if(r==0){
            return omegaQ00_i[0];
        } else
            return omegaQ00_ext(r,deltaJ20_00);
    } else{
        GSL_ERROR_VAL ("omegaQ00: No interpolations supplied: run comp_l100_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::domegaQ00(double r) {
    if(mag_star_O2_l1_status>=2){
        if(r<=star->R){
            return domegaQ00_of_r.f(r);
        }else {
            return domegaQ00_ext(r,deltaJ20_00);
        }
    }else{
        GSL_ERROR_VAL ("domegaQ00: No interpolations supplied: run comp_l100_int().", GSL_FAILURE,GSL_FAILURE);
    }
}