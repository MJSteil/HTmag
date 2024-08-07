//
// Mag_Star_O2 metric perturbations
// Created by M. J. Steil on 2017.02.27.
//

#include "../../include/Mag_Star_O2.hpp"

// Analytic exterior solutions

double Mag_Star_O2::m00_ext(double r, double deltaM0in) {
    double M = star->M;
    double M2 = M*M;
    double M3 = M2*M;

    double y = 1.-2.*M/r;
    double logy = log(y);

    double mu0mu0 = gsl_pow_2(mag_star->mu0);
    double Qe0Qe0 = gsl_pow_2(mag_star->Qe0);

    return deltaM0in + (mu0mu0*(-0.5625 - 0.375*logy + y*(0.75 + (-1.125 - 0.75*logy)*logy + y*(0.375 + 1.875*logy
                        + (-0.75 - 0.375*logy + 0.1875*y)*y))))/(M3*(-1. + y*(3. + (-3. + y)*y)))
                     + (Qe0Qe0*(-0.25 + 0.25*y))/M;
}

double Mag_Star_O2::n00_ext(double r, double deltaM0in) {
    double M = star->M;
    double M2 = M*M;
    double M4 = M2*M2;

    double y = 1.-2.*M/r;
    double logy = log(y);

    double mu0mu0 = gsl_pow_2(mag_star->mu0);
    double Qe0Qe0 = gsl_pow_2(mag_star->Qe0);

    return ((-0.5 + 0.5*y)*deltaM0in/(M*y))
                    + (mu0mu0*(-0.65625 - 0.1875*logy + y*(2.25 + (0.5625 + 0.375*logy)*logy + y*(-2.4375
                        + 0.1875*logy + (0.75 - 0.5625*logy + 0.09375*y)*y))))/(M4*y*(1. + (-2. + y)*y))
                    + (Qe0Qe0*(0.125 + (-0.25 + 0.125*y)*y))/(M2*y);
}

// Complete metric potentials

double Mag_Star_O2::m00(double r) {
    if(mag_star_O2_l0_status>=2){
        if(r<=star->R){
            return m00_of_r.f(r);
        }else {
            return m00_ext(r,deltaM0);
        }
    }else{
        GSL_ERROR_VAL ("m00: No interpolations supplied: run comp_l00_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::dm00(double r) {
    if(mag_star_O2_l0_status>=2){
        if(r<=star->R){
            return m00_of_r.df(r);
        }else {
            GSL_ERROR_VAL ("dm00: Not implemented for r>R.", GSL_FAILURE,GSL_FAILURE);
        }
    }else{
        GSL_ERROR_VAL ("dm00: No interpolations supplied: run comp_l00_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::n00(double r) {
    if(mag_star_O2_l0_status>=2){
        if(r<=star->R){
            return -h00_of_r.f(r)+2./3*mag_star->cjphi0*mag_star->aphi0(r)+ch00;
        }else {
            return n00_ext(r,deltaM0);
        }
    }else{
        GSL_ERROR_VAL ("dn00: No interpolations supplied: run comp_l00_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::dn00(double r) {
    if(mag_star_O2_l0_status>=2){
        if(r<=star->R){
            return -h00_of_r.df(r)+2./3*mag_star->cjphi0*mag_star->daphi0(r);
        }else {
            GSL_ERROR_VAL ("dn00: Not implemented for r>R.", GSL_FAILURE,GSL_FAILURE);
        }
    }else{
        GSL_ERROR_VAL ("dn00: No interpolations supplied: run comp_l00_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::h00(double r) {
    if(mag_star_O2_l0_status>=2){
        if(r<=star->R){
            return h00_of_r.f(r);
        }else {
            return 0;
        }
    }else{
        GSL_ERROR_VAL ("h00: No interpolations supplied: run comp_l00_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::dh00(double r) {
    if(mag_star_O2_l0_status>=2){
        if(r<=star->R){
            return h00_of_r.df(r);
        }else {
            return 0;
        }
    }else{
        GSL_ERROR_VAL ("h00: No interpolations supplied: run comp_l00_int().", GSL_FAILURE,GSL_FAILURE);
    }
}
