//
// Mag_Star_O2 metric perturbations
// Created by M. J. Steil on 2017.02.27.
//

#include "../../include/Mag_Star_O2.hpp"

// Derived quantites

double Mag_Star_O2::xi00(double r) {
    if(mag_star_O2_l0_status>=2){
        if(r>0){
            return 2./(star->dNu(r))*h00(r);
        }else{
            return 0;
        }
    }else{
        GSL_ERROR_VAL ("xi00: No interpolations supplied: run comp_l00_int().", GSL_FAILURE,GSL_FAILURE);
    }
}


double Mag_Star_O2::xi20(double r){
    if(mag_star_O2_l2_status>=2){
        if(r>0){
            return 2./(star->dNu(r))*h20(r);
        }else{
            return 0;
        }
    }else{
        GSL_ERROR_VAL ("xi20: No interpolations supplied: run comp_l20_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::rbar(double r, double theta){
    if (r > 0) {
        return r + BB()*( xi00(r) + xi20(r) * (1. - 1.5 * gsl_pow_2(sin(theta))) );
    } else {
        return 0;
    }
}

double Mag_Star_O2::h(double r, double theta){
    if (r >= 0) {
        return star->h(r) + BB()*( h00(r) + h20(r) * (1. - 1.5 * gsl_pow_2(sin(theta))) );
    } else {
        return 0;
    }
}