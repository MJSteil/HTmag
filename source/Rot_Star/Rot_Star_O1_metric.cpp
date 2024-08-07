//
// Metric potentials
// Created by M. J. Steil on 2017.02.27.
//

#include "../../include/Rot_Star_O1.hpp"

double Rot_Star_O1::omegabar0_ext(double r, int type){
    double f;
    switch ( type )
    {
        case 0: // omegabar0(r)
            f = 1-2*J0/gsl_pow_3(r);
            break;
        case 1: // domegabar0(r)
            f = 6*J0/gsl_pow_4(r);
            break;
        case 2: // ddomegabar(r)
            f = -24*J0/gsl_pow_5(r);
            break;
        default:
            GSL_ERROR_VAL ("omegabar0_ext: No valid type: 0:omegabar0(r), 1:domegabar0(r), 2: ddomegabar(r).", GSL_FAILURE,GSL_FAILURE);
    }
    return f;
};

double Rot_Star_O1::omegabar0(double r) {
    if(rot_star_status<2) {
        GSL_ERROR_VAL ("omegabar0(r): No interpolations supplied: run comp_omegabar0_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r==0.){
        return omegabar0C;
    }else if(r<=star->R){
        return omegabar0_of_r.f(r);
    } else {
        return omegabar0_ext(r, 0);
    };
}

double Rot_Star_O1::domegabar0(double r) {
    if(rot_star_status<2) {
        GSL_ERROR_VAL ("domegabar0(r): No interpolations supplied: run comp_omegabar0_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<star->R){
        return domegabar0_of_r.f(r);
    } else {
        return omegabar0_ext(r, 1);
    };
}

double Rot_Star_O1::omega0(double r) {
    if(rot_star_status<2) {
        GSL_ERROR_VAL ("omega0(r): No interpolations supplied: run comp_omegabar0_int().",GSL_FAILURE,GSL_FAILURE);
    }
    return 1 - omegabar0(r);

}

double Rot_Star_O1::domega0(double r) {
    if(rot_star_status<2) {
        GSL_ERROR_VAL ("domega0(r): No interpolations supplied: run comp_omegabar0_int().",GSL_FAILURE,GSL_FAILURE);
    }
    return -domegabar0(r);
}


