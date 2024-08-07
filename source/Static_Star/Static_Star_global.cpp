//
// Global quantities
// Created by M. J. Steil on 2017.02.24.
//

#include "../../include/Static_Star.hpp"

void Static_Star::comp_B(int print) {
    if (static_star_status<3) {
        GSL_ERROR_VAL ("comp_B: Background star with interpolations required.", GSL_FAILURE,);
    }

    gsl_function_pp_fkt comp_B_dBdr = gsl_function_pp_fkt([&](double r)->double{
        return  r*r * nbar(r) * sqrt(expLambda(r));
    });

    // Integration Method
    gsl_integration comp_B_int(comp_B_dBdr.gsl_fkt(),{0,R},comp_B_params[0],comp_B_params[1],"odeiv2",(size_t)comp_B_params[2],6,print);

    B   = M_4PI*comp_B_int.F;
    MB  = B*star_eos->mB;
    EB = MB-M;

    if(print){
        printf("|=> Baryon number B: %.16E | Baryonic Mass MB: %.16E MS\n",B,MB/MSkm);
    }
}

double Static_Star::M_error(int print, double err_rel, double err_abs) {
    double dM = comp_MGm2(print,err_rel,err_abs);
    return abs(1.-dM/M);
}
