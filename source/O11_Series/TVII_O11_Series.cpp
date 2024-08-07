//
// Created by M. J. Steil on 2017.05.12.
//

#include "../../include/TVII_O11_Series.hpp"

TVII_O11_Series::TVII_O11_Series(Static_Star *tvii_in, int mu_in, int print) {
    tvii = tvii_in;

    M = &tvii->M;
    R = &tvii->R;
    Z = &tvii->Z;

    muI = mu_in;
}


double TVII_O11_Series::omegabar0c( int order) {
    double f =0;
    for(int i = 0; i<= order; i++){
        f =f + womegabar0c[muI][i]*gsl_pow_int((*Z),i);
    }
    return f;
}

double TVII_O11_Series::J0(int order) {
    double I0 = (*R)*(*R)*(*M);

    double f =0;
    for(int i = 0; i<= order; i++){
        f =f + wI[muI][i]*gsl_pow_int((*Z),i-1); // Z**(i-1) because we divide out one Z with I0
    }

    return I0*f;
}

double TVII_O11_Series::mu0(int order) {
    double mu0f = -1.*(*R)*(*R)*(*R);

    double f =0;
    for(int i = 0; i<= order; i++){
        f =f + muvec[muI][i]*gsl_pow_int((*Z),i);
    }

    return mu0f*f;
}

double TVII_O11_Series::B0(int order) {

    double f =0;
    for(int i = 0; i<= order; i++){
        f =f + B0vec[muI][i]*gsl_pow_int((*Z),i);
    }

    return -f;
}

double TVII_O11_Series::BR(int order) {

    double f =0;
    for(int i = 0; i<= order; i++){
        f =f + BPvec[muI][i]*gsl_pow_int((*Z),i);
    }

    return -f;
}


double TVII_O11_Series::Qe00(int order) {
/*    double Qef = gsl_pow_5((*R));

    double f =0;
    for(int i = 0; i<= order; i++){
        f =f + QeInfi[i]*gsl_pow_int((*Z),i);
    }*/

    double Y = 1.-2.*(*Z);
    double logY = log(Y);

    double g =( (*M)*(*M)*(*M)*mu0(order)*(0.6 + 0.4*logY + (-0.8 + 0.2*Y)*Y)
      + J0(order)*mu0(order)*(-0.6 - 0.1*logY + Y*(-4.1 - 2.7*logY + Y*(3.9 - 3.3*logY + (0.9 + 0.1*logY - 0.1*Y)*Y)))
    )/((*M)*(-1. + Y*(-9. - 6.*logY + Y*(9. - 6.*logY + Y))));

    return g;//Qef*f;
}
