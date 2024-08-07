//
// Created by M. J. Steil on 2017.05.12.
//

#include "../../include/ISS_O11_Series.hpp"

ISS_O11_Series::ISS_O11_Series(Static_Star *iss_in, int print) {
    iss = iss_in;

    M = &iss->M;
    R = &iss->R;
    Z = &iss->Z;
}

ISS_O11_Series::ISS_O11_Series(Static_Star_ISS *iss_in, int print) {
    iss_ana = iss_in;

    M = &iss_ana->M;
    R = &iss_ana->R;
    Z = &iss_ana->Z;
}

double ISS_O11_Series::omegabar0(double r, int d, int order) {
    double s= r/(*R);

    double wi_res[3];

    double f =0;

    if(r<=(*R)){
        for(int i = 0; i<= order; i++){
            gsl_poly_eval_derivs (wi[i], (size_t)i*2+1, s, wi_res, 3);
            f = f + wi_res[d]*gsl_pow_int((*Z),i)/gsl_pow_int((*R),d);
        }
    }else{
        if(d==0){
            f = 1.-2.*J0(order)/gsl_pow_3(r);
        }else if(d==1){
            f = 6.*J0(order)/gsl_pow_4(r);
        }else if(d==2){
            f = -24.*J0(order)/gsl_pow_5(r);
        }
    }
    return f;
}

double ISS_O11_Series::J0(int order) {
    double I0 = 2./5.*(*R)*(*R)*(*M);

    double f =0;
    for(int i = 0; i<= order; i++){
        f =f + wI[i]*gsl_pow_int((*Z),i-1); // Z**(i-1) because we divide out one Z with I0
    }

    return I0*f;
}

double ISS_O11_Series::mu0(int order) {
    double mu0f = -1./5.*(*R)*(*R)*(*R);

    double f =0;
    for(int i = 0; i<= order; i++){
        f =f + aMu[i]*gsl_pow_int((*Z),i);
    }

    return mu0f*f;
}

double ISS_O11_Series::B0(int order) {

    double f =0;
    for(int i = 0; i<= order; i++){
        f =f + aB0[i]*gsl_pow_int((*Z),i);
    }

    return f;
}

double ISS_O11_Series::aphi0(double r, int d, int order) {
    double s= r/(*R);

    double ai_res[2];

    double f =0;

    if(r<=(*R)){
        for(int i = 0; i<= order; i++){
            gsl_poly_eval_derivs (ai[i], (size_t)i*2+3, s, ai_res,2);
            f = f + ai_res[d]*gsl_pow_int((*Z),i);
        }
    }else{
        for(int i = 0; i<= order; i++){
            gsl_poly_eval_derivs (aEi[i], 20, s, ai_res, 2);

            if(d==0){
                f = f + ai_res[0]*gsl_pow_int((*Z)/s,i);
            }else if(d==1){
                f = f + -1.*i/s*ai_res[0]*gsl_pow_int((*Z)/s,i)+ai_res[1]*gsl_pow_int((*Z)/s,i);
            }

        }
    }
    return f*(*R)*(*R)/gsl_pow_int((*R),d);
}

double ISS_O11_Series::Btet0(double r, double theta, int i, int order) {
    switch (i){
        case 1: //Btet_r
            if(r>0){
                return -2*aphi0(r,0,order)/gsl_pow_2(r)*cos(theta);
            }else{
                return B0(order)*cos(theta);
            }
        case 2: //Btet_theta
            if(r>0){
                double f=0;
                double s=r/(*R);

                if(r<=(*R)){
                    for(int j = 0; j<= order; j++){

                        f = f + gsl_poly_eval(aBthetaI[j], 41, s)*gsl_pow_int((*Z),j);
                    }
                }else{
                    for(int j = 0; j<= order; j++){

                        f = f + gsl_poly_eval(aBthetaE[j], 20, s)*gsl_pow_int((*Z)/s,j)/(s*s);
                    }
                }

                return f*sin(theta);
            }else{
                return -B0(order)*sin(theta);
            }

        case 3: //Btet_phi
            return 0;
        default:
            GSL_ERROR_VAL ("Btet0(r,theta,i,order): Invalid input parameters i=1: B_r, i=2: B_theta, i=3: B_phi.",GSL_FAILURE,GSL_FAILURE);
    }
}

double ISS_O11_Series::Qe00(int order) {
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
