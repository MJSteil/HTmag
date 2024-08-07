//
// Created by M. J. Steil on 2017.07.14.
//

#include "../../include/Mag_Rot_Star_O2O1.hpp"

double Mag_Rot_Star_O2O1::S0(double r) {
    double omegabar0 = rot_star->omegabar0(r);
    double domegabar0 = rot_star->domegabar0(r);

    double expLambda = star->expLambda(r);
    double dLambda = star->dLambda(r);
    double dNu = star->dNu(r);

    if(r>0){
        return r*r/star->j(r)*(
                r*r*mag_star2->dn00(r)*domegabar0
                + expLambda*mag_star2->m00(r)*( 4*(dNu+dLambda)*omegabar0 + (r*dLambda-1)*domegabar0 )
                + expLambda*mag_star2->dm00(r)*r*domegabar0
                + 16*M_PI*expLambda*r*r*mag_star2->h00(r)*( star->drhodh(r) + star->P(r) + star->rho(r) )*omegabar0
        );
    }else{
        return 0;
    }

}

double Mag_Rot_Star_O2O1::S2(double r) {
    double omegabar0 = rot_star->omegabar0(r);
    double domegabar0 = rot_star->domegabar0(r);

    double expLambda = star->expLambda(r);
    double dLambda = star->dLambda(r);
    double dNu = star->dNu(r);

    if(r>0){
        return r*r/star->j(r)*(
                r*r*domegabar0*( mag_star2->dn20(r) - 4./5.*mag_star2->dv20(r) )
                + 0.2*expLambda*mag_star2->m20(r)*( 4*(dNu+dLambda)*omegabar0 + (r*dLambda-1)*domegabar0 )
                + 0.2*expLambda*mag_star2->dm20(r)*r*domegabar0
                + 16/5.*M_PI*expLambda*r*r*mag_star2->h20(r)*( star->drhodh(r) + star->P(r) + star->rho(r) )*omegabar0
        );
    }else{
        return 0;
    }

}

double Mag_Rot_Star_O2O1::S1(double r) {
    double expLambda = star->expLambda(r);

    if(r>0&&r<=star->R){
        return 16/5./star->j(r)*rot_star->omegabar0(r)*( expLambda*gsl_pow_2(mag_star1->aphi0(r)) + r*r*gsl_pow_2(mag_star1->daphi0(r)) );
    }else if(r==0. ){
        return 0;
    } else{
        return 16/5./star->j(r)*(
                -( expLambda*gsl_pow_2(mag_star1->aphi0(r)) + r*r*gsl_pow_2(mag_star1->daphi0(r)) )*rot_star->omega0(r)
                -1.5*expLambda*mag_star1->at20(r)*mag_star1->aphi0(r)
                +0.25*r*r*( 5*mag_star1->dat00(r)  -mag_star1->dat20(r) )*mag_star1->daphi0(r)
        );
    }

}

double Mag_Rot_Star_O2O1::S3(double r) {

    double expLambda = star->expLambda(r);

    if(r>0&&r<=star->R){
        return 8/15./star->j(r)*rot_star->omegabar0(r)*( 4*expLambda*gsl_pow_2(mag_star1->aphi0(r)) - r*r*gsl_pow_2(mag_star1->daphi0(r)) );
    }else if(r==0.){
        return 0;
    }else{
        return 8/15./star->j(r)*(
            ( -4.*expLambda*gsl_pow_2(mag_star1->aphi0(r)) + r*r*gsl_pow_2(mag_star1->daphi0(r)) )*2/3.*rot_star->omega0(r)
            -4*expLambda*mag_star1->at20(r)*mag_star1->aphi0(r)
            +r*r*mag_star1->dat20(r)*mag_star1->daphi0(r)
        );
    }

}