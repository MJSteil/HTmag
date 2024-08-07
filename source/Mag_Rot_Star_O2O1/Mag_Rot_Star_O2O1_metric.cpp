//
// Mag_Rot_Star_O2O1 metric perturbations
// Created by M. J. Steil on 2017.04.05.
//

#include "../../include/Mag_Rot_Star_O2O1.hpp"

double Mag_Rot_Star_O2O1::W10_ext(double r, double deltaJ0in, double part) {
    //region Variables
    double M = star->M;
    double M2 = M*M;
    double M3 = M2*M;
    double M4 = M3*M;
    double M5 = M4*M;
    double M6 = M5*M;
    double M7 = M6*M;

    double y = 1.-2.*M/r;
    double y2 = y*y;
    double y3 = y2*y;
    double y4 = y3*y;
    double y5 = y4*y;
    double y6 = y5*y;

    double logy = log(y);
    double logylogy = logy*logy;
    double ym1 = y-1.;
    double ym1_2 = ym1*ym1;
    double ym1_3 = ym1_2*ym1;
    double ym1_4 = ym1_3*ym1;

    double mu0 = mag_star1->mu0;
    double mu0mu0 = gsl_pow_2(mu0);

    double J0= rot_star->J0;

    double QE20 = mag_star1->QEr00;

    double QM20 = mag_star2->Qm0;
    //endregion

    double W10_ext_H = 0.;
    double W10_ext_P = 0.;

    if(deltaJ0in!=0.){
        W10_ext_H = -deltaJ0in*ym1_3/4./M3; // = 2*deltaJ0in/r**3
    }

    if(part!=0.){
        W10_ext_P = -3*J0*QM20*(ym1*(3. - 3.*y2 + (1. + 4.*y + y2)*logy))/(8.*M6)
                     -QE20*mu0*(ym1_2*(11. - 69.*y - 15.*y2 + y3) + 12.*y*(-7 + 2*y + 5*y2)*logy - 36*y*(1 + y)*logylogy)/(8.*M6*ym1_2)
                     +J0*mu0mu0*(    ym1_2*(10133. - 24240.*y - 3900.*y2 + 3400.*y3 + 135.*y4 + 72.*y5)
                                     - 60.*(-54 + 337*y + 19*y2 - 334*y3 + 19*y4 + 7*y5 + 6*y6)*logy
                                     + 360.*(1- 26*y - 13*y2 + y3 - 4*y4 + y5)*logylogy
        )/(3200.*M7*ym1_2);
    }
    return W10_ext_H + W10_ext_P;
}

double Mag_Rot_Star_O2O1::dW10_ext(double r, double deltaJ0in, double part) {
    //region Variables
    double M = star->M;
    double M2 = M*M;
    double M3 = M2*M;
    double M4 = M3*M;
    double M5 = M4*M;
    double M6 = M5*M;
    double M7 = M6*M;

    double y = 1.-2.*M/r;
    double y2 = y*y;
    double y3 = y2*y;
    double y4 = y3*y;
    double y5 = y4*y;
    double y6 = y5*y;

    double logy = log(y);
    double logylogy = logy*logy;
    double ym1 = y-1.;
    double ym1_2 = ym1*ym1;
    double ym1_3 = ym1_2*ym1;
    double ym1_4 = ym1_3*ym1;

    double mu0 = mag_star1->mu0;
    double mu0mu0 = gsl_pow_2(mu0);

    double J0= rot_star->J0;

    double QE20 = mag_star1->QEr00;

    double QM20 = mag_star2->Qm0;
    //endregion

    double dW10_ext_H = 0.;
    double dW10_ext_P = 0.;

    if(deltaJ0in!=0.){
        dW10_ext_H = (-3*ym1_4*deltaJ0in)/(8.*M4); // = -6*deltaJ0in/r**4
    }

    if(part!=0.){
        dW10_ext_P = (-3*ym1_2*J0*QM20*(-1 + 9*y2 + 3*logy*y*(-1 + 2*y + y2) - 8*y3))/(16.*M7*y)
                     +(-3*mu0*QE20*(12*logylogy*(1 + 3*y) + ym1_2*(51 + 7*y - 11*y2 + y3) + 4*logy*(13 + 3*y - 21*y2 + 5*y3)))/(16.*M7*ym1)
                     +(3*J0*mu0mu0*(6*logylogy*y*(24 + 52*y - 3*y2 + 17*y3 - 13*y4 + 3*y5)
                                    + ym1_2*(-54 + 687*y + 28*y2 - 332*y3 + 148*y4 - 3*y5 + 6*y6)
                                    + logy*(-12 + 553*y - 24*y5*y2 + 219*y2 - 1170*y3 + 470*y4 - 63*y5 + 27*y6)))/(320.*M7*M*y*ym1);

    }
    return dW10_ext_H + dW10_ext_P;
}



double Mag_Rot_Star_O2O1::W30_ext(double r, double W3asymp0in, double part){
    //region Variables
    double M = star->M;
    double M2 = M*M;
    double M3 = M2*M;
    double M4 = M3*M;
    double M5 = M4*M;
    double M6 = M5*M;
    double M7 = M6*M;

    double y = 1.-2.*M/r;
    double y2 = y*y;
    double y3 = y2*y;
    double y4 = y3*y;
    double y5 = y4*y;
    double y6 = y5*y;

    double logy = log(y);
    double logylogy = logy*logy;
    double ym1 = y-1.;
    double ym1_2 = ym1*ym1;
    double ym1_3 = ym1_2*ym1;
    double ym1_4 = ym1_3*ym1;

    double mu0 = mag_star1->mu0;
    double mu0mu0 = gsl_pow_2(mu0);

    double J0= rot_star->J0;

    double QE20 = mag_star1->QEr00;

    double QM20 = mag_star2->Qm0;
    //endregion

    double W30_ext_H = 0.;
    double W30_ext_P = 0.;

    if(W3asymp0in!=0.){
        W30_ext_H = (-7*W3asymp0in*(-6 - 125*y - 60*logy*y*(1 + 2*y) + 80*y2 + 60*y3 - 10*y4 + y5))/(64.*M5*ym1_2);
    }

    if(part!=0.){
        W30_ext_P = -(mu0*QE20*(6*logy*y*(32 + logy*(3 - 12*y) + 58*y + 15*y2) + (-3 + y*(-481 + y*(-129 + (-18 + y)*y)))*ym1))/(24.*M6*ym1_2)
                    +(J0*QM20*(12*logy*(-6 - 30*logy*y*(1 + 2*y) + y*(49 + 2*y*(214 + 3*y*(12 + (-4 + y)*y)))) - (-6 + y*(3667 + y*(3315 + y*(-837 + 161*y))))*ym1))/(192.*M6*ym1_2)
                    + (3*J0*logy*mu0mu0*(6 + 131*y + 51*y2 - 9*y3 + y4)*log(1. - y))/(20.*M7*ym1)
                    + (3*J0*mu0mu0*(-6 - 125*y + 60*logy*y*(1 + 2*y) + 80*y2 + 60*y3 - 10*y4 + y5)*real(lisk_workspace->Li2(y)))/(20.*M7*ym1_2)
                    + (-18*J0*mu0mu0*y*(1 + 2*y)*real(lisk_workspace->Li3(y)))/(M7*ym1_2)
                    + -(J0*mu0mu0*(59460 + 1148279*y + 432*y6*y + 7200*logy*logylogy*y*(1 + 2*y) - 359942*y2 - 972960*y3 + 142780*y4 + 1080*logylogy*(2 + 33*y - 6*y2 - 18*y3 + 2*y4) - 19975*y5 + 480*M_PI_SQARE*(-6 - 125*y + 80*y2 + 60*y3 - 10*y4 + y5) + 1926*y6 - 60*logy*(-234 + (-7186 + 480*M_PI_SQARE)*y + 10*(-2045 + 96*M_PI_SQARE)*y2 - 6480*y3 + 1062*y4 - 159*y5 + 36*y6) - 345600*y*zeta3 - 691200*y2*zeta3))/(19200.*M7*ym1_2);
    }

    return W30_ext_H + W30_ext_P;
}

double Mag_Rot_Star_O2O1::dW30_ext(double r, double W3asymp0in, double part){
    //region Variables
    double M = star->M;
    double M2 = M*M;
    double M3 = M2*M;
    double M4 = M3*M;
    double M5 = M4*M;
    double M6 = M5*M;
    double M7 = M6*M;

    double y = 1.-2.*M/r;
    double y2 = y*y;
    double y3 = y2*y;
    double y4 = y3*y;
    double y5 = y4*y;
    double y6 = y5*y;

    double logy = log(y);
    double logylogy = logy*logy;
    double ym1 = y-1.;
    double ym1_2 = ym1*ym1;
    double ym1_3 = ym1_2*ym1;
    double ym1_4 = ym1_3*ym1;

    double mu0 = mag_star1->mu0;
    double mu0mu0 = gsl_pow_2(mu0);

    double J0= rot_star->J0;

    double QE20 = mag_star1->QEr00;

    double QM20 = mag_star2->Qm0;
    //endregion

    double dW30_ext_H = 0.;
    double dW30_ext_P = 0.;

    if(W3asymp0in!=0.){
        dW30_ext_H = (-7*W3asymp0in*(197 + 25*y + 60*logy*(1 + 5*y) - 300*y2 + 100*y3 - 25*y4 + 3*y5))/(128.*M6*ym1);
    }

    if(part!=0.){
        dW30_ext_P = -(mu0*QE20*(-676 + 70*y + 18*logylogy*(-1 + 7*y) + 591*y2 + 55*y3 + 6*logy*(-38 - 118*y - 69*y2 + 15*y3) - 43*y4 + 3*y5))/(48.*M7*ym1)
                     + (J0*QM20*(72 - 4321*y + 360*logylogy*y*(1 + 5*y) - 7517*y2 + 16728*y3 - 6992*y4 + 2441*y5 + 12*logy*y*(23 - 845*y - 336*y2 + 168*y3 - 78*y4 + 18*y5) - 411*y6))/(384.*M7*y*ym1)
                     + (3*J0*logy*mu0mu0*(-197 - 222*y + 78*y2 - 22*y3 + 3*y4)*log(1 - y))/(40.*M7*M)
                     + (3*J0*mu0mu0*(197 + 25*y - 60*logy*(1 + 5*y) - 300*y2 + 100*y3 - 25*y4 + 3*y5)*real(lisk_workspace->Li2(y)))/(40.*M7*M*ym1)
                     + (9*J0*mu0mu0*(1 + 5*y)*real(lisk_workspace->Li3(y)))/(M7*M*ym1)
                     + (J0*mu0mu0*(14040 - 2520*y6*y - 2160*y6*y2 + 7200*logylogy*logy*y*(1 + 5*y) + 120*(-31309 + 1200*M_PI_SQARE)*y3 - 40*(-27289 + 1200*M_PI_SQARE)*y4 - 1080*logylogy*y*(-57 - 41*y + 94*y2 - 26*y3 + 4*y4) + 25*(-12487 + 480*M_PI_SQARE)*y5 + 60*logy*(72 + (8482 - 480*M_PI_SQARE)*y + 144*y6*y + (40682 - 2400*M_PI_SQARE)*y2 + 22848*y3 - 7128*y4 + 2367*y5 - 645*y6) + (59781 - 1440*M_PI_SQARE)*y6 + y*(1684319 - 94560*M_PI_SQARE - 345600*zeta3) - 5*y2*(-244847 + 2400*M_PI_SQARE + 345600*zeta3)))/(38400.*M7*M*y*ym1);

    }

    return dW30_ext_H + dW30_ext_P;
}

double Mag_Rot_Star_O2O1::W10(double r) {
    if (mag_rot_star_O2O1_l1_status >= 2) {
        if(r==0.){
            return v10_i[0]*star->j(r);
        } else if (r <= star->R) {
            return v10_of_r.f(r)*star->j(r);
        } else {
            return W10_ext(r, deltaJ0, 1);
        }
    } else {
        GSL_ERROR_VAL ("W10: No interpolations supplied: run comp_l10_int().", GSL_FAILURE, GSL_FAILURE);
    }
}

double Mag_Rot_Star_O2O1::dW10(double r) {
    if (mag_rot_star_O2O1_l1_status >= 2) {
        if (r <= star->R&&r>0) {
            return u10_of_r.f(r)*star->j(r)/gsl_pow_4(r);
        } else if(r > star->R) {
            return dW10_ext(r, deltaJ0, 1);
        } else{
            return 0;
        }
    } else {
        GSL_ERROR_VAL ("dW10: No interpolations supplied: run comp_l10_int().", GSL_FAILURE, GSL_FAILURE);
    }
}

double Mag_Rot_Star_O2O1::W30(double r) {
    if (mag_rot_star_O2O1_l3_status >= 2) {
        if(r==0.){
            return v30_i[0]*star->j(r);
        } else if (r <= star->R) {
            return v30_of_r.f(r)*star->j(r);
        } else {
            return W30_ext(r, W3asymp0, 1);
        }
    } else {
        GSL_ERROR_VAL ("W30: No interpolations supplied: run comp_l30_int().", GSL_FAILURE, GSL_FAILURE);
    }
}

double Mag_Rot_Star_O2O1::dW30(double r) {
    if (mag_rot_star_O2O1_l3_status >= 2) {
        if (r <= star->R&&r>0) {
            return u30_of_r.f(r)*star->j(r)/gsl_pow_4(r);
        } else if(r > star->R) {
            return dW30_ext(r, W3asymp0, 1);
        } else{
            return 0;
        }
    } else {
        GSL_ERROR_VAL ("dW30: No interpolations supplied: run comp_l30_int().", GSL_FAILURE, GSL_FAILURE);
    }
}