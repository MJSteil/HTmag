//
// Global parameters of Mag_Star_O2
// Created by M. J. Steil on 2017.02.27.
//

#include "../../include/Mag_Star_O2.hpp"

double Mag_Star_O2::BB() {
    return gsl_pow_2(mag_star->eA);
}

double Mag_Star_O2::deltaM() {
    return deltaM0*BB();
}

double Mag_Star_O2::deltaB(int print) {
    return MBR/star->star_eos->mB*BB();
}

double Mag_Star_O2::deltaMB(int print) {
    return MBR*BB();
}

double Mag_Star_O2::deltaR_l0() {
    return xi00(star->R)*BB();
}

double Mag_Star_O2::deltaA() {
    return 8*M_PI*deltaR_l0()*star->R;
}

double Mag_Star_O2::A_r(double r) {
    return M_4PI*(2.*xi00(r)*BB()+r)*r;
}


double Mag_Star_O2::Qm() {
    return BB()*Qm0;
}

double Mag_Star_O2::rsph(double r, double theta) {
    if (r > 0) {
        return r + BB()*( xi00(r) + (xi20(r) + r * (v20(r) - n20(r))) * (1. - 1.5 * gsl_pow_2(sin(theta))) );
    } else {
        return 0;
    }
}

double Mag_Star_O2::Bmax() {
    double R = star->R;
    return sqrt(abs(-R/( xi00(R) + (xi20(R) + R * (v20(R) - n20(R))) )));
}

double Mag_Star_O2::Bth(double thin) {
    double R = star->R;

    return sqrt(abs(2.*R*thin/( ( R*(v20(R)-n20(R)) -2*xi00(R) + xi20(R) ) )));
}


double Mag_Star_O2::e() {
    double R = star->R;

    return sqrt(3.*(n20(R)-v20(R)-xi20(R)/R)*BB());
}

double Mag_Star_O2::epsilon() {
    double R = star->R;

    return -1.5*(-n20(R)+v20(R)+xi20(R)/R)*BB();
}

double Mag_Star_O2::epsilonErr() {
    double R = star->R;

    return abs(1-epsilon()/(1-rsph(star->R,0)/rsph(star->R,M_PI/2)));
}

void Mag_Star_O2::eSummary() {
    double R = star->R;

    double RpRe = 1 + 1.5*(-n20(R)+v20(R)+xi20(R)/R)*BB();

    printf("e: %.8E \n",e());
    printf("epsilon: %.8E (%.1E)\n",epsilon(),epsilonErr());
    printf("RpRe: %.8E \n",RpRe);
}

double Mag_Star_O2::deltaJ20() {
    return deltaJ20_00*mag_star->Qe()*mag_star->B0;
}

double Mag_Star_O2::deltaMQ() {
    if(mag_star->electric_l0==1){
        return gsl_pow_2(mag_star->Qe())/star->R*(2.-3.*star->Z)/(2.-2.*star->Z);
    }else{
        return 0;
    }
}
