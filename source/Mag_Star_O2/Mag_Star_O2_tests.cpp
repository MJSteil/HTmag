//
// Mag_Star_O2 test methods
// Created by m. J. Steil on 2016.10.10.
//

#include "../../include/Mag_Star_O2.hpp"


// l=2 tests

void Mag_Star_O2::test_isoSurface() {
    double r = 12;
    double theta = 0;
    double rb = rbar(r,theta);

    double H = (star->h(rb)+BB()*(h00(rb) + h20(rb)*(1. - 1.5 * gsl_pow_2(sin(theta)))))/star->h(r)-1;

    printf("theta=0:\t\t %.16E \n",H);

    theta = 0.5*M_PI;
    rb = rbar(r,theta);
    H = (star->h(rb)+BB()*(h00(rb) + h20(rb)*(1. - 1.5 * gsl_pow_2(sin(theta)))))/star->h(r)-1;

    printf("theta=0.5*M_PI:\t %.16E \n",H);
}

void Mag_Star_O2::test_Bernoulli(double rel_r, double theta) {
    double R = star->R;
    double r = rel_r*R;
    double X = sin(theta)*sin(theta);
    double P2= 1-1.5*X;

    double dBO0 = star->dh(r)+star->dNu(r)/2;
    double BO0 = star->h(r)+log(star->expNu(r))/2;
    double cNu = star->h(0)+log(star->expNu(0))/2;

    printf("dBernoulli O(B**0) @ r=%.2ER: %.16E \n",r/R,dBO0);
    printf("Bernoulli O(B**0) @ r=%.2ER: %.16E \n",r/R,1-BO0/cNu);

    double dBO2 = dh00(r)+dh20(r)*P2+dn00(r)+dn20(r)*P2-mag_star->cjphi0*mag_star->daphi0(r)*X;
    double BO2 = h00(r)+h20(r)*P2+n00(r)+n20(r)*P2-mag_star->cjphi0*mag_star->aphi0(r)*X;
    printf("dBernoulli O(B**2) @ (r=%.2E R,theta=%.2E M_PI): %.16E \n",r/R,theta/M_PI,dBO2);
    printf("Bernoulli O(B**2) @ (r=%.2E R,theta=%.2E M_PI): %.16E \n",r/R,theta/M_PI,1-BO2/ch00);
}