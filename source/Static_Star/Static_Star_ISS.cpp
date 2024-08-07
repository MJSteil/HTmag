//
// Created by M. J. Steil on 2016.09.09.
//

#include "../../include/Static_Star_ISS.hpp"

using namespace std;
using namespace units;

Static_Star_ISS::Static_Star_ISS() {

}

Static_Star_ISS::Static_Star_ISS(double rhoCin) {
    rhoC = rhoCin;
    eos_irf eos_in(rhoCin);
    star_eos = &eos_in;
}

Static_Star_ISS::Static_Star_ISS(double Min /*[km]*/, double Rin /*[km]*/) {
    setMR(Min,Rin);
}

void Static_Star_ISS::sethCrhoC(double hC /*[1]*/, double rhoC /*[km**-2]*/) {

    R = sqrt((3.-9.*exp(hC)+6.*exp(2.*hC))/(2.*M_PI*rhoC))/(3.*exp(hC)-2.);
    M = 4.*M_PI/3.*gsl_pow_3(R)*rhoC;

    setMR(M,R);
}

void Static_Star_ISS::setMR(double Min /*[km]*/, double Rin /*[km]*/) {
    M = Min;
    R = Rin;
    rhoC = Min/(4.*M_PI/3.*gsl_pow_3(Rin));

    eos_irf eos_in(rhoC);
    star_eos = &eos_in;

    setAux();
}

void Static_Star_ISS::setZ(double Zin) {
    R = sqrt(Zin/(M_4PI/3.*rhoC));
    M = Zin*R;

    eos_irf eos_in(rhoC);
    star_eos = &eos_in;

    setAux();
}



void Static_Star_ISS::setAux() {
    Z=M/R;

    hC = log( ( 3. - 6.*Z +sqrt(1.-2.*Z) )/( 4.-9.*Z ) );
    PC= star_eos->P(hC);

    twoMR=2*M/R;
    twoMR3=2*M/gsl_pow_3(R);

    MB = 3./4.*R*(asin(sqrt(twoMR))/sqrt(twoMR)-sqrt(1-twoMR));
}

double Static_Star_ISS::P(double r) {
    if(r<R){
        return rhoC*(sqrt(1-twoMR3*gsl_pow_2(r)) - sqrt(1-twoMR))/(3*sqrt(1-twoMR) -sqrt(1-twoMR3*gsl_pow_2(r)));
    }else{
        return 0;
    }
}

double Static_Star_ISS::dP(double r) {
    if(r<R){
        return -3.*M*M/(M_PI*gsl_pow_6(R))* r*sqrt(1-twoMR)/sqrt(1-twoMR3*r*r)/gsl_pow_2(sqrt(1-twoMR3*r*r)-3*sqrt(1-twoMR));
    }else{
        return 0;
    }
}

double Static_Star_ISS::rho(double r) {
    if(r<=R){
        return rhoC;
    }else{
        return 0;
    }
}

double Static_Star_ISS::drhodh(double r) {
    return 0;
}

double Static_Star_ISS::h(double r) {
    if(r<=R){
        double nu = log(0.25*gsl_pow_2(sqrt(1-twoMR3*gsl_pow_2(r))-3*sqrt(1-twoMR)));
        return 0.5*(-nu+log(1-twoMR));
    }else{
        return 0;
    }
}

double Static_Star_ISS::dh(double r) {
    if(r<=R){
        return 2*M*r/(-2*M*r*r+R*R*R*(1-3*sqrt(1-twoMR3*r*r)*sqrt(1-twoMR)));
    }else{
        return 0;
    }
}

double Static_Star_ISS::nbar(double r) {
    return rho(r)/star_eos->mB;
}

double Static_Star_ISS::dnbardh(double r) {
    return 0;
}

// Static_Star metric potentials (g_tt=-exp(nu), g_rr=exp(lambda)) and their derivatives

double Static_Star_ISS::expLambda(double r) {

    if(r<=R&&r>0){
        return 1./(1-2*M/gsl_pow_3(R)*gsl_pow_2(r));
    }else if (r==0){
        return 1.;
    }else{
        return 1./(1-2*M/r);
    }
}

double Static_Star_ISS::dLambda(double r) {
    if(r<=R){
        return -4.*M*r/(2.*M*r*r-gsl_pow_3(R));
    }else{
        return 2.*M/(2.*M*r-r*r);
    }
}

double Static_Star_ISS::expNu(double r) {
    if(r<=R&&r>0){
        return 0.25*gsl_pow_2(sqrt(1-twoMR3*gsl_pow_2(r))-3*sqrt(1-twoMR));
    }else if (r==0){
        return 0.25*gsl_pow_2(1-3*sqrt(1-twoMR));
    }else{
        return 1-2*M/r;
    }
}

double Static_Star_ISS::dNu(double r) {
    if(r<=R){
      return -2*dh(r);
    }else{
        return -2*M/(2*M*r-r*r);
    }
}

double Static_Star_ISS::m(double r) {
    if(r<=R){
        return M*gsl_pow_3(r/R);
    }else{
        return M;
    }
}

double Static_Star_ISS::dm(double r){
    if(r<=R){
        return 3*M*gsl_pow_2(r/R)/R;
    }else{
        return 0;
    }
}

void Static_Star_ISS::show_limits() {
    printf("Schwarzschild/Buchdahl-Limit - M/R=4/9:\n");
    double zr = 4./9.;
    double Rr = sqrt(zr/(M_4PI/3.*rhoC));
    double Mr = zr*Rr;
    double MBr = 3./4.*Rr*(asin(sqrt(2*zr))/sqrt(2*zr)-sqrt(1-2*zr));
    printf("|-> M=%.17g MS | R=%.17g km | M=%.17g MS  \n", Mr/MSkm,Rr,MBr/MSkm);

    printf("'Causality'-Limit - M/R=3/8:\n");
    zr = 3./8.;
    Rr = sqrt(zr/(M_4PI/3.*rhoC));
    Mr = zr*Rr;
    MBr = 3./4.*Rr*(asin(sqrt(2*zr))/sqrt(2*zr)-sqrt(1-2*zr));
    printf("|-> M=%.17g MS | R=%.17g km | M=%.17g MS  \n", Mr/MSkm,Rr,MBr/MSkm);

    printf("Tr(T_fluid)=0<=>P=1/3rho-Limit - M/R=5/18:\n");
    zr = 5./18.;
    Rr = sqrt(zr/(M_4PI/3.*rhoC));
    Mr = zr*Rr;
    MBr = 3./4.*Rr*(asin(sqrt(2*zr))/sqrt(2*zr)-sqrt(1-2*zr));
    printf("|-> M=%.17g MS | R=%.17g km | M=%.17g MS  \n", Mr/MSkm,Rr,MBr/MSkm);


    printf("\n");
}

void Static_Star_ISS::show_discontinuity() {
    printf("[lambda']=%.8E\n",8*M_PI*R*expLambda(R)*rhoC);
    printf("[nu'']=%.8E\n",8*M_PI*(1+R*dNu(R)/2.)*expLambda(R)*rhoC);
}



