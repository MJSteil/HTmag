//
// Created by M. J. Steil on 2016.07.06.
//

#include "../../include/Rot_Star_O1.hpp"
Rot_Star_O1::Rot_Star_O1() {

}

Rot_Star_O1::Rot_Star_O1(Static_Star * star_in, string name) {
    rot_star_status=0;
    star = star_in;
}

Rot_Star_O1::Rot_Star_O1(Static_Star * star_in, double f_in, int print) {
    rot_star_status=0;
    star = star_in;

    set_f(f_in);

    comp(1,print);
}

void Rot_Star_O1::set_f(double f_in) {
    f = f_in;
    Omega = M_2PI*f;
}

void Rot_Star_O1::info() {
    printf("Rot_Star_O1\n");
    printf("|-> f = %.8E Hz | Omega = %.8E rad/s | P = %.8E s\n",f/cHzkm,Omega/cHzkm,1/f*cHzkm);
    printf("|-: %d nodes (%d steps,n/s=%.2g)\n",n,steps,1.*n/steps);
    printf("|=> omega(0)/Omega = %.8E ,omega(R)/Omega = %.8E \n",omega0(0),omega0(star->R));
    printf("|=> J = %.8E Msol**2 = %.8E kg m**2 s**-1\n",J()/gsl_pow_2(MSkm),J()/ckgm2s1km2);
    printf("|=> I = %.8E Msol**3 = %.8E M*R**2 =%.16E kg m**2 \n",I/gsl_pow_3(MSkm),I/star->M/gsl_pow_2(star->R),I/ckgm2km3);
    printf("\n");
}