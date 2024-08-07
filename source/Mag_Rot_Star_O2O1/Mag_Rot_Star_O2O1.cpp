//
// Created by M. J. Steil on 2017.04.05.
//

#include "../../include/Mag_Rot_Star_O2O1.hpp"

Mag_Rot_Star_O2O1::Mag_Rot_Star_O2O1() {
}

Mag_Rot_Star_O2O1::Mag_Rot_Star_O2O1(Mag_Star_O2 * mag_star2_in) {
    mag_rot_star_O2O1_l1_status=0;
    //mag_rot_star_O2O1_l3_status=0;

    star = mag_star2_in->star;
    rot_star = mag_star2_in->mag_star->rot_star;
    mag_star1 = mag_star2_in->mag_star;
    mag_star2 = mag_star2_in;

    lisk_workspace = new LiSK::LiSK<complex<double>>(3);
}

void Mag_Rot_Star_O2O1::info() {
    printf("Mag_Rot_Star_O2O1\n");
    printf("|-> f = %.8E Hz | Omega = %.8E rad/s | P = %.8E s\n",rot_star->f/cHzkm,rot_star->Omega/cHzkm,1/rot_star->f*cHzkm);
    printf("|-> J = %.8E Msol**2 = %.8E kg m**2 s**-1\n",rot_star->J()/gsl_pow_2(MSkm),rot_star->J()/ckgm2s1km2);
    printf("|-> mu = %.8E km**2 = %.8E Msol**2 = %8E  Am**2\n",mag_star1->mu(),mag_star1->mu()/gsl_pow_2(MSkm),mag_star1->mu()/cAm2km2);

    printf("|=> deltaJ = %.8E Msol**2 = %.8E kg m**2 s**-1 = (%.4E)*J\n",deltaJ()/gsl_pow_2(MSkm),deltaJ()/ckgm2s1km2,deltaJ()/rot_star->J());
    printf("|=> W3asymp = %.8E Msol**4 = %.8E kg m**4 s**-1 \n",W3asymp()/gsl_pow_4(MSkm), W3asymp()/ckgm2s1km2*1E6);
    printf("\n");
}