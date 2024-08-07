//
// Created by m. J. Steil on 2016.08.11.
//

#include "../../include/Mag_Star_O2.hpp"

Mag_Star_O2::Mag_Star_O2(){

}

Mag_Star_O2::Mag_Star_O2(Mag_Star_O1 * mag_star) {
    mag_star_O2_l0_status=0;
    mag_star_O2_l2_status=0;
    this->mag_star = mag_star;
    this->star = mag_star->star;
}

void Mag_Star_O2::info() {
    //comp_deltaB0(0);

    printf("Mag_Star_O2:\n");
    printf("|=> deltaM  = %.8E MS = %.8E*M_0 \t=> M=\t%.16E MS\n",deltaM()/MSkm,deltaM()/star->M,(star->M+deltaM())/MSkm);
    printf("|=> deltaMB = %.8E MS = %.8E*MB_0 \t=> MB=\t%.16E MS \n",deltaMB()/MSkm,deltaMB()/star->MB,(star->MB+deltaMB())/MSkm);
    printf("|=> deltahC = %.8E <=> hC = %.8E\n",h00_i[0]*BB(),h00_i[0]*BB()+star->hC);
    printf("|=> QM = %.8E kg m**2 \n",Qm()/ckgm2km3);

    printf("|=> deltaA=%.8E km**2 = %.8E*A_0 \n",deltaA(),deltaA()/(M_4PI*star->R*star->R));
    printf("|=> R_mean:\t %.8E (=Rs+x: x=%.8E) km\n",star->R+deltaR_l0(),deltaR_l0());

    if(mag_star_O2_l2_status>1){
        printf("|=> R_e:\t %.8E (=Rs+x: x= %.8E) km\n",rsph(star->R,M_PI/2),rsph(star->R,M_PI/2)-star->R);
        printf("|=> R_p:\t %.8E (=Rs+x: x=%.8E) km\n",rsph(star->R,0),rsph(star->R,0)-star->R);
        printf("|=> epsilon:\t %.4E (rel. error= %.4E)\n",epsilon(),epsilonErr());
    }

    printf("\n");

    if(mag_star->electric_l0==1){
        printf("|=> deltaJ_QB = %.8E Msol**2 = %.8E kg m**2 s**-1 <- (%.4E GT|%.4E C)\n",deltaJ20()/gsl_pow_2(MSkm),deltaJ20()/ckgm2s1km2,mag_star->B0/cGTkm,mag_star->Qe()/cCkm);
        printf("|=> deltaMQ = %.8E Msol \n",deltaMQ());
    }

}
