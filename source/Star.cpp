//
// Created by M. J. Steil on 2017.04.26.
//

#include "../include/Star.hpp"

Star::Star(eos &eos_in) {
    eos *typeEoS = &eos_in;
    star = Static_Star(*typeEoS);

    rot_star = Rot_Star_O1(&star);

    mag_star_o1 = Mag_Star_O1(&star);
    mag_star_o1.set_rot(&rot_star);

    mag_star_o2 = Mag_Star_O2(&mag_star_o1);

    mag_rot_star_o2o1 = Mag_Rot_Star_O2O1(&mag_star_o2);
}

Star::Star(eos &eos_in, double hC, double B0, double f, double Q, int level) {
    eos *typeEoS = &eos_in;
    star = Static_Star(*typeEoS);

    rot_star = Rot_Star_O1(&star);

    mag_star_o1 = Mag_Star_O1(&star);
    mag_star_o1.set_rot(&rot_star);

    mag_star_o2 = Mag_Star_O2(&mag_star_o1);

    mag_rot_star_o2o1 = Mag_Rot_Star_O2O1(&mag_star_o2);

    comp(hC, B0, f, Q, level);
}

void Star::comp(double hC, double B0, double f, double Q, int level) {
    if(level>=1){
        star.comp(hC,0);
    }

    if(level>=2&&f>0.){
        rot_star.set_f(f);
        rot_star.comp(1,0);
    }

    if(level>=3){
        mag_star_o1.set_B0Q0(B0,Q);
        if(mag_star_o1.magnetic_l1){
            mag_star_o1.comp(B0,1,0);
        }
        if(f>0){
            mag_star_o1.set_rot(&rot_star);
        }
    }

    if(level>=4){
        mag_star_o2.comp_l00(1,0);
        mag_star_o2.comp_l00_int();
        if(mag_star_o1.electric_l0==1){
            mag_star_o2.comp_l100(1,0);
            mag_star_o2.comp_l100_int();
        }
        mag_star_o2.comp_l20(1,0);
        mag_star_o2.comp_l20_int();
    }

    if(level>=5){
        mag_rot_star_o2o1.comp_l10(1,0);
        mag_rot_star_o2o1.comp_l10_int();

        mag_rot_star_o2o1.comp_l30(1,0);
        mag_rot_star_o2o1.comp_l30_int();
    }
}

void Star::Jrot_vs_JQB() {
    double Irot = rot_star.I;
    double JQB = mag_star_o2.deltaJ20();

    printf("%.8E \n", JQB/Irot/cHzkm);
}


