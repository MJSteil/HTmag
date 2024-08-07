//
// Created by M. J. Steil on 2017.04.26.
//

#ifndef HTMAG_STAR_HPP
#define HTMAG_STAR_HPP

#include "Static_Star.hpp"
#include "Rot_Star_O1.hpp"
#include "Mag_Star_O1.hpp"
#include "Mag_Star_O2.hpp"
#include "Mag_Rot_Star_O2O1.hpp"

class Star {
public:
    Star(eos &eos_in);
    Star(eos &eos_in, double hC, double B0, double f, double Q=0, int level =5);

    Static_Star star = Static_Star();

    Rot_Star_O1 rot_star = Rot_Star_O1();
    Mag_Star_O1 mag_star_o1 = Mag_Star_O1();

    Mag_Star_O2 mag_star_o2 = Mag_Star_O2();

    Mag_Rot_Star_O2O1 mag_rot_star_o2o1 = Mag_Rot_Star_O2O1();

    void comp(double hC, double B0, double f, double Q=0, int level=5);

    void Jrot_vs_JQB();
    vector<double> errors;
};


#endif //HTMAG_STAR_HPP
