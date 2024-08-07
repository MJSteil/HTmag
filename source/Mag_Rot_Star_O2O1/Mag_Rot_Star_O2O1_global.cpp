//
// Global parameters of Mag_Rot_Star_O2O1
// Created by M. J. Steil on 2017.04.05.
//

#include "../../include/Mag_Rot_Star_O2O1.hpp"

double Mag_Rot_Star_O2O1::deltaJ() {
    return deltaJ0*gsl_pow_2(mag_star1->B0)*rot_star->Omega;
}

double Mag_Rot_Star_O2O1::W3asymp() {
    return W3asymp0*gsl_pow_2(mag_star1->B0)*rot_star->Omega;
}
