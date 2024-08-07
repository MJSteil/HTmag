//
// Global quantities
// Created by M. J. Steil on 2017.02.27.
//

#include "../../include/Rot_Star_O1.hpp"

double Rot_Star_O1::J() {
    return Omega*J0;
}

double Rot_Star_O1::J_error(int print, double err_rel, double err_abs) {
    return abs(1-comp_I_intT(print,err_rel,err_abs)/J0);
}