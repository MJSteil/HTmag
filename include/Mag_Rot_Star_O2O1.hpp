//
// Created by M. J. Steil on 2017.04.05.
//

#ifndef HTMAG_MAG_ROT_STAR_O2O1_HPP
#define HTMAG_MAG_ROT_STAR_O2O1_HPP

#include "Mag_Star_O2.hpp"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>

#include "../lisk/LiSK/lisk.hpp" // https://github.com/MJSteil/LiSK forked from https://github.com/DSNJTurtle/LiSK
#include <gsl/gsl_sf_zeta.h>

class Mag_Rot_Star_O2O1 {
public:
    Mag_Rot_Star_O2O1();
    Mag_Rot_Star_O2O1 (Mag_Star_O2 *);

    Static_Star * star;
    Rot_Star_O1 * rot_star;
    Mag_Star_O1 * mag_star1;
    Mag_Star_O2 * mag_star2;

    double S0(double r);
    double S2(double r);

    // l=1 computations
    int mag_rot_star_O2O1_l1_status = -1;
    unsigned int nl1;
    double deltaJ0, cHl1;

    void comp_l10(int mode = 1, int print = 1, vector<double> comp_params = {});
    double S1(double r);

    vector<double> comp_l1_params = {1e-8, 1E-3, 1E-16, 1E-16, 1.0, 1.0, 1E+4, 1E+1,1e+10, 1e+10,5E3,1E-4};
    void comp_l10_print_params();

    vector<double> rl1_i, u10_i, du10_i, v10_i, dv10_i, u1H0_i, du1H0_i, v1H0_i, dv1H0_i;
    gsl_hermite_spline_obj u10_of_r, v10_of_r, u1H0_of_r, v1H0_of_r;
    void comp_l10_int();

    void comp_deltaI_int(int print = 1);
    double deltaIext();


    // l=3 computations
    int mag_rot_star_O2O1_l3_status = -1;
    unsigned int nl3;
    double W3asymp0, cHl3;

    void comp_l30(int mode = 1, int print = 1, vector<double> comp_params = {});
    double S3(double r);

    vector<double> comp_l3_params = {1e-8, 1E-3, 1E-16, 1E-16, 1.0, 1.0, 1E+4, 1E+1,1e+10, 1e+10,5E3,1E-4};
    void comp_l30_print_params();

    vector<double> rl3_i, u30_i, du30_i, v30_i, dv30_i, u3H0_i, du3H0_i, v3H0_i, dv3H0_i;
    gsl_hermite_spline_obj u30_of_r, v30_of_r, u3H0_of_r, v3H0_of_r;
    void comp_l30_int();

    // Metric perturbations
    double W10_ext(double r, double deltaJ0in, double part);
    double dW10_ext(double r, double deltaJ0in, double part);
    double W10(double r);
    double dW10(double r);

    double W30_ext(double r, double W3asymp0in, double part);
    double dW30_ext(double r, double W3asymp0in, double part);
    double W30(double r);
    double dW30(double r);

    // Global parameters
    double deltaJ();
    double W3asymp();
    void info();

private:
    double zeta3 = gsl_sf_zeta(3.);
    LiSK::LiSK<complex<double>> *lisk_workspace;
};


#endif //HTMAG_MAG_ROT_STAR_O2O1_HPP
