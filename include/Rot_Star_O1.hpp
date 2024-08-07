//
// Created by M. J. Steil on 2016.07.06.
//

#ifndef NEUTRONSTAR_ROT_STAR_HPP
#define NEUTRONSTAR_ROT_STAR_HPP

#include "Static_Star.hpp"


using namespace std;
using namespace gsl_wrapper;
using namespace units;

class Rot_Star_O1 {
public:
    Rot_Star_O1();
    Rot_Star_O1(Static_Star *,string name = "");
    Rot_Star_O1(Static_Star *,double f /*Hz*/, int print = 1);

    int rot_star_status = -1; // -1:no Background star, 0: Valid background star, 1: Valid omegabar0 computation, 2: Interpolations generated

    Static_Star * star;

    double Omega;
    double f;
    double J0;
    double I;
    double omegabar0C;

    unsigned int n, steps;
    vector<double> r_i, omegabar0_i, domegabar0_i, ddomegabar0_i;
    gsl_hermite_spline_obj omegabar0_of_r, domegabar0_of_r;

    void comp_omegabar0(int mode = 0, int print = 1, vector<double> comp_params = {});
    vector<double> comp_omegabar0_params = {1e-15, 1E-3, 1E-15, 1E-12, 1.0, 1.0, 1e-30, 1e-30, 1E5};
    void comp_omegabar0_print_params();
    void comp_omegabar0_int();

    void set_f(double f0);
    void comp(int mode = 0, int print = 1, vector<double> comp_params = {});
    // Global parameters
    double J();

    // Metric potentials
    double omegabar0_ext(double r, int type);

    double omegabar0(double r);
    double domegabar0(double r);
    double omega0(double r);
    double domega0(double r);

    // Auxiliary
    void info();

    void comp_I_int(int print = 1, double err_rel=1E-14, double err_abs=1E-10);
    double comp_I_intT(int print = 1, double err_rel=1E-14, double err_abs=1E-10);
    double J_error(int print = 1, double err_rel=1E-14, double err_abs=1E-10);
};


#endif //NEUTRONSTAR_ROT_STAR_HPP
