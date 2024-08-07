//
// Created by M. J. Steil on 2016.07.04.
//

#ifndef NEUTRONSTAR_MAG_STAR_O1_HPP
#define NEUTRONSTAR_MAG_STAR_O1_HPP

#include "Static_Star.hpp"
#include "Rot_Star_O1.hpp"

using namespace std;
using namespace units;
using namespace gsl_wrapper;

class Mag_Star_O1 {
public:
    Mag_Star_O1();
    Mag_Star_O1(Static_Star *);
    Mag_Star_O1(Static_Star *, Rot_Star_O1 *);
    Mag_Star_O1(Static_Star *, double B0, double Qe=0, int mode = 1, int print = 1, vector<double> comp_cjphi0_range_in = {},vector<double> comp_params = {});

    int magnetic_l1 = 0;
    int electric_l0 = 0;
    int electric_l2 = 0;

    Static_Star * star;

    Rot_Star_O1 * rot_star;
    double Omega();
    double J0();

    double B0  = 0.;    // Magnetic factor Btet(0,0,1)=Br00*eA
    double Qe0 = 0.;    // Electric factor Qe = Qe0*eA


    /**
     * Computes B0 and Qe0 to realise a Mag_Star with Btet(0,0,1)=Br00 and Qe = Qe.
     * @param Br00 [km**-1] Central magnetic field e.g. 500*cGTkm;
     * @param Qe [1] Global electric l=0 Charge e.g. 1E10*cC;
     */
    void set_B0Q0(double Br00 /*[km**-1]*/, double Qe = 0/*[1]*/);

    double eA = 1.;  // Order A scale factor

    // Magnetic
    int mag_star_O1_status = -1; // -1:no Background star, 0: Valid background star, 1: Valid aphi0_full computation, 2: Interpolations generated

    double mu0 = 1;
    double cjphi0 = 0;
    double DeltadaphiRs;

    unsigned int n;
    vector<double> r_i, aphi0_i, daphi0_i, ddaphi0_i;
    gsl_hermite_spline_obj aphi0_of_r, daphi0_of_r;

    // Computations

    void comp_aphi0(double comp_aphi_cjphi0, int mode = 0, int print = 1,vector<double> comp_params = {});
    void comp_cjphi0(int print = 1, vector<double> comp_cjphi0_range_in = {},vector<double> comp_params = {});
    vector<double> comp_cjphi0_range = {-10,0,0};
    vector<double> comp_aphi0_params = {1e-15, 1E-3, 1E-16, 1E-12, 1.0, 1.0, 1e-30, 1e-30, 1E5};
    void comp_aphi0_print_params();

    void comp_aphi0_int();

    void comp(double B0, int mode = 1, int print = 1, vector<double> comp_params = {}, vector<double> comp_cjphi0_range_in = {});


    // Monopole Electric

    // Quadrupole Electric
    int rot_mag_star = 0;

    // Global parameters

    double mu();
    double cjphi();
    double QEr();

    double Qe();
    double QeS();
    double QeI();
    double Qe_error(int print = 1 , double err_rel=1E-15, double err_abs=1E-15);

    void set_rot(Rot_Star_O1 * rot_star_in);

    // EM zonal harmonic components

    double aphi0_ext(double r);
    double daphi0_ext(double r);
    double ddaphi0_ext(double r);

    double aphi0(double r);
    double daphi0(double r);
    double jphi0(double r);

    double Catoi = 0;
    double Cat00(double r);
    double QEr00 =0;

    double at00_ext(double r);
    double dat00_ext(double r);

    double aQ0(double r,int ie=0);
    double daQ0(double r,int ie=0);

    double at20_ext(double r);
    double dat20_ext(double r);

    double at00(double r);
    double dat00(double r);
    double at20(double r);
    double dat20(double r);

    double jt00(double r);
    double jt20(double r);

    // EM tetrad components

    /*!
     * Tetrad components of the four potential Atet={A_t,A_r,A_theta,A_phi}_tet
     * @param r radius [km]
     * @param theta polar angle [rad]
     * @param i i-th component: i=1: _r, i=2: _theta, i=3: _phi
     * @param j j-th derivative: j=0: no derivative, j=1: first derivative
     * @return [d^j/dq^j A_i]tet
     */
    double Atet(double r, double theta, int i = 3, int j = 0);

    double Jtet(double r, double theta, int i = 3);             // [J_i]tet
    double Btet(double r, double theta, int i);
    double Etet(double r, double theta, int i, int ie=0);

    double sigma_e(double theta);
    double dQei0(double r);

    double L1(double r, double theta);
    double L2(double r, double theta);

    // Auxiliary
    void info();
    void info_E();
    double s_B0(int print=1);


    double Jtet_2D_plot(double r, int i, int l);
    double Btet_2D_plot(double r, int i);
    double Etet_2D_plot(double r,int i, int l);
    double Atet_2D_plot(double r,int i, int l);

    double comp_QeI(int print =1, double err_rel=1E-14, double err_abs=1E-10);
    void comp_QeS(int print =1);

};

#endif //NEUTRONSTAR_MAG_STAR_O1_HPP
