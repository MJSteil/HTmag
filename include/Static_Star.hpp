//
// Created by M. J. Steil on 2017.02.08.
//

#ifndef HTMAG_STATIC_STAR_HPP
#define HTMAG_STATIC_STAR_HPP

#include <iostream>
#include <vector>

#include "../../units/units.hpp" // https://github.com/MJSteil/units
#include "../../gsl_wrapper/gsl_wrapper.hpp" // https://github.com/MJSteil/gsl_wrapper
#include "../../eos/eos.hpp" // https://github.com/MJSteil/eos

using namespace std;
using namespace units;
using namespace gsl_wrapper;

typedef double (*fptr)(double , void *);

class Static_Star {
public:
    // Constructors
    Static_Star();
    Static_Star(eos &eos_in, string name = "");
    Static_Star(eos &eos_in, double hC_in, int print = 1, vector<double> comp_params = {});

    // Inital conditions
    eos *star_eos;
    string name="";
    double hC;

    // Derived central (r=0) Quantities
    double PC, rhoC, nbarC, drhodhC, dnbardhC, csC;

    // Global parameters
    double M, R, Z, B, MB, EB, nuO;
    int static_star_status = -1; // -1:no EoS, 0:EoS, 1:comp_TOV successful, 2:Data, 3:Interpolations

    // Computational parameters
    // Raw integration gradients
    unsigned int n, count, fails;
    double r0;
    vector<double> a_i, rr_i, drrda_i, z_i, dzda_i;
    gsl_hermite_spline_obj rr_of_a,z_of_a;

    // Derived gradients
    vector<double> r_i, h_i, dhdr_i, dzdr_i;
    gsl_hermite_spline_obj h_of_r;
    gsl_hermite_spline_obj z_of_r;

    // TOV equation in a=h-hC
    double tov_series(double h, int i, int order = 1);
    void comp_tov (double hC_in, int mode = 1, int print = 1, vector<double> comp_params = {});
    void comp_tov_print_params();
    vector<double> comp_tov_params = {1E-16, 1E-3, 1.E-15, 1.E-15, 1., 1., 1, 1.E-3, 5E4};
    const gsl_odeiv2_step_type *comp_tov_step_type = gsl_odeiv2_step_rk8pd;

    // TOV eq. in different variables
    void comp_tov_inr (double hC_in);
    void comp_tov_inP (double hC_in);
    void comp_tov_inh_rm (double hC_in);

    void comp_tov_int();
    void comp(double hC_in, int print = 1, vector<double> comp_params = {});

    // Search for configurations with specific parameters
    double target_fkt(double hC, double target, int type);

    /** Compute specific background star.
     * @param type: 0=Mmax, 1=M_target, 2=R_target, 3=M/R_target, 4=MB_target, 5=Mmin, 6=Rmin, 7=Rmax
     * @param target
     * @param range
     * @param nMax
     * @param abs_error_threshold
     * @param rel_error_threshold
     * @param print
     */
    void target(int type, double target, vector<double> range, int nMax,double abs_error_threshold, double rel_error_threshold, int print);

    /** Compute specific background star using void target(...).
     * @param type: 0=Mmax, 1=M_target, 2=R_target, 3=M/R_target, 4=MB_target, 5=Mmin, 6=Rmin, 7=Rmax
     * @param target
     * @param range
     * @param nMax
     * @param abs_error_threshold
     * @param rel_error_threshold
     * @param print
     * @return target_hC
     */
    double target_hC(int type, double target_int, vector<double> range, int nMax,double abs_error_threshold, double rel_error_threshold, int print);

    // Integrals
    void comp_B(int print = 1);
    vector<double> comp_B_params = {1E-15,1E-10,5E4};

    void comp_MGm1(int print = 1);
    double comp_MGm2(int print = 1, double err_rel=1E-16, double err_abs=1E-20);
    double M_error(int print = 0, double err_rel=1E-16, double err_abs=1E-20);

    void comp_MR(double hC_min, double hC_max,double, double, int &type);

    // Thermodynamic properties
    double h(double r);
    double dh(double r);

    double P(double r);
    double dP(double r);

    double rho(double r);
    double drhodh(double r);

    double nbar(double r);
    double dnbardh(double r);

    double cs(double r);

    // Metric potentials
    double expLambda(double r);
    double dLambda(double r);
    double expNu(double r);
    double dNu(double r);
    double m(double r);
    double dm(double r);

    double j(double r);

    // Auxiliary
    void info();
    void save(string filename);
    void test(double *array);

/*    void comp_tov2 (double hC_in, int mode = 1, int print = 1);
    */
};


#endif //HTMAG_STATIC_STAR_HPP