//
// Created by m. J. Steil on 2016.08.11.
//

#ifndef NEUTRONSTAR_MAG_STAR_O2_HPP
#define NEUTRONSTAR_MAG_STAR_O2_HPP

#include "Static_Star.hpp"
#include "Mag_Star_O1.hpp"

class Mag_Star_O2 {
public:
    Mag_Star_O2();
    Mag_Star_O2(Mag_Star_O1 *); // -1:no Background star, 0: Valid background star

    Static_Star * star;
    Mag_Star_O1 * mag_star;

    // l=0 computations
    int mag_star_O2_l0_status = -1;
    unsigned int nl0;
    double deltaM0, ch00, ch0r0, m0R, h0R, MBR;

    vector<double> rl0_i, m00_i, dm00_i, h00_i, dh00_i, MB0_i, dMB0_i;

    gsl_hermite_spline_obj m00_of_r, h00_of_r;

    void comp_l00(int mode = 1, int print = 1, vector<double> comp_params = {});
    void comp_l00_odeiv(int mode = 1, int print = 1, vector<double> comp_params = {});
    void comp_l00_match(int print = 1);
    void comp_l00_MB0(int print = 1, int mode = 0);

    void comp_l00_print_params();
    vector<double> comp_l0_params = {1e-10, 1E-3, 1E-16, 1E-12, 1.0, 1.0, 1e-30, 1e-30, 1e10, 5E4, 0.};
    int deltaM0_C = 1;

    void comp_l00_int();

    void comp_deltaB0(int print = 1);
    double deltaB0, deltaMB0;

    // l=1 computations
    int mag_star_O2_l1_status = -1;
    unsigned int nl1;
    double deltaJ20_00;

    vector<double> rl1_i, omegaQ00_i, domegaQ00_i, ddomegaQ00_i;

    gsl_hermite_spline_obj omegaQ00_of_r, domegaQ00_of_r;

    void comp_l100(int mode = 1, int print = 1, vector<double> comp_params = {});
    void comp_l100_print_params();
    vector<double> comp_l100_params = {1e-20, 1E-3, 1E-16, 1E-12, 1.0, 1.0, 1e-30, 1e-30, 5E4};

    void comp_l100_int();

    // l=2 computations
    int mag_star_O2_l2_status = -1;
    unsigned int nl2;
    double Qm0, cn2H00;

    vector<double> rl2_i, n20_i, dn20_i, y20_i, dy20_i, n2H0_i, dn2H0_i, y2H0_i, dy2H0_i;
    gsl_hermite_spline_obj n20_of_r, y20_of_r, n2H0_of_r, y2H0_of_r;

    void comp_l20(int mode = 1, int print = 1, vector<double> comp_params = {});
    void comp_l20_print_params();
    vector<double> comp_l2_params = {1e-10, 1E-3, 1E-12, 1E-8, 1.0, 1.0, 1e-1, 1e-1,1e+10, 1e+10,5E4};

    void comp_l20_int();

    // Global parameters
    double BB();

    double deltaM();
    double deltaMQ();

    /**
     * Komar mass shift integral - EM sources
     * @param print
     * @return
     */
    double deltaM_EM(int print = 1);
    /**
     * Komar mass shift integral - fluid sources
     * @param print
     * @return
     */
    double deltaM_F(int print = 1);

     /**
     * Komar mass shift integral
     * @param print
     * @return | 1 - (deltaM_EM()+deltaM_F())/deltaM()|
     */
    double deltaM_error(int print = 1);

    double deltaB(int print = 1);
    double deltaMB(int print = 1);

    double deltaR_l0();
    double deltaA();
    double A_r(double r);

    double Qm();
    double rsph(double r, double theta);
    double e();
    double epsilon();
    double epsilonErr();

    double Bmax();
    double Bmax_hc();
    double Bth(double thin = .05);

    double deltaJ20();
    double deltaJ20_int(int print = 1);
    double deltaJ20_error(int print = 1);

    // Metric perturbations
    double m00_ext(double r, double deltaMin);
    double n00_ext(double r, double deltaMin);

    double omegaQ00_ext(double r, double deltaJ20_00in);
    double domegaQ00_ext(double r, double deltaJ20_00in);

    double n20_ext(double r, double Qm0in, double part = 1.);
    double dn20dr_ext(double r, double Qm0in, double part = 1.);
    double v20_ext(double r, double Qm0in, double part = 1.);
    double y20_ext(double r, double Qm0in, double part = 1.);

    double m00(double r);
    double dm00(double r);
    double n00(double r);
    double dn00(double r);

    double h00(double r);
    double dh00(double r);
    double xi00(double r);

    double omegaQ00(double r);
    double domegaQ00(double r);

    double n20(double r);
    double dn20(double r);
    double y20(double r);
    double v20(double r);
    double dv20(double r);
    double m20(double r);
    double dm20(double r);

    double h20(double r);
    double dh20(double r);
    double xi20(double r);


    double rbar(double r, double theta);
    double h(double r, double theta);

    // Tests

//    double deltaM0_int_error(double fin);
//    void opt_deltaM0();

    void test_isoSurface();
    void test_Bernoulli(double rel_r, double theta);
    void eSummary();

    // Auxiliary
    void info();

};


#endif //NEUTRONSTAR_MAG_STAR_O2_HPP
