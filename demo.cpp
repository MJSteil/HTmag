/**
 * @file demo.cpp
 * @author M. J. Steil
 * @date 2022.07.26
 * @brief Demo file for Static_Star and Rot_Star
 * @details
 */
#include <iostream>
#include <string>
#include <sys/time.h>

#include "../matplotlib-cpp/matplotlibcpp.h" // https://github.com/MJSteil/matplotlib-cpp

#include "../eos/eos_poly.hpp"
#include "../eos/eos_table.hpp"
#include "../eos/eos_irf.hpp"
#include "../eos/eos_tvii.hpp"

#include "include/Static_Star.hpp"
#include "include/Rot_Star_O1.hpp"

using namespace std;
using namespace units;
using namespace gsl_wrapper;

#pragma clang diagnostic push
#pragma ide diagnostic ignored "IncompatibleTypes"
#pragma ide diagnostic ignored "ConstantConditions"

constexpr bool plotQ = true;
constexpr bool interactiveQ = false;

int main() {
    //region Setup main: start timer, cout.precision(10);
    struct timeval start{}, end{};
    gettimeofday(&start, nullptr);
    cout.precision(10);
    //endregion

    //region EoS setup
    eos_table DD2("../debug/eos/eos_HSDD2.d","cgs");
    DD2.eos::info();
    printf("\n");

    eos_table SFHo("../debug/eos/eos_HSSFHo.d","cgs");
    SFHo.eos::info();
    printf("\n");

    eos_table FSG("../debug/eos/eos_HSFSG.d","cgs");
    FSG.eos::info();
    printf("\n");

    eos_poly Poly(0.05,2.0);
    Poly.eos::info();
    printf("\n");

    eos_irf IRF(2.5069593735181091E-04);
    IRF.eos::info();
    printf("\n");

    eos_tvii TVII(6.2673984337952710E-04,0.1,1.);
    TVII.eos::info();
    printf("\n");

    // Central enthalpies of configurations with specific compactnesses/masses
    //                            0                   M/R=0.01               M/R=0.1                M/R=0.15                  M/R=0.20                   Mmax
    double hC[6][6] = {{2.9546772585658378E-02, 3.3485404922152982E-02, 1.1830856162385392E-01, 1.8178962523458420E-01, 2.7023932076852242E-01, 7.1359566739112523E-01}, // DD2
                       {2.7274388028397274E-02, 3.1732280238905022E-02, 1.2413791302091727E-01, 1.9224688467166381E-01, 2.8799397306439012E-01, 7.2216439845584857E-01}, // SFHo
                       {3.0837169932198309E-02, 3.4390175991622439E-02, 1.2388193026998041E-01, 1.9983187294638846E-01, 3.2100088228750107E-01, 4.8313333633466871E-01}, // FSG
                       {1.0023129052436517E-03, 1.0236881640623766E-02, 1.3150661736031935E-01, 2.3802194607153093E-01, 4.1319532491808797E-01, 4.9255023739938592E-01}, // Poly
                       {3.5559218444199368E-08, 5.0892003156061526E-03, 6.0830199469453478E-02, 1.0271325016717041E-01, 1.5723552817887479E-01, 6.9314718055994529E-01}, // IRF
                       {1.0000000000000000E-04, 8.9272270945435044E-03, 1.0960958925361219E-01, 1.8907843385208778E-01, 2.9815224518724281E-01, 5.4957457263891252E-01}};// TVII

    double hCtab_max[6] = { DD2.hmax, SFHo.hmax, FSG.hmax,3.9, 6.9314718055994529E-01*1.1,5.4957457263891252E-01*1.1};

    if constexpr (plotQ){
        vector<double> hi,PiDD2,rhoiDD2,PiSFHo,rhoiSFHo,PiFSG,rhoiFSG,PiPoly,rhoiPoly;
        double lh = -10;
        while (lh < log(1.2)) {
            const double h = exp(lh);
            hi.emplace_back(h);

            PiDD2.emplace_back(DD2.P(h)/cMeVfm3km2);
            rhoiDD2.emplace_back(DD2.rho(h)/cMeVfm3km2);

            PiSFHo.emplace_back(SFHo.P(h)/cMeVfm3km2);
            rhoiSFHo.emplace_back(SFHo.rho(h)/cMeVfm3km2);

            PiFSG.emplace_back(FSG.P(h)/cMeVfm3km2);
            rhoiFSG.emplace_back(FSG.rho(h)/cMeVfm3km2);

            PiPoly.emplace_back(Poly.P(h)/cMeVfm3km2);
            rhoiPoly.emplace_back(Poly.rho(h)/cMeVfm3km2);
            lh+=0.1;
        }

        matplotlibcpp::named_loglog("DD2",hi,PiDD2);
        matplotlibcpp::named_loglog("SFHo",hi,PiSFHo);
        matplotlibcpp::named_loglog("FSG",hi,PiFSG);
        matplotlibcpp::named_loglog("Poly(0.05|2.0)",hi,PiPoly);
        matplotlibcpp::legend();
        matplotlibcpp::grid(true);
        matplotlibcpp::xlabel("h");
        matplotlibcpp::ylabel("P (MeV fm^-3)");
        matplotlibcpp::save("Poverh.pdf");
        interactiveQ ? matplotlibcpp::show() : void();
        matplotlibcpp::clf();


        matplotlibcpp::named_loglog("DD2",rhoiDD2,PiDD2);
        matplotlibcpp::named_loglog("SFHo",rhoiSFHo,PiSFHo);
        matplotlibcpp::named_loglog("FSG",rhoiFSG,PiFSG);
        matplotlibcpp::named_loglog("Poly(0.05|2.0)",rhoiPoly,PiPoly);
        matplotlibcpp::legend();
        matplotlibcpp::grid(true);
        matplotlibcpp::xlabel("rho (MeV fm^-3)");
        matplotlibcpp::ylabel("P (MeV fm^-3)");
        matplotlibcpp::save("Poverrho.pdf");
        interactiveQ ? matplotlibcpp::show() : void();
        matplotlibcpp::clf();
    }

    //endregion

    //region Static_Star setup
    Static_Star dd2_star(DD2);
    Static_Star sfho_star(SFHo);
    Static_Star fsg_star(FSG);
    Static_Star poly_star(Poly);
    Static_Star irf_star(IRF);
    Static_Star tvii_star(TVII);

    vector<Static_Star> star = {dd2_star,sfho_star,fsg_star,poly_star,irf_star,tvii_star};
    //endregion

    //region Example computation of a single Static_Star
    printf("\n ************** Example computation of a single Static_Star ************** \n\n");
    dd2_star.comp(hC[0][3],1);
    if(plotQ) {
        vector<double> ri, P, m, nu, lambda, rie, nue, lambdae;
        {
            double r = 0;
            while ( r < dd2_star.R) {
                ri.emplace_back(r);
                P.emplace_back(dd2_star.P(r) / cMeVfm3km2);
                m.emplace_back(dd2_star.m(r) / MSkm);
                nu.emplace_back(dd2_star.expNu(r));
                lambda.emplace_back(dd2_star.expLambda(r));
                r += 0.1;
            }
        }
        {
            const double r = dd2_star.R;
            ri.emplace_back(r);
            P.emplace_back(dd2_star.P(r) / cMeVfm3km2);
            m.emplace_back(dd2_star.m(r) / MSkm);
            nu.emplace_back(dd2_star.expNu(r));
            lambda.emplace_back(dd2_star.expLambda(r));
        }
        {
            double r = dd2_star.R;
            while (r < 2 * dd2_star.R) {
                rie.emplace_back(r);
                nue.emplace_back(dd2_star.expNu(r));
                lambdae.emplace_back(dd2_star.expLambda(r));
                r += 0.1;
            }
        }

        matplotlibcpp::plot(ri, P);
        matplotlibcpp::xlabel("r (km)");
        matplotlibcpp::ylabel("P (MeV fm^-3)");
        matplotlibcpp::grid(true);
        matplotlibcpp::save("dd2_static_star_Pofr.pdf");
        (interactiveQ) ? matplotlibcpp::show() : void();
        matplotlibcpp::clf();

        matplotlibcpp::plot(ri, m);
        matplotlibcpp::xlabel("r (km)");
        matplotlibcpp::ylabel("m (Msol)");
        matplotlibcpp::grid(true);
        matplotlibcpp::save("dd2_static_star_mofr.pdf");
        (interactiveQ) ? matplotlibcpp::show() : void();
        matplotlibcpp::clf();

        matplotlibcpp::named_plot("exp(nu)", ri, nu, "r");
        matplotlibcpp::named_plot("exp(lambda)", ri, lambda, "b");
        matplotlibcpp::plot(rie, nue, "r--");
        matplotlibcpp::plot(rie, lambdae, "b--");
        matplotlibcpp::xlabel("r (km)");
        matplotlibcpp::legend();
        matplotlibcpp::grid(true);
        matplotlibcpp::save("dd2_static_star_metricofr.pdf");
        (interactiveQ) ? matplotlibcpp::show() : void();
        matplotlibcpp::clf();
    }
    //endregion


    //region Example computation of maximum mass Static_Star and full Mass-Radius curve
    printf("\n ************** Example computation of maximum mass Static_Star ************** \n\n");
    dd2_star.target(0,0.0,{0.2*DD2.hmax,0.3*DD2.hmax,0.7*DD2.hmax},100,1E-10,1E-14,2);

    printf("\n ************** Example computation of full Mass-Radius Curve ************** \n\n");
    if (plotQ) {
        const double Mmax = dd2_star.M / MSkm;
        const double MBmax = dd2_star.MB / MSkm;
        const double RMmax = dd2_star.R;
        vector<double> M, MB, R;
        for (double lh = -13; lh < log(1.2); lh += 0.025) {
            const double h = exp(lh);
            dd2_star.comp(h, 0);
            M.emplace_back(dd2_star.M / MSkm);
            dd2_star.comp_B(0);
            MB.emplace_back(dd2_star.MB / MSkm);
            R.emplace_back(dd2_star.R);
        }
        matplotlibcpp::named_semilogx("MG", R, M, "b");
        matplotlibcpp::named_semilogx("MB", R, MB, "r");
        matplotlibcpp::semilogx(vector<double>{RMmax}, vector<double>{Mmax}, "b*");
        matplotlibcpp::semilogx(vector<double>{RMmax}, vector<double>{MBmax}, "r*");
        matplotlibcpp::grid(true);
        matplotlibcpp::ylabel("M (Msol)");
        matplotlibcpp::xlabel("R (km)");
        matplotlibcpp::legend();
        matplotlibcpp::save("dd2_MR.pdf");
        (interactiveQ) ? matplotlibcpp::show() : void();
        matplotlibcpp::clf();
    }
    //endregion

    //region Example computation of a single Static_Star
    printf("\n ************** Example computation of a single Rot_Star_O1 ************** \n\n");
    dd2_star.comp(hC[0][3],1);
    Rot_Star_O1 dd2_rotStar(&dd2_star);
    dd2_rotStar.set_f(0.01*cHzkm);
    dd2_rotStar.comp(1,0);
    dd2_rotStar.info();

    if(plotQ){
        vector<double> ri,omega0,rie,omega0e;
        for (double r = 0; r < dd2_star.R; r += 0.1) {
            ri.emplace_back(r);
            omega0.emplace_back(dd2_rotStar.omega0(r));
        }
        const double r = dd2_star.R;
        ri.emplace_back(r);
        omega0.emplace_back(dd2_rotStar.omega0(r));

        for (double r = dd2_star.R; r < 2 * dd2_star.R; r += 0.1) {
            rie.emplace_back(r);
            omega0e.emplace_back(dd2_rotStar.omega0(r));
        }

        matplotlibcpp::plot(ri, omega0, "r");
        matplotlibcpp::plot(rie, omega0e, "r--");
        matplotlibcpp::xlabel("r (km)");
        matplotlibcpp::ylabel("omega/Omega");
        matplotlibcpp::grid(true);
        matplotlibcpp::save("dd2_rotStar_metricofr.pdf");
        (interactiveQ) ? matplotlibcpp::show() : void();
        matplotlibcpp::clf();
    }
    //endregion

    //region End main: stop timer, return 0;
    gettimeofday(&end, nullptr);
    double delta = static_cast<double>((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
    printf("\n\nDone in %.4E s\n",delta);

    return 0;
    //endregion
}
#pragma clang diagnostic pop