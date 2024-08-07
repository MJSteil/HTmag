//
// Global parameters of Mag_Star_O1
// Created by M. J. Steil on 2017.02.24.
//

#include "../../include/Mag_Star_O1.hpp"

double Mag_Star_O1::mu() {
    return eA*mu0;
}

double Mag_Star_O1::cjphi() {
    return eA*cjphi0;
}

double Mag_Star_O1::QEr() {
    return eA*Omega()*QEr00;
}

double Mag_Star_O1::Qe() {
    return eA*Qe0;
}

double Mag_Star_O1::QeS() {

    double R = star->R;
    double M = star->M;

    return (Qe0-(R*R*R-2*J0())*mu0*(2.*M*(M-R)+(2.*M-R)*R*log(1.-2.*M/R))/(2.*M*M*M*R*(2*M-R))*Omega())*eA;
}

double Mag_Star_O1::QeI() {

    return (Qe()-QeS());
}

double Mag_Star_O1::Qe_error(int print,double err_rel, double err_abs) {
    double DeltaQ = abs(1+comp_QeI(print,err_rel,err_abs)/QeS());
    if(print) printf("Delta Q = %.4E\n",DeltaQ);
    return DeltaQ;
}

double Mag_Star_O1::comp_QeI(int print,double err_rel, double err_abs) {

    double R = star->R;
    double M = star->M;

    gsl_function_pp_fkt comp_Qe0_dQe0 = gsl_function_pp_fkt([&](double r)->double{
        //return M_4PI/3.*r*r*sqrt(star->expLambda(r)/star->expNu(r))*(3.*jt00(r)-2.*rot_star->omega0(r)*jphi0(r));
        return dQei0(r);
    });

    // Integration Method
    gsl_integration comp_Qe0_int(comp_Qe0_dQe0.gsl_fkt(),{0,star->R},err_rel,err_abs,"odeiv2",(int)1E6,6,print==2);

    if(print){

        printf("|=>QeI: %.12E C (QeI+QeS=%.12E) \n",comp_Qe0_int.F*B0*rot_star->Omega/cCkm,comp_Qe0_int.F*B0*rot_star->Omega/cCkm+QeS()/cCkm);
    }
    return comp_Qe0_int.F*B0*rot_star->Omega;
}

void Mag_Star_O1::comp_QeS(int print) {


    double R = star->R;
    double M = star->M;

    gsl_function_pp_fkt comp_QeS_dQs = gsl_function_pp_fkt([&](double r)->double{
        return 2.*M_PI*sigma_e(r)*sin(r)*gsl_pow_2(star->R);
    });

    // Integration Method
    gsl_integration comp_Qe0_int(comp_QeS_dQs.gsl_fkt(),{0,M_PI},1E-15,1E-6,"qag",print);


    if(print){
        if(electric_l2){
            printf("|=>QsA1: %.12E C\n",(Qe0-electric_l2*2./3*R*R*sqrt(1/star->expLambda(R)/star->expNu(R))*rot_star->omegabar0(R)*daphi0(R)*Omega())*eA/cCkm);
        }
        printf("|=>QsA2: %.12E C\n",(Qe0-(R*R*R-2*J0())*mu0*(2.*M*(M-R)+(2.*M-R)*R*log(1.-2.*M/R))/(2.*M*M*M*R*(2*M-R))*Omega())*eA/cCkm);
        printf("|=>Qs: %.12E C\n",comp_Qe0_int.F/cCkm);
    }

}

double Mag_Star_O1::s_B0(int print){

    gsl_function_pp_fkt comp_s_B0_fkt = gsl_function_pp_fkt([&](double s)->double{
        return -Jtet(s*star->R,M_PI/2.,3);//daphi0(s*star->R);
    });

    //gsl_rootfinder comp_s_B0(comp_s_B0_fkt, 200, 1e-16,0, {0.3,0,1.},print);
    gsl_minimizer comp_s_B0(comp_s_B0_fkt.gsl_fkt(), 200, 1e-16,0, {0.01,.5,.99},gsl_min_fminimizer_brent,print);

    return comp_s_B0.x_min;
}