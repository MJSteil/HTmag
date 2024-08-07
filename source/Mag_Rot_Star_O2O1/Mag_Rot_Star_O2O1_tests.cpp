//
// Mag_Rot_Star_O2O1 test methods
// Created by M. J. Steil on 2017.04.05.
//

#include "../../include/Mag_Rot_Star_O2O1.hpp"

void Mag_Rot_Star_O2O1::comp_deltaI_int(int print) {


    // gsl_function via template wrapper gsl_function_pp
    auto dI1f = [=](double r)->double{
        double aphi0 = mag_star1->daphi0(r);
        double daphi0 = mag_star1->daphi0(r);

        double omegabar0 = rot_star->omegabar0(r);
        double expLambda = star->expLambda(r);
        double expNu = star->expNu(r);
        double P = star->P(r);
        double rho= star->rho(r);

        return 8./3.*M_PI*sqrt(expLambda/expNu)*gsl_pow_4(r)*(
                -W10(r)*(P+rho)
                +(P+rho+star->drhodh(r))*(mag_star2->h00(r)-1./5.*mag_star2->h20(r))*omegabar0
                +omegabar0/M_PI/5.*(aphi0*aphi0/gsl_pow_4(r)+daphi0*daphi0/r/r/expLambda)
                -0.2/r*omegabar0*(P+rho)*( 5*r*(mag_star2->h00(r)-mag_star2->h20(r))
                                           +expLambda*(mag_star2->m20(r)-5*mag_star2->m00(r))
                                           +4*r*mag_star2->v20(r) )

        );
    };
    gsl_function_pp<decltype(dI1f)> dI1fp(dI1f);
    gsl_function *dI1 = static_cast<gsl_function*>(&dI1fp);

    printf("{M->%.8E,y->%.8E,mu->%.8E,cat2H2->%.8E,cat0H2->%.8E,J->%.8E,QM->%.8E}\n",star->M,1-2*star->M/star->R,mag_star1->mu0,0.,0.,rot_star->J0,mag_star2->Qm0);


    // Integration Method
    gsl_integration comp_I1_int(dI1,{0,star->R},100,1E-16,"qag",print);



    double R = star->R;
    double j = star->j(R);
    double domega = rot_star->domega0(R);
    printf("I = %.8E \n", -1./6./j*gsl_pow_4(R)*domega/gsl_pow_3(MSkm));
    double dIW1= -1./6./j*gsl_pow_4(R)* dW10(R);
    double dIS = +1/30./j*gsl_pow_3(R)*domega*(
                                                      5*R*(mag_star2->h00(R) -mag_star2->h20(R))
                                                      +star->expLambda(R)*(5*mag_star2->m00(R)-mag_star2->m20(R))
                                                      +4*R*mag_star2->v20(R)
    );
    double dIxi = 1./15./j*(star->dNu(R)+star->dLambda(R))*R*R*R*rot_star->omegabar0(R)*(5*mag_star2->xi00(R)-mag_star2->xi20(R));


    printf("deltaI = %.8E => rel_error=%.8E \n", dIW1+dIS+deltaIext()+dIxi, 1-(dIW1+dIS+deltaIext()+dIxi)/deltaJ0);
    printf("deltaIi = %.8E \n", dIW1+dIS+dIxi);
    printf("dIW1= %.4E ,dIS= %.4E ,dIxi= %.4E \n", dIW1,dIS,dIxi);
}

double Mag_Rot_Star_O2O1::deltaIext() {
    double M = star->M;
    double MMM = M*M*M;
    double MMMM = MMM*M;

    double y = 1.-2.*M/star->R;
    double yy = y*y;
    double yyy = yy*y;
    double yyyy = yyy*y;
    double yyyyy = yyyy*y;
    double yyyyyy = yyyyy*y;

    double logy = log(y);
    double logylogy = logy*logy;
    double ym1 = y-1.;
    double ym1ym1 = ym1*ym1;
    double ym1ym1ym1 = ym1ym1*ym1;

    double mu0 = mag_star1->mu0;
    double mu0mu0 = gsl_pow_2(mu0);
    double cat0H2 = 0;
    double cat2H2 = 0;
    double J0 = rot_star->J0;

    double T1 = 0.5*cat0H2*mu0*(3. + 2.*logy - 4.*y + yy)/M/ym1ym1;
    double T2 = -0.2*cat2H2*mu0*(12.*logylogy*(1. + 3.*y)
                                 + ym1ym1*(51. + 7.*y - 11.*yy + yyy)
                                 + 4.*logy*(13. + 3.*y - 21.*yy + 5.*yyy))/gsl_pow_5(ym1);
    double T3 = -0.15*J0*mu0mu0/MMMM/gsl_pow_5(ym1)*(-2.*logylogy*(11. + 39.*y - 3.*yy + yyy)
                                                     + ym1ym1*(-99. - 24.*y + 34.*yy - 8.*yyy + yyyy)
                                                     - 1.*logy*(99. + 45.*y - 202.*yy + 66.*yyy - 9.*yyyy + yyyyy));

    return T1+T2+T3;
}