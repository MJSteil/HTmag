//
// Electro-magnetic fields of Mag_Star_O1
// Created by M. J. Steil on 2017.02.21.
//

#include "../../include/Mag_Star_O1.hpp"
// EM planar harmonic components

double Mag_Star_O1::aphi0_ext(double r) {
    double f;
    double M = star->M;
    double y = 1.-2.*M/r;
    double logy = log(y);

    f = (mu0*(2.25 + 1.5*logy + (-3. + 0.75*y)*y))/(M*(1. + (-2. + y)*y)); //aphi0_e(y)

    return f;
}

double Mag_Star_O1::daphi0_ext(double r) {
    double f;
    double M = star->M;
    double y = 1.-2.*M/r;
    double logy = log(y);

    f = (mu0*(-0.75 + (-1.5*logy + 0.75*y)*y))/(M*M*(-1. + y)*y); // daphi0_e/dr(y)


    return f;
}

double Mag_Star_O1::ddaphi0_ext(double r) {
    double f;
    double M = star->M;
    double y = 1.-2.*M/r;
    double logy = log(y);

    f = (mu0*(-0.375 + y*(1.5 + (-1.125 + 0.75*logy)*y)))/(M*M*M*y*y);// ddaphi0_e(y)

    return f;
}


double Mag_Star_O1::aphi0(double r) {
    if(magnetic_l1){
        if(r<=star->R){
            if(mag_star_O1_status<2) {
                GSL_ERROR_VAL ("aphi0(r): No interpolations supplied: run comp_aphi_int().",GSL_FAILURE,GSL_FAILURE);
            }
            return aphi0_of_r.f(r);
        } else {
            return aphi0_ext(r);
        };
    }else{
        return 0;
    }
}

double Mag_Star_O1::daphi0(double r) {
    if(magnetic_l1){
        if(r<=star->R){
            if(mag_star_O1_status<2) {
                GSL_ERROR_VAL ("daphi0(r): No interpolations supplied: run comp_aphi_int().",GSL_FAILURE,GSL_FAILURE);
            }
            return aphi0_of_r.df(r);
        } else {
            return magnetic_l1*daphi0_ext(r);
        };
    }else{
        return 0;    }

}

double Mag_Star_O1::jphi0(double r) {
    if(magnetic_l1){
        if(r<=star->R){
            if(mag_star_O1_status<2) {
                GSL_ERROR_VAL ("jphi0(r): No interpolations supplied: run comp_aphi_int().",GSL_FAILURE,GSL_FAILURE);
            }
            return cjphi0*r*r*(star->P(r)+star->rho(r));
        } else {
            return 0;
        };
    }else{
        return 0;
    }
}



double Mag_Star_O1::at00_ext(double r) {
    if(electric_l2){
        double f;
        double M = star->M;
        double M4 = gsl_pow_4(M);
        double M5 = M*M4;

        double y = 1.-2.*M/r;
        double logy = log(y);

        f = (J0()*mu0*(0.625 + 0.25*logy + (-0.5 + 0.5*logy - 0.125*y)*y))/M4;

        return f;
    } else{
        return 0;
    };
}

double Mag_Star_O1::dat00_ext(double r) {
    if(electric_l2){
        double f;
        double M = star->M;
        double M4 = gsl_pow_4(M);
        double M5 = M*M4;

        double y = 1.-2.*M/r;
        double logy = log(y);

        f = J0()*mu0*(0.125 + y*(-0.25 + 0.25*logy + y*(-0.5*logy + (0.25 + 0.25*logy - 0.125*y)*y)))/(M5*y);

        return f;
    } else{
        return 0;
    };
}

double Mag_Star_O1::at20_ext(double r) {
    if(electric_l2){
        double f;
        double M = star->M;
        double M3 = gsl_pow_3(M);
        double M4 = M*M3;

        double y = 1.-2.*M/r;
        double logy = log(y);

        f = (QEr00*(2.5 + y*(22.5 + 15.*logy + (-22.5 + 15.*logy - 2.5*y)*y)))/(M3*(1. + (-2. + y)*y))
            + (J0()*mu0*(-1.5 - 0.25*logy + y*(-10.25 - 6.75*logy + y*(9.75 - 8.25*logy + (2.25 + 0.25*logy - 0.25*y)*y)))
              )/(M4*(1. + (-2. + y)*y));

        return f;
    } else {
        return 0;
    };
}

double Mag_Star_O1::dat20_ext(double r) {
    if(electric_l2){
        double f;
        double M = star->M;
        double M3 = gsl_pow_3(M);
        double M4 = M*M3;
        double M5 = M*M4;

        double y = 1.-2.*M/r;
        double logy = log(y);

        f = (QEr00*(-21.25 - 7.5*logy + y*(11.25 - 22.5*logy + (11.25 - 1.25*y)*y)))/(M4*(-1. + y))
            + (J0()*mu0*(0.125 + y*(9.875 + 3.625*logy + y*(-3.875 + 11.625*logy + y*(-7.625 - 0.375*logy
            + (1.75 + 0.125*logy - 0.25*y)*y)))))/(M5*(-1. + y)*y);

        return f;
    } else {
        return 0;
    };
}


double Mag_Star_O1::aQ0(double r,int ie) {
    if(r<=star->R){
        return 0;
    } else if(r>star->R||ie){
        return -electric_l0*Qe0/r;
    } else {
        return 0;
    };
}

double Mag_Star_O1::daQ0(double r,int ie) {
    if(r<=star->R&&ie-1){
        return 0;
    } else if(r>star->R||ie){
        return electric_l0*Qe0/r/r;
    } else {
        return 0;
    };
}

double Mag_Star_O1::at00(double r) {
    if(r<=star->R&&electric_l2){
        return 2./3.*aphi0(r);
    } else {
        return electric_l2*at00_ext(r);
    }
}

double Mag_Star_O1::Cat00(double r) {
    return Catoi*(r<=star->R);
}

double Mag_Star_O1::dat00(double r) {
    if(r<=star->R&&electric_l2){
        return 2./3.*daphi0(r);
    } else {
        return electric_l2*dat00_ext(r);
    }
}


double Mag_Star_O1::at20(double r) {
    if(r<=star->R&&electric_l2){
        return -2./3.*aphi0(r);
    } else {
        return electric_l2*at20_ext(r);
    }
}

double Mag_Star_O1::dat20(double r) {
    if(r<=star->R&&electric_l2){
        return -2./3.*daphi0(r);
    } else {
        return electric_l2*dat20_ext(r);
    }
}


double Mag_Star_O1::jt00(double r) {
    if(r<=star->R&&electric_l2){
        double omegabar0 = rot_star->omegabar0(r);
        double domegabar0 = rot_star->domegabar0(r);
        double expLambda = star->expLambda(r);
        double dNu = star->dNu(r);

        return 2./3.*jphi0(r) + ( - aphi0(r)*omegabar0/r + daphi0(r)/expLambda*( omegabar0*(-2+r*dNu) - r*domegabar0 )/2. )/(3.*M_PI*r);
    } else {
        return 0;
    };
}

double Mag_Star_O1::jt20(double r) {
    if(r<=star->R&&electric_l2){
        double omegabar0 = rot_star->omegabar0(r);
        double domegabar0 = rot_star->domegabar0(r);
        double expLambda = star->expLambda(r);
        double dNu = star->dNu(r);

        return -2./3.*jphi0(r) + ( - 2.*aphi0(r)*omegabar0/r + daphi0(r)/expLambda*( omegabar0*(2-r*dNu) + r*domegabar0 )/2. )/(3.*M_PI*r);
    } else {
        return 0;
    };
}

// EM tetrad components

double Mag_Star_O1::Atet(double r,double theta,int i,int j) {
    switch (i){
        case 0: //Atet_t
                return sqrt(star->expNu(r))*(
                                            (at00(r)+at20(r)*(1-1.5*gsl_pow_2(sin(theta))))*Omega()+// l_2
                                            aQ0(r)+ // l_0
                                            Catoi*(r<star->R)
                )*eA;

        case 1: //Atet_r
            return 0;
        case 2: //Atet_theta
            return 0;
        case 3: //Atet_phi
            if(j==0){
                if(r>0){
                    return -eA*aphi0(r)*sin(theta)/r;
                }else{
                    return 0;
                }
            }else{
                if(r>0){
                    return -eA*daphi0(r)*sin(theta)/r;
                }else{
                    return -eA*(-1.)*sin(theta);
                }
            }
        default:
            GSL_ERROR_VAL ("Atet(r,theta,i,j): Invalid input parameters. Run with i=0: A_t, i=1: A_r, i=2: A_theta, i=3: A_phi; j=0 for Aphi, j!=0 dAphi.",GSL_FAILURE,GSL_FAILURE);
    }
}



double Mag_Star_O1::Jtet(double r,double theta,int i) {
    switch (i){
        case 0: //Jtet_t
            if(r>0){
                return B0*Omega()/sqrt(star->expNu(r))*(jt00(r)+jt20(r)*(1-1.5*gsl_pow_2(sin(theta))));
            }else{

                return 0;//B0*Omega()/sqrt(star->expNu(1E-10))*(jt00(1E-10)+jt20(1E-10)*(1-1.5*gsl_pow_2(sin(theta))));
            }
        case 1: //Jtet_r
            return 0;
        case 2: //Jtet_theta
            return 0;
        case 3: //Jtet_phi
            if(r>0){
                return -B0*jphi0(r)*sin(theta)/r;
            }else{
                return 0;
            }
        default:
            GSL_ERROR_VAL ("Jtet(r,theta,i): Invalid input parameters. Run with i=1: J_r, i=2: J_theta, i=3: J_phi.",GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O1::Btet(double r,double theta,int i) {
    switch (i){
        case 1: //Btet_r
            if(r>0){
                return -2*B0*aphi0(r)/gsl_pow_2(r)*cos(theta);
            }else{
                return B0*cos(theta);
            }
        case 2: //Btet_theta
            if(r>0){
                return B0/sqrt(star->expLambda(r))/r*daphi0(r)*sin(theta);
            }else{
                return -B0*sin(theta);
            }

        case 3: //Btet_phi
            return 0;
        default:
            GSL_ERROR_VAL ("Btet(r,theta,i): Invalid input parameters i=1: B_r, i=2: B_theta, i=3: B_phi.",GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O1::Etet(double r,double theta,int i, int ie) {
    switch (i) {
        case 1: //Etet_r
            if (r <= star->R&&ie-1&&r!=0&&electric_l2) {
                return rot_star->omegabar0(r)*daphi0(r)*gsl_pow_2(sin(theta))/sqrt(star->expNu(r)*star->expLambda(r))*Omega()*B0;
            } else if(r>star->R||ie) {
                double Xtheta = gsl_pow_2(sin(theta));
                if(electric_l2){
                    return (daQ0(r,ie) + ( dat00_ext(r) + dat20_ext(r)*(1-1.5*Xtheta) - rot_star->omega0(r)*daphi0(r)*Xtheta )*Omega())*eA;
                           /*/sqrt(star->expNu(r)*star->expLambda(r)) //=1 for r>R */;
                }else{
                    return (daQ0(r,ie))*eA/*/sqrt(star->expNu(r)*star->expLambda(r)) //=1 for r>R */;
                }

            } else {
                return 0;
            }
        case 2: //Etet_theta
            if (r <= star->R&&r!=0&&electric_l2) {
                return rot_star->omegabar0(r)*aphi0(r)/r/sqrt(star->expNu(r))*2.*sin(theta)*cos(theta)*Omega()*B0;
            } else if(r>star->R&&electric_l2) {
                return -( 3.*at20_ext(r) + 2.*aphi0(r)*rot_star->omega0(r) )/r/sqrt(star->expNu(r))*sin(theta)*cos(theta)*Omega()*B0;
            } else {
                return 0;
            }

        case 3: //Etet_phi
            return 0;
        default:
            GSL_ERROR_VAL ("Etet(r,theta,i): Invalid input parameters i=1: E_r, i=2: E_theta, i=3: E_phi.",
                           GSL_FAILURE, GSL_FAILURE);
    }

}

double Mag_Star_O1::L1(double r, double theta) {
    return Btet(r,theta,1)*Btet(r,theta,1) + Btet(r,theta,2)*Btet(r,theta,2)
            - ( Etet(r,theta,1)*Etet(r,theta,1) + Etet(r,theta,2)*Etet(r,theta,2) )
           + r*sin(theta)*rot_star->omega0(r)*rot_star->Omega/sqrt(star->expNu(r))*
            ( Btet(r,theta,1)*Etet(r,theta,2) - Btet(r,theta,2)*Etet(r,theta,1) );
}

double Mag_Star_O1::L2(double r, double theta) {
    return 2.*( Btet(r,theta,1)*Etet(r,theta,1) + Btet(r,theta,2)*Btet(r,theta,2) );
}

double Mag_Star_O1::Atet_2D_plot(double r,int i,int l) {
    switch (i){
        case 0: //Atet_t
            if(l==0){
                return Atet(r,asin(sqrt(2./3.)),0)/cTkm*1E3;
            }else{
                return (Atet(r,0,0)-Atet(r,asin(sqrt(2./3.)),0))/cTkm*1E3;
            }

        case 1: //Atet_r
            return 0;
        case 2: //Atet_theta
            return 0;
        case 3: //Atet_phi/sin(theta)
            return Atet(r,M_PI/2,3,l)/cTkm*1E3;
        default:
            GSL_ERROR_VAL ("Atet_2D_plot(r,i,j): Invalid input parameters. Run with i=0: A_t, i=1: A_r, i=2: A_theta, i=3: A_phi; j=0 for Aphi, j!=0 dAphi.",GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O1::Jtet_2D_plot(double r, int i, int l) {

    switch (i){
        case 0: //Jtet_t
            if(l==0){ //l=0
                return Jtet(r,asin(sqrt(2./3.)),0)/cA*1e-6;
            }else{  // l=2
                return ( Jtet(r,0,0) - Jtet(r,asin(sqrt(2./3.)),0) )/cA*1e-6;
            }
        case 3: //Jtet_phi/sin(theta)
            if(r>0){
                return Jtet(r,M_PI/2.,3)/cA*1e-6;
            }else{
                return 0;
            }
        default:
            GSL_ERROR_VAL ("Jtet_2D_plot(r,i): Invalid input parameters. Run with i=0: J_t_l, i=3: J_phi.",GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O1::Btet_2D_plot(double r,int i) {
    switch (i){
        case 1: //Btet_r/cos(theta);
            return Btet(r,0,1)/cTkm;
        case 2: //Btet_theta/sin(theta)
            return Btet(r,M_PI/2.,2)/cTkm;
        default:
            GSL_ERROR_VAL ("Btet_2D_plot(r,i): Invalid input parameters i=1: B_r, i=2: B_theta.",GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O1::Etet_2D_plot(double r,int i,int l) {

    if(1) {
        switch (i) {
            case 1: //Etet_r/gsl_pow_2(sin(theta)) or /(1-1.5*gsl_pow_2(sin(theta))) for r>R [V/m]
                if(r<= star->R){
                    return Etet(r,M_PI/2.,1)/(cV)*1E-3;
                }else{
                    if(l==0){
                        return Etet(r,asin(sqrt(2./3.)),1)/(cV)*1E-3;
                    }else{
                        return (Etet(r,0,1)-Etet(r,asin(sqrt(2./3.)),1))/(cV)*1E-3;;
                    }
                }

            case 2: //Etet_theta/(sin(theta)*cos(theta)) [V/m]
                return 2.*Etet(r,M_PI/4.,2)/(cV)*1E-3;
            default:
                GSL_ERROR_VAL ("Etet_2D_plot(r,i): Invalid input parameters i=1: E_r_j, i=2: E_theta.",
                               GSL_FAILURE, GSL_FAILURE);
        }
    }else{
        return 0;
    }
}

double Mag_Star_O1::sigma_e(double theta){
    return (Etet(star->R,theta,1, 1)-Etet(star->R,theta,1, 0))/M_4PI;
}

double Mag_Star_O1::dQei0(double r){
    if(r>0){
        return -M_4PI*r*r*sqrt(star->expLambda(r)/star->expNu(r))*(jt00(r)-2./3.*rot_star->omega0(r)*jphi0(r));
    }else{
        return 0;
    }

}

