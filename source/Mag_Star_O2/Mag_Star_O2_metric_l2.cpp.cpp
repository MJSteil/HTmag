//
// Mag_Star_O2 metric perturbations
// Created by M. J. Steil on 2017.02.27.
//

#include "../../include/Mag_Star_O2.hpp"

// Analytic exterior solutions

double Mag_Star_O2::n20_ext(double r, double Qm0in, double part) {
    double M = star->M;
    double MMM = M*M*M;
    double MMMM = MMM*M;

    double y = 1.-2.*M/r;
    double yy = y*y;
    double yyy = yy*y;
    double yyyy = yy*yy;
    double logy = log(y);
    double logylogy = logy*logy;
    double ym1 = y-1.;
    double ym1ym1 = ym1*ym1;

    double mu0mu0 = gsl_pow_2(mag_star->mu0);

    double n20_ext_H = 0.;
    double n20_ext_P = 0.;

    if(Qm0in!=0.){
        n20_ext_H = -Qm0in*5.*( -1. + 8.*(y - yyy) + yyyy + 12.*yy*logy )/(16.*MMM*ym1ym1*y);
    }

    if(part!=0.){
        n20_ext_P = 3.*mu0mu0/(4.*MMMM*ym1*y)*(
                ( -7. + 33.*y + 37.*yy - 3*yyy )/(8.) +
                ( 1. - 24.*yy - 8*yyy + yyyy )*logy/(4.*ym1)+
                ( yy*logylogy )/(ym1)
        );
    }

    return n20_ext_H + n20_ext_P;
}

double Mag_Star_O2::dn20dr_ext(double r, double Qm0in, double part) {
    double M = star->M;
    double MMM = M*M*M;
    double MMMM = MMM*M;

    double y = 1.-2.*M/r;
    double yy = y*y;
    double yyy = yy*y;
    double yyyy = yy*yy;
    double logy = log(y);
    double logylogy = logy*logy;
    double ym1 = y-1.;
    double ym1ym1 = ym1*ym1;

    double mu0mu0 = gsl_pow_2(mag_star->mu0);

    double dn20dr_ext_H = 0.;
    double dn20dr_ext_P = 0.;

    if(Qm0in!=0.){
        dn20dr_ext_H = -Qm0in*5.*( -1. + 3.*y - 28.*(yy-yyy) -3.*yyyy + y*yyyy - 12.*yy*(y+1.)*logy )/(32.*MMMM*ym1*yy);
    }

    if(part!=0.){
        dn20dr_ext_P = 3.*mu0mu0/(8.*MMMM*M)*(
                -( 5. - 14.*y + 118.*yy + 10*yyy + yyyy )/(8.*yy) +
                ( 1. - 3.*y + 16.*yy +48.*yyy -3.*yyyy + y*yyyy )*logy/(4.*ym1*yy)
                -( (1+y)*logylogy )/(ym1)
        );
    }

    return dn20dr_ext_H + dn20dr_ext_P;
}

double Mag_Star_O2::y20_ext(double r, double Qm0in, double part) {
    double M = star->M;
    double MMM = M*M*M;
    double MMMM = MMM*M;

    double y = 1.-2.*M/r;
    double yy = y*y;
    double yyy = yy*y;
    double yyyy = yy*yy;
    double logy = log(y);
    double logylogy = logy*logy;
    double ym1 = y-1.;
    double ym1ym1 = ym1*ym1;

    double mu0mu0 = gsl_pow_2(mag_star->mu0);

    double y20_ext_H = 0.;
    double y20_ext_P = 0.;

    if(Qm0in!=0.){
        y20_ext_H = Qm0in*5.*( -1. -9*y +9.*yy + yyy - 6.*y*(1.+y)*logy )/(16.*MMM*ym1*y);
    }

    if(part!=0.){
        y20_ext_P = 3.*mu0mu0/(16.*MMMM)*(
                ( 6. + 45.*y + 8*yy + yyy)/(2.*y) +
                -( 1. + 12.*y + 14.*yy +3*yyy )*logy/(ym1*y)+
                ( -1. + 3*y )*logylogy/(ym1)
        );
    }

    return y20_ext_H + y20_ext_P;
}

double Mag_Star_O2::v20_ext(double r, double Qm0in, double part) {
    double M = star->M;
    double MMM = M*M*M;
    double MMMM = MMM*M;

    double y = 1.-2.*M/r;
    double yy = y*y;
    double yyy = yy*y;
    double yyyy = yy*yy;
    double logy = log(y);
    double logylogy = logy*logy;
    double ym1 = y-1.;
    double ym1ym1 = ym1*ym1;

    double mu0mu0 = gsl_pow_2(mag_star->mu0);

    double v20_ext_H = 0.;
    double v20_ext_P = 0.;

    if(Qm0in!=0.){
        v20_ext_H = Qm0in*5.*( -1. -9*y +9.*yy + yyy - 6.*y*(1.+y)*logy )/(16.*MMM*ym1*y);
    }

    if(part!=0.){
        v20_ext_P = 3.*mu0mu0/(16.*MMMM)*(
                ( 7. + 62.*y + 7.*yy )/(2.*y) +
                -( 1. + 22.*y + 22.*yy + yyy )*logy/(ym1*y)+
                ( 3. + 2.*y + 3*yy )*logylogy/(ym1ym1)
        );
    }

    return v20_ext_H + v20_ext_P;
}


// Complete metric potentials

double Mag_Star_O2::n20(double r){
    if(mag_star_O2_l2_status>=2){
        if(r<=star->R&&r>0){
            return n20_of_r.f(r);
        }else if(r>star->R){
            return n20_ext(r,Qm0,1);
        }else{
            return 0;
        }
    }else{
        GSL_ERROR_VAL ("n20: No interpolations supplied: run comp_l20_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::dn20(double r){
    if(mag_star_O2_l2_status>=2){
        if(r<=star->R&&r>0){
            return n20_of_r.df(r);
        }else if(r>star->R){
            return dn20dr_ext(r,Qm0,1);
        }else{
            return 0;
        }
    }else{
        GSL_ERROR_VAL ("dn20: No interpolations supplied: run comp_l20_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::v20(double r){
    if(mag_star_O2_l2_status>=2){
        if(r<=star->R&&r>0){
            double aphi0 = mag_star->aphi0(r);
            double daphi0 = mag_star->daphi0(r);
            double expLambda = star->expLambda(r);

            return y20_of_r.f(r) + 2./3./r*( aphi0*aphi0/r + aphi0*daphi0/expLambda ) + 1./6./expLambda*daphi0*daphi0;
        }else if(r>star->R){
            return v20_ext(r,Qm0,1);
        }else{
            return 0;
        }
    }else{
        GSL_ERROR_VAL ("y20: No interpolations supplied: run comp_l20_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::y20(double r){
    if(mag_star_O2_l2_status>=2){
        if(r<=star->R&&r>0){
            return y20_of_r.f(r);
        }else if(r>star->R){
            return y20_ext(r,Qm0,1);
        }else{
            return 0;
        }
    }else{
        GSL_ERROR_VAL ("y20: No interpolations supplied: run comp_l20_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::dv20(double r){
    if(mag_star_O2_l2_status>=2){
        if(r<=star->R&&r>0){
            double aphi0 = mag_star->aphi0(r);
            double daphi0 = mag_star->daphi0(r);
            double jphi0 = mag_star->jphi0(r);
            double expLambda = star->expLambda(r);
            double dLambda = star->dLambda(r);
            double dNu = star->dLambda(r);

            return y20_of_r.df(r) + aphi0*daphi0/3./expLambda/r/r*( 6.*expLambda -r*(dNu+dLambda) - 2.)
                   - M_4PI*jphi0*(2.*aphi0-r*daphi0)/3./r - daphi0*daphi0*(r*dNu-4.)/expLambda/6./r;
        }else if(r>star->R){
            GSL_ERROR_VAL ("dy20: Not implemented for r>R.", GSL_FAILURE,GSL_FAILURE);
        }else{
            return 0;
        }
    }else{
        GSL_ERROR_VAL ("dy20: No interpolations supplied: run comp_l20_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::m20(double r){
    if(mag_star_O2_l2_status>=2){
        double daphi0 = mag_star->daphi0(r);
        double expLambda = star->expLambda(r);

        return ( -n20(r) + 2./3./expLambda*daphi0*daphi0 )*r/expLambda;
    }else{
        GSL_ERROR_VAL ("m20: No interpolations supplied: run comp_l20_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::dm20(double r){
    if(mag_star_O2_l2_status>=2){
        double aphi0 = mag_star->daphi0(r);
        double daphi0 = mag_star->daphi0(r);
        double jphi0 = mag_star->jphi0(r);
        double expLambda = star->expLambda(r);
        double dLambda = star->dLambda(r);
        double dNu = star->dLambda(r);

        if(r>0){
            return ( n20(r)*(-1+r*dLambda) - r*dn20(r)  - 2./3./expLambda*daphi0*daphi0*(r*(dLambda+dNu)-1)
                     + 8.*daphi0*(aphi0-M_2PI*r*r*jphi0)/3./r)/expLambda;
        }else{
            return 0;
        }

    }else{
        GSL_ERROR_VAL ("dm20: No interpolations supplied: run comp_l20_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::h20(double r){
    if(mag_star_O2_l2_status>=2){
        if(r>0){
            return -n20(r)- 2/3.*(mag_star->cjphi0)*(mag_star->aphi0(r));;
        }else{
            return 0;
        }
    }else{
        GSL_ERROR_VAL ("h20: No interpolations supplied: run comp_l20_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

double Mag_Star_O2::dh20(double r){
    if(mag_star_O2_l2_status>=2){
        if(r>0){
            return -dn20(r)- 2/3.*(mag_star->cjphi0)*(mag_star->daphi0(r));;
        }else{
            return 0;
        }
    }else{
        GSL_ERROR_VAL ("dh20: No interpolations supplied: run comp_l20_int().", GSL_FAILURE,GSL_FAILURE);
    }
}

