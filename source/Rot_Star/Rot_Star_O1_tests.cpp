//
// Rot_Star_O1 test methods
// Created by M. J. Steil on 2017.03.01.
//

#include "../../include/Rot_Star_O1.hpp"

void Rot_Star_O1::comp_I_int(int print, double err_rel, double err_abs) {
    //region Method for Moment of inertia explicit Integration -> I1 and I2 || Alternative to direct value from omegabar0 boundary condition
    /* int print = 0: Do not print text to console
     *             1: Default - Print text to console
     * */
    //endregion

    if (this->rot_star_status<2) {
        GSL_ERROR_VAL ("comp_I2: Background star with interpolations required.", GSL_FAILURE,);
    }

    // gsl_function via template wrapper gsl_function_pp
    auto dI1f = [=](double r)->double{
        return sqrt(star->expLambda(r)/star->expNu(r))*gsl_pow_4(r)*(star->P(r)+star->rho(r))* omegabar0(r);
    };
    gsl_function_pp<decltype(dI1f)> dI1fp(dI1f);
    gsl_function *dI1 = static_cast<gsl_function*>(&dI1fp);

    auto dI2f = [=](double r)->double{
        return sqrt(1/star->expLambda(r)/star->expNu(r))*gsl_pow_3(r)*(star->dNu(r)+star->dLambda(r))* omegabar0(r);
    };
    gsl_function_pp<decltype(dI2f)> dI2fp(dI2f);
    gsl_function *dI2 = static_cast<gsl_function*>(&dI2fp);

    // Integration Method
    gsl_integration comp_I1_int(dI1,{0,star->R},err_rel,err_abs,"odeiv2",(int)1E6,6,print);
    gsl_integration comp_I2_int(dI2,{0,star->R},err_rel,err_abs,"odeiv2",(int)1E6,6,print);

    double I1 = 8./3.*M_PI*comp_I1_int.F;
    double I2 = 1./3.*comp_I2_int.F;

    if(print){
        printf("Moment of inertia I1: %.16E MS**3 (1-I1/I0=%.4E)\nMoment of inertia I2: %.16E MS**3 (1-I2/J0=%.4E)\n",I1/gsl_pow_3(MSkm),1-I1/I,I2/gsl_pow_3(MSkm),1-I2/I);
    }
}

double Rot_Star_O1::comp_I_intT(int print, double err_rel, double err_abs) {
    //region Method for Moment of inertia explicit Integration -> I1 and I2 || Alternative to direct value from omegabar0 boundary condition
    /* int print = 0: Do not print text to console
     *             1: Default - Print text to console
     * */
    //endregion

    if (this->rot_star_status<2) {
        GSL_ERROR_VAL ("comp_I2: Background star with interpolations required.", GSL_FAILURE,GSL_FAILURE);
    }

    gsl_function_pp_fkt dI1 = gsl_function_pp_fkt([&](double r)->double{
        return sqrt(star->expLambda(r)/star->expNu(r))*gsl_pow_4(r)*(star->P(r)+star->rho(r))* omegabar0(r);
    });

    // Integration Method
    gsl_integration comp_I1_int(dI1.gsl_fkt(),{0,star->R},err_rel,err_abs,"odeiv2",(int)1E6,6,print==2);

    double I1 = 8./3.*M_PI*comp_I1_int.F;

    if(print){
        printf("Moment of inertia I1: %.16E MS**3 (1-I1/I0=%.4E)\n",I1/gsl_pow_3(MSkm),1-I1/I);
    }
    return I1;
}