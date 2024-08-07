//
// Metric potentials
// Created by M. J. Steil on 2017.02.16.
//

#include "../../include/Static_Star.hpp"

// Static_Star metric potentials (g_tt=-exp(nu), g_rr=exp(lambda)) and their derivatives

double Static_Star::expLambda(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("expLambda(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("expLambda(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R&&r>0){
        return r/(r-2.*z_of_r.f(r)*r);
    }else if (r==0){
        return 1.;
    }else{
        return 1./(1-2*M/r);
    }
}

double Static_Star::dLambda(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("dLambda(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("dLambda(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R){
        return 2.*z_of_r.df(r)/(1.-2.*z_of_r.f(r));
    }else{
        return 2.*M/(2.*M*r-r*r);
    }
}

double Static_Star::expNu(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("expNu(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("expNu(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R&&r>0){
        return exp(-2.*h_of_r.f(r)+nuO);
    }else if (r==0){
        return exp(-2.*hC+nuO);
    }else{
        return 1-2*M/r;
    }
}

double Static_Star::dNu(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("dNu(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("dNu(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R){
        return -2.*h_of_r.df(r);
    }else{
        return -2*M/(2*M*r-r*r);
    }
}

double Static_Star::j(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("j(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("j(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R){
        return sqrt(expLambda(r)*expNu(r));
    }else{
        return 1;
    }
}

double Static_Star::m(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("m(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("m(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R){
        return z_of_r.f(r)*r;
    }else{
        return M;
    }
}

double Static_Star::dm(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("dm(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("dm(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R){
        return z_of_r.df(r)*r + z_of_r.f(r); // dm/dr = r*dz/dr + z
    }else{
        return 0;
    }
}