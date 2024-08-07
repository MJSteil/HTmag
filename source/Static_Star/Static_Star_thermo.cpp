//
// Created by M. J. Steil on 2017.02.16.
//

#include "../../include/Static_Star.hpp"

// Static_Star thermodynamic properties

double Static_Star::h(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("h(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("h(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R&&r>0){
        return h_of_r.f(r);
    }else if (r==0){
        return hC;
    }else{
        return 0;
    }

}

double Static_Star::dh(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("dh(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("dh(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R){
        if(r<=R){
            return h_of_r.df(r);
        }else{
            return 0;
        }
    }else{
        return 0;
    }
}

double Static_Star::P(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("P(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("P(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R&&r>0){
        return star_eos->P(h(r));
    }else if (r==0){
        return PC;
    }else{
        return 0;
    }

}

double Static_Star::dP(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("dP(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("dP(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R){
        return star_eos->dPdh(h(r)) * dh(r); // dP/dh * dh/dr = dP/dr

    }else{
        return 0;
    }
}

double Static_Star::rho(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("rho(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("rho(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R&&r>0){
        return star_eos->rho(h(r));
    }else if (r==0){
        return rhoC;
    }else{
        return 0;
    }
}

double Static_Star::drhodh(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("drhodh(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("drhodh(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R){
        return star_eos->drhodh(h(r));
    }else{
        return 0;
    }
}

double Static_Star::nbar(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("nbar(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("nbar(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R&&r>0){
        return star_eos->nbar(h(r));
    }else if (r==0){
        return nbarC;
    }else{
        return 0;
    }
}

double Static_Star::dnbardh(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("dnbardh(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("dnbardh(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R){
        return star_eos->dnbardh(h(r));
    }else{
        return 0;
    }

}

double Static_Star::cs(double r) {
    //region Errors: No interpolations or negative radius
    if(static_star_status<3) {
        GSL_ERROR_VAL ("cs(r): No interpolations supplied: run comp_tov_int().",GSL_FAILURE,GSL_FAILURE);
    }
    if(r<0) {
        GSL_ERROR_VAL ("cs(r): Called with negative r.",GSL_EDOM,GSL_EDOM);
    }
    //endregion

    if(r<=R){
        return star_eos->cs(h(r),0);
    }else{
        return 0;
    }

}