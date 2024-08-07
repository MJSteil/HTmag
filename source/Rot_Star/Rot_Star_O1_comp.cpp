//
// Main computational methods of Rot_Star_O1
// Created by M. J. Steil on 2017.02.24.
//

#include "../../include/Rot_Star_O1.hpp"

void Rot_Star_O1::comp_omegabar0(int mode, int print, vector<double> comp_params) {
    //region Method for omegabar0 Integration => omegabar0[r]
    /*      double f = Rotation frequency in Hz
     *      int mode = 0: Do not copy data to Rot_Star_O1
     *                 1: Default - Copy data to Rot_Star_O1
     *      int print = 0: Do not print text to console
     *                  1: Default - Print text to console
     * */
    //endregion

    //region Framedragging equation for omegabar0 {domegabar0/dr,d2omegabar/dr2}
    function<double (double, const double *, double *, int)> domegabar_eq = [&](double r , const double * omegabar0, double *domegabar0, int i)->double{
        double df_i_dr;
        if (star->static_star_status == 3) {
            switch (i)
            {
                case 0: // domega0/dr
                    df_i_dr = omegabar0[1];
                    break;
                case 1: // d2omega0/d2r
                    df_i_dr = ((r*omegabar0[1]+4*omegabar0[0])*(star->dNu(r)+star->dLambda(r))-8*omegabar0[1])/(2.*r);
                    break;
                default:
                    GSL_ERROR_VAL ("comp_omegabar0_eq:", GSL_FAILURE,GSL_FAILURE);
                    df_i_dr = 0;
            }
            return df_i_dr;
        }else {
            // ODE Static_Star Error
            GSL_ERROR ("domegabar_eq: No fully computed background star supplied: initialize with a pointer to a star with static_star_status=3.",GSL_FAILURE);
            return GSL_FAILURE;
        }

    };

    gsl_odeiv2_ode_system omegabar_eq(2);
    omegabar_eq.df= &domegabar_eq;
    struct gsl_odeiv2_ode_struct comp_omegabar0_ode = {&omegabar_eq};
    //endregion

    // Computational Parameters
    if(comp_params.empty()-1){
        comp_omegabar0_params = comp_params;
    }


    //region Setup of gsl_odeiv2_solver and gsl_odeiv2_data
    // Computational Parameters
    const double comp_omegabar0_delta_r_0   = comp_omegabar0_params[0];    // [km], Initial step size
    const double comp_omegabar0_delta_r_max = comp_omegabar0_params[1]*star->R;         // [km], Maximal step size: for hmax<.5 rk4 runs with same number of steps as rk8pd
    const double comp_omegabar0_rMax = star->R; // [km], Maximum value for radial integration

    double comp_omegabar0_odeiv2_control_params[4]   = {comp_omegabar0_params[2],comp_omegabar0_params[3],comp_omegabar0_params[4],comp_omegabar0_params[5]};  // [error_rel, error_abs, error_ay, error_day]
    double comp_omegabar0_error_scale[2]             = {comp_omegabar0_params[6], comp_omegabar0_params[7]};     // Error scale for gsl_odeiv2_solver

    unsigned long comp_omegabar0_nMax = (unsigned long)comp_omegabar0_params[8]; // Maximum number of Iteration steps of gsl_odeiv2_evolve_apply


    // gsl_odeiv2_solver setup
    gsl_odeiv2_solver comp_omegabar0_solver(2,mode,print);
    comp_omegabar0_solver.set_parameters(comp_omegabar0_delta_r_0,
                                        comp_omegabar0_delta_r_max,
                                        comp_omegabar0_rMax,
                                        comp_omegabar0_nMax,
                                        (double *)comp_omegabar0_odeiv2_control_params,
                                        (double *)comp_omegabar0_error_scale
    );
    comp_omegabar0_solver.set_system(&comp_omegabar0_ode);

    // gsl_odeiv2_data data container
    gsl_odeiv2_data comp_omegabar0_data(4,comp_omegabar0_nMax);
    //endregion

    //region 0. and 1. explicit step
    double omegabar0[2] = {(1+8./5.*M_PI*(star->PC+star->rhoC)*comp_omegabar0_delta_r_0*comp_omegabar0_delta_r_0),16./5.*M_PI*(star->PC+star->rhoC)*comp_omegabar0_delta_r_0};
    double domegabar0[2] = {16./5.*M_PI*(star->PC+star->rhoC)*comp_omegabar0_delta_r_0,16./5.*M_PI*(star->PC+star->rhoC)};
    double comp_omegabar0_r = comp_omegabar0_delta_r_0;
    double comp_omegabar0_delta_r = comp_omegabar0_delta_r_0;

    unsigned int comp_omegabar0_n=0;

    comp_omegabar0_data.put(comp_omegabar0_n,{0,0,1,0,0,16./3.*M_PI*(star->PC+star->rhoC)});
    comp_omegabar0_n++;
    comp_omegabar0_data.put(comp_omegabar0_n,{comp_omegabar0_r,comp_omegabar0_delta_r,omegabar0[0],omegabar0[1],domegabar0[0],domegabar0[1]});
    comp_omegabar0_n++;
    //endregion


    //region Main comp_omegabar0_solver.comp loop
    comp_omegabar0_solver.set_IC(comp_omegabar0_n,
                                comp_omegabar0_r,
                                comp_omegabar0_delta_r,
                                omegabar0,
                                domegabar0,
                                &comp_omegabar0_data);// Inital conditions for integration
    comp_omegabar0_solver.comp();// Integration
    comp_omegabar0_solver.get_nrh(comp_omegabar0_n,comp_omegabar0_r,comp_omegabar0_delta_r);
    //endregion

    //region Rescale solution
    double comp_omegabar0_OmegaScale = 1/(omegabar0[0]+star->R/3.*omegabar0[1]);

    omegabar0[0] = omegabar0[0]*comp_omegabar0_OmegaScale;
    omegabar0[1] = omegabar0[1]*comp_omegabar0_OmegaScale;
    //endregion

    if(print){printf("|-: omegabar0 integration successfull@ n=%.d - r=%.8E - h=%.8E - omegabar0=%.8E- domegabar0=%.8E \n",
                     comp_omegabar0_n,comp_omegabar0_r,comp_omegabar0_delta_r,
                     omegabar0[0],omegabar0[1]);};

    n=comp_omegabar0_n;
    steps =(int)comp_omegabar0_solver.odeiv2_evolve->failed_steps+n;
    J0=gsl_pow_4(star->R)*omegabar0[1]/6.;
    I=J0;
    omegabar0C = comp_omegabar0_OmegaScale*comp_omegabar0_data.data[2][0];

    //region Paste computation Data to star if(mode==1)
    if(comp_omegabar0_solver.solver_status==5&&mode==1){
        rot_star_status=1;
        comp_omegabar0_data.copy(0,&r_i);

        comp_omegabar0_data.scale(comp_omegabar0_OmegaScale,2);
        comp_omegabar0_data.scale(comp_omegabar0_OmegaScale,3);
        comp_omegabar0_data.scale(comp_omegabar0_OmegaScale,4);
        comp_omegabar0_data.copy(2,&omegabar0_i);
        comp_omegabar0_data.copy(3,&domegabar0_i);
        comp_omegabar0_data.copy(4,&ddomegabar0_i);
    }
    //endregion
}

void Rot_Star_O1::comp_omegabar0_print_params() {
    printf("|-: comp_omegabar0_params = {h0, hmax, error_rel, error_abs, error_ay, error_day, error_scale[2], nMax}\n");
    printf("|-: {");
    for(int i=0;i<8;i++){
        printf("%.4E,",comp_omegabar0_params[i]);
    }
    printf("%.4E}\n",comp_omegabar0_params[8]);
}

void Rot_Star_O1::comp_omegabar0_int() {
    if (rot_star_status < 1) {
        GSL_ERROR_VAL ("comp_omegabar0_int: no valid input data - comp_omegabar0 has not been run successfully.", GSL_FAILURE,);
    } else {
        rot_star_status = 2;

        omegabar0_of_r = gsl_hermite_spline_obj(r_i,omegabar0_i,domegabar0_i,n);
        domegabar0_of_r = gsl_hermite_spline_obj(r_i,domegabar0_i,ddomegabar0_i,n);
    }
}

void Rot_Star_O1::comp(int mode, int print, vector<double> comp_params) {
    comp_omegabar0(mode,print,comp_params);
    if(mode==1){
        comp_omegabar0_int();
    }
    if(print==1||print==2){
        info();
    }
}