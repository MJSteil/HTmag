//
// Created by M. J. Steil on 2017.05.08.
//
#include "../../include/Mag_Star_O2.hpp"

void Mag_Star_O2::comp_l100(int mode, int print, vector<double> comp_params) {

    //region Framedragging equation for omegaQ {domegaQ/dr,d2omegaQ/dr2}
    function<double (double, const double *, double *, int)> domegaQ_eq = [&](double r , const double * omegaQ, double *domegaQ, int i)->double{
        double df_i_dr;

        double dnuPdlambda;

        if (star->static_star_status == 3) {
            switch (i)
            {
                case 0: // domegaQ00/dr
                    df_i_dr = omegaQ[1];
                    break;
                case 1: // d2omegaQ00/d2r
                    dnuPdlambda = star->dLambda(r)+star->dNu(r);
                    df_i_dr = 2.*dnuPdlambda/r*omegaQ[0] + 0.5*( -8. + r*dnuPdlambda )/r*omegaQ[1];
                    break;
                default:
                    GSL_ERROR_VAL ("domegaQ_eq:", GSL_FAILURE,GSL_FAILURE);
                    df_i_dr = 0;
            }
            return df_i_dr;
        }else {
            // DGL Static_Star Error
            GSL_ERROR ("domegaQ_eq: No fully computed background star supplied: initialize with a pointer to a star with static_star_status=3.",GSL_FAILURE);
            return GSL_FAILURE;
        }

    };

    gsl_odeiv2_ode_system omegaQ_eq(2);
    omegaQ_eq.df= &domegaQ_eq;
    struct gsl_odeiv2_ode_struct comp_l100_ode = {&omegaQ_eq};
    //endregion

    // Computational Parameters
    if(comp_params.empty()-1){
        comp_l100_params = comp_params;
    }

    const double comp_l100_delta_r_0 = comp_l100_params[0];   // [km], Initial step size
    const double comp_l100_delta_r_max = comp_l100_params[1]*star->R;   // [km], Maximal step size
    const double comp_l100_rMax = star->R;    // [km], Maximum value for radial integration
    double comp_l100_odeiv2_control_params[4] = {comp_l100_params[2],comp_l100_params[3],comp_l100_params[4],comp_l100_params[5]};  // [error_rel, error_abs, error_ay, error_day]
    double comp_l100_error_scale[2]           = {comp_l100_params[6], comp_l100_params[7]};       // Error scale for *comp_l0_control

    unsigned long comp_l100_nMax = (unsigned long)comp_l100_params[8];   // Maximum number of Iteration steps of gsl_odeiv2_evolve_apply


    // gsl_odeiv2_solver setup
    gsl_odeiv2_solver comp_l100_solver(2,mode,print);
    comp_l100_solver.set_parameters(comp_l100_delta_r_0,
                                         comp_l100_delta_r_max,
                                         comp_l100_rMax,
                                         comp_l100_nMax,
                                         (double *)comp_l100_odeiv2_control_params,
                                         (double *)comp_l100_error_scale
    );
    comp_l100_solver.set_system(&comp_l100_ode);

    // gsl_odeiv2_data data container
    gsl_odeiv2_data comp_l100_data(4,comp_l100_nMax);
    //endregion

    //region 0. and 1. explicit step
    double omegaQ[2] = {(1+8./5.*M_PI*(star->PC+star->rhoC)*comp_l100_delta_r_0*comp_l100_delta_r_0),16./5.*M_PI*(star->PC+star->rhoC)*comp_l100_delta_r_0};
    double domegaQ[2] = {16./5.*M_PI*(star->PC+star->rhoC)*comp_l100_delta_r_0,16./5.*M_PI*(star->PC+star->rhoC)};
    double comp_l100_r = comp_l100_delta_r_0;
    double comp_l100_delta_r = comp_l100_delta_r_0;

    unsigned int comp_l100_n=0;

    comp_l100_data.put(comp_l100_n,{0,0,1,0,0,16./3.*M_PI*(star->PC+star->rhoC)});
    comp_l100_n++;
    comp_l100_data.put(comp_l100_n,{comp_l100_r,comp_l100_delta_r,omegaQ[0],omegaQ[1],domegaQ[0],domegaQ[1]});
    comp_l100_n++;
    //endregion


    //region Main comp_l100_solver.comp loop
    comp_l100_solver.set_IC(comp_l100_n,
                                 comp_l100_r,
                                 comp_l100_delta_r,
                                 omegaQ,
                                 domegaQ,
                                 &comp_l100_data);// Inital conditions for integration
    comp_l100_solver.comp();// Integration
    comp_l100_solver.get_nrh(comp_l100_n,comp_l100_r,comp_l100_delta_r);
    //endregion

    //region Rescale solution
    double R = comp_l100_r;

    double omegaQ00_ext_P_R  =  omegaQ00_ext(comp_l100_r,0);
    double domegaQ00_ext_P_R = domegaQ00_ext(comp_l100_r,0);

    double comp_l100_omegaQScale = ( 3.*omegaQ00_ext_P_R + domegaQ00_ext_P_R*R )/( omegaQ[1]*R + 3.*omegaQ[0] );
    deltaJ20_00 = -gsl_pow_4(R)*( omegaQ00_ext_P_R*omegaQ[1] - domegaQ00_ext_P_R*omegaQ[0] )/( 2.*omegaQ[1]*R+ 6.*omegaQ[0] );

    omegaQ[0] = omegaQ[0]*comp_l100_omegaQScale;
    omegaQ[1] = omegaQ[1]*comp_l100_omegaQScale;
    //endregion

    if(print){printf("|-: l100 integration successfull@ n=%.d - r=%.8E - h=%.8E - omegaQ=%.8E- domegaQ=%.8E \n",
                     comp_l100_n,comp_l100_r,comp_l100_delta_r,
                     omegaQ[0],omegaQ[1]);};

    nl1=comp_l100_n;

    //region Paste computation Data to star if(mode==1)
    if(comp_l100_solver.solver_status==5&&mode==1){
        mag_star_O2_l1_status=1;
        comp_l100_data.copy(0,&rl1_i);

        comp_l100_data.scale(comp_l100_omegaQScale,2);
        comp_l100_data.scale(comp_l100_omegaQScale,3);
        comp_l100_data.scale(comp_l100_omegaQScale,4);
        comp_l100_data.copy(2,&omegaQ00_i);
        comp_l100_data.copy(3,&domegaQ00_i);
        comp_l100_data.copy(4,&ddomegaQ00_i);
    }
    //endregion
}

void Mag_Star_O2::comp_l100_int() {
    if (mag_star_O2_l1_status < 1) {
        GSL_ERROR_VAL ("comp_l100_int: no valid input data - comp_l100 has not been run successfully.", GSL_FAILURE,);
    } else {
        mag_star_O2_l1_status = 2;

        omegaQ00_of_r = gsl_hermite_spline_obj(rl1_i,omegaQ00_i,domegaQ00_i,nl1);
        domegaQ00_of_r = gsl_hermite_spline_obj(rl1_i,domegaQ00_i,ddomegaQ00_i,nl1);
    }
}

void Mag_Star_O2::comp_l100_print_params() {
    printf("|-: comp_l100_params = {h0, hmax, error_rel, error_abs, error_ay, error_day, error_scale[2], nMax}\n");
    printf("|-: {");
    for(int i=0;i<8;i++){
        printf("%.4E,",comp_l100_params[i]);
    }
    printf("%.4E}\n",comp_l100_params[8]);
}

// l=1 tests
double  Mag_Star_O2::deltaJ20_error(int print){
    return  abs(1-deltaJ20_int(print)/deltaJ20());
}

double Mag_Star_O2::deltaJ20_int(int print) {

    gsl_function_pp_fkt test_deltaJ20_dJ = gsl_function_pp_fkt([&](double r)->double{
        double r4 = r*r*r*r;

        double expLambda = star->expLambda(r);
        double expNu = star->expNu(r);
        double P = star->P(r);
        double rho = star->rho(r);

        if(r>0){
            return -8./3.*sqrt(expLambda/expNu)*M_PI*r4*(P+rho)*omegaQ00(r);
        }else{
            return 0;
        }
    });

    //region Integration
    gsl_integration test_deltaJ20_J_i(test_deltaJ20_dJ.gsl_fkt(),{0,star->R},1E-16,1E-16,"odeiv2",(int)1E6,6,print==2);

    double M = star->M;
    double R = star->R;
    double test_deltaJ20_J_e_ana = -mag_star->mu0*( 2.*M*(M+R) + R*R*log(1-2.*M/R) )/(4*M*M*M);

    double deltaJ20_int = (test_deltaJ20_J_e_ana +test_deltaJ20_J_i.F)*mag_star->Qe()*mag_star->B0;
    //endregion

    //region Printing
    if(print)printf("deltaJ20_int=%.16E Msol**2, |1-deltaJ20_int/deltaJ20()|=%.4E\n",deltaJ20_int/MSkm/MSkm,1-deltaJ20_int/deltaJ20());
    if(print)printf("deltaJ20_int=%.16E kg m**2 s**-1 = (int,ext)=(%.2f, %.2f) \n",deltaJ20_int/ckgm2s1km2,
           test_deltaJ20_J_i.F/(test_deltaJ20_J_e_ana+test_deltaJ20_J_i.F),test_deltaJ20_J_e_ana/(test_deltaJ20_J_e_ana+test_deltaJ20_J_i.F));
    //endregion

    return deltaJ20_int;

}