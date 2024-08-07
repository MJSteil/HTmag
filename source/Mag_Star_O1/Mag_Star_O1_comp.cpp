//
// Main computational methods of Mag_Star_O1
// Created by M. J. Steil on 2017.02.21.
//

#include "../../include/Mag_Star_O1.hpp"

void Mag_Star_O1::comp_aphi0(double comp_aphi_cjphi0, int mode, int print, vector<double> comp_params) {
    //region Method for aphi Integration => aphi0[r]
    /*      double comp_aphi_B0 = Radial Component of central magnetic field in km^-1
     *      double comp_aphi_cjphi0 = Factor between cjphi and B0: cjphi=cjphi0*B0
     *      int mode = 0: Do not copy data to Rot_Star_O1
     *                 1: Default - Copy data to Rot_Star_O1
     *      int print = 0: Do not print text to console
     *                  1: Default - Print text to console
     * */
    //endregion

    //region Maxwell equation for aphi0 {d2aphi/dr2,daphi0/dr} integrated with B0=1
    function<double (double, const double *, double *, int)> daphi_eq = [&](double r , const double * aphi, double *daphi, int i)->double{
        double df_i_dr;
        double dLambda, dNu, expLambda, P, rho, rr;

        if (star->static_star_status == 3) {
            switch (i)
            {
                case 0: // daphi0/dr
                    df_i_dr = aphi[1];
                    break;
                case 1: // d2aphi/d2r
                    dLambda      = star->dLambda(r);
                    dNu          = star->dNu(r);
                    expLambda    = star->expLambda(r);
                    P            = star->P(r);
                    rho          = star->rho(r);
                    rr           = r*r;

                    df_i_dr = 0.5*( dLambda - dNu )*aphi[1] + expLambda*( 2*aphi[0]/rr - M_4PI*cjphi0*rr*( P + rho ));
                    break;
                default:
                    GSL_ERROR_VAL ("daphi_eq:", GSL_FAILURE,GSL_FAILURE);
                    df_i_dr = 0;
            }
            return df_i_dr;
        }else {
            // DGL Static_Star Error
            GSL_ERROR ("daphi_eq: No fully computed background star supplied: initialize with a pointer to a star with static_star_status=3.",GSL_FAILURE);
            return GSL_FAILURE;
        }

    };

    gsl_odeiv2_ode_system aphi_eq(2);
    aphi_eq.df= &daphi_eq;
    struct gsl_odeiv2_ode_struct comp_aphi_ode = {&aphi_eq};
    //endregion

    // Computational Parameters
    if(comp_params.empty()-1){
        comp_aphi0_params = comp_params;
    }

    const double comp_aphi_delta_r_0 = comp_aphi0_params[0];   // [km], Initial step size
    const double comp_aphi_delta_r_max = comp_aphi0_params[1]*star->R;   // [km], Maximal step size
    const double comp_aphi_rMax = star->R;    // [km], Maximum value for radial integration

    double comp_aphi_odeiv2_control_params[4]   = {comp_aphi0_params[2],comp_aphi0_params[3],comp_aphi0_params[4],comp_aphi0_params[5]};  // [error_rel, error_abs, error_ay, error_day]
    double comp_aphi_error_scale[2]             = {comp_aphi0_params[6], comp_aphi0_params[7]};     // Error scale for gsl_odeiv2_solver

    unsigned long comp_aphi_nMax = (unsigned long)comp_aphi0_params[8]; // Maximum number of Iteration steps of gsl_odeiv2_evolve_apply

    mu0     = 1;
    cjphi0  = comp_aphi_cjphi0;

    // gsl_odeiv2_solver setup
    gsl_odeiv2_solver comp_aphi_solver(2,mode,print);
    comp_aphi_solver.set_parameters(comp_aphi_delta_r_0,
                                    comp_aphi_delta_r_max,
                                    comp_aphi_rMax,
                                    comp_aphi_nMax,
                                    (double *)comp_aphi_odeiv2_control_params,
                                    (double *)comp_aphi_error_scale
    );
    comp_aphi_solver.set_system(&comp_aphi_ode);

    // gsl_odeiv2_data data container
    gsl_odeiv2_data comp_aphi_data(4,comp_aphi_nMax);

    double aphi[2] = {-1/2.*gsl_pow_2(comp_aphi_delta_r_0),-1*comp_aphi_delta_r_0};
    double daphi[2]= {-1*comp_aphi_delta_r_0,-1};
    double comp_aphi_r = comp_aphi_delta_r_0;
    double comp_aphi_delta_r = comp_aphi_delta_r_0;

    unsigned int comp_aphi_n=0;

    // Central Values
    comp_aphi_data.put(comp_aphi_n,{0,0,0,0,0,0});
    comp_aphi_n++;
    comp_aphi_data.put(comp_aphi_n,{comp_aphi_r,comp_aphi_delta_r,aphi[0],aphi[1],daphi[0],daphi[1]});
    comp_aphi_n++;

    //region Main comp_aphi_solver.comp loop
    comp_aphi_solver.set_IC(comp_aphi_n,
                            comp_aphi_r,
                            comp_aphi_delta_r,
                            aphi,
                            daphi,
                            &comp_aphi_data);// Inital conditions for integration
    comp_aphi_solver.comp();// Integration
    comp_aphi_solver.get_nrh(comp_aphi_n,comp_aphi_r,comp_aphi_delta_r);
    //endregion

    mu0=aphi[0]/ aphi0_ext(star->R);
    DeltadaphiRs=aphi[1]- daphi0_ext(star->R);
    n=comp_aphi_n;

    if(print){
        printf("|-: aphi integration successfull@ n=%.d - r=%.16E - h=%.8E - aphi0=%.8E- daphi0=%.8E - daphi0=%.8E(from evolve)\n|-: DeltadaphiRs=%.8E \n",
               comp_aphi_n,comp_aphi_r,comp_aphi_delta_r, aphi[0],aphi[1],comp_aphi_solver.odeiv2_evolve->dydt_out[0],DeltadaphiRs);
    };

    //region Paste computation Data to star if(mode==1)
    if(mode==1){
        if(comp_aphi_solver.solver_status==5){
            // Assign temporary method data to star
            comp_aphi_data.copy(0,&r_i);

            comp_aphi_data.copy(2,&aphi0_i);
            comp_aphi_data.copy(3,&daphi0_i);
            comp_aphi_data.copy(3,&ddaphi0_i);
        }
    }
    //endregion

}

void Mag_Star_O1::comp_aphi0_print_params() {
    printf("|-: comp_aphi0_params = {h0, delta_r_max, error_rel, error_abs, error_ay, error_day, error_scale[2], nMax}\n");
    printf("|-: {");
    for(int i=0;i<8;i++){
        printf("%.4E,",comp_aphi0_params[i]);
    }
    printf("%.4E}\n",comp_aphi0_params[8]);
}



void Mag_Star_O1::comp_cjphi0(int print,vector<double> comp_cjphi0_range_in, vector<double> comp_params) {

    if(comp_cjphi0_range_in.empty()-1){
        comp_cjphi0_range = comp_cjphi0_range_in;
    }

    // gsl_function via template wrapper gsl_function_pp
    auto comp_cjphi0_lambda = [=](double f0)->double{
        comp_aphi0(f0, 0, 0, comp_params);
        return DeltadaphiRs;
    };

    gsl_function_pp<decltype(comp_cjphi0_lambda)> comp_cjphi0_fkt_pp(comp_cjphi0_lambda);
    gsl_function *comp_cjphi0_fkt = static_cast<gsl_function*>(&comp_cjphi0_fkt_pp);

    double Z = star->M/star->R;
    vector<double> comp_cjphi0_range_Z = {comp_cjphi0_range[0]/Z,comp_cjphi0_range[1]/Z,comp_cjphi0_range[2]/Z};

    gsl_rootfinder star_target(comp_cjphi0_fkt, 200, 1e-16,0, comp_cjphi0_range_Z,print);
    cjphi0 = star_target.x_min;
    mag_star_O1_status = 1;
}

void Mag_Star_O1::comp_aphi0_int() {
    if (mag_star_O1_status < 1) {
        GSL_ERROR_VAL ("comp_aphi0_int: no valid input data - comp_aphi0 has not been run successfully.", GSL_FAILURE,);
    } else {
        mag_star_O1_status = 2;
        aphi0_of_r = gsl_hermite_spline_obj(r_i,aphi0_i,daphi0_i,n);
        daphi0_of_r = gsl_hermite_spline_obj(r_i,daphi0_i,ddaphi0_i,n);
    }
}

void Mag_Star_O1::comp(double B0_in, int mode, int print,vector<double> comp_params , vector<double> comp_cjphi0_range_in) {
    //B0 = B0_in;
    set_B0Q0(B0_in);

    comp_cjphi0(print==2,comp_cjphi0_range_in,comp_params);
    comp_aphi0(cjphi0, mode, print==2,comp_params);

    if(mode){
        comp_aphi0_int();
    }
    if(print==2||print==1){
        info();
    }

}