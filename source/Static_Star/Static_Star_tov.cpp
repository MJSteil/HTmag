//
// Created by M. J. Steil on 2017.02.16.
//

#include "../../include/Static_Star.hpp"

// Static_Star TOV equation in a=h-hC

double Static_Star::tov_series(double a, int i, int order) {
    // => rr    = crr1*a*(1 + 2*crr2*a)
    // => drrda = crr1*(1 + 2*crr2*a)
    // => z     = cm1*(1 + cm2*a)
    // => dzda  = cm1*(1 + 2*cm2*a)

    double crr1, crr2, cz1, cz2;
    crr1 = 3./(2.*M_PI*(3.*PC+rhoC));
    crr2 = (3*drhodhC+15*PC-5*rhoC)/(10.*(3*PC+rhoC))*order;
    cz1  = 2.*rhoC/(3.*PC+rhoC);
    cz2  = -(5.*rhoC*(-3.*PC+rhoC)+3*drhodhC*(6.*PC+rhoC))/(10*rhoC*(3.*PC+rhoC))*order;

    double df_i_dr;
    switch (i)
    {
        case 0: //rr(a)
            df_i_dr = crr1*a*(1+crr2*a);
            break;
        case 1: //drr(a)
            df_i_dr = crr1*(1+2*crr2*a);
            break;
        case 2: //z(a)
            df_i_dr = cz1*a*(1+cz2*a);
            break;
        case 3: //dz(a)
            df_i_dr = cz1*(1+2*cz2*a);
            break;
        default:
            GSL_ERROR_VAL ("tov_series:", GSL_EDOM,GSL_EDOM);
            df_i_dr = 0;
    }
    return df_i_dr;
}

void Static_Star::comp_tov(double hC_in, int mode, int print,vector<double> comp_params){

    if(star_eos->type == "Tolman VII solution/TVII EoS"){
        star_eos->set_params({hC_in});
    }

    // Initial values at r=0 <=> a=0
    hC = hC_in;
    PC = star_eos->P(hC);
    rhoC = star_eos->rho(hC);
    nbarC = star_eos->nbar(hC);
    drhodhC = star_eos->drhodh(hC);
    dnbardhC = star_eos->dnbardh(hC);
    csC = star_eos->cs(hC);

    // Computational Parameters
    if(comp_params.empty()-1){
        comp_tov_params = comp_params;
    }

    double comp_tov_odeiv2_control_params[4]   = {comp_tov_params[2],comp_tov_params[3],comp_tov_params[4],comp_tov_params[5]};  // [error_rel, error_abs, error_ay, error_day]
    double comp_tov_error_scale[2]             = {comp_tov_params[6],comp_tov_params[7]};   // Error scale for gsl_odeiv2_solver
    unsigned long comp_tov_nMax                = (unsigned long)comp_tov_params[8];         // Maximum number of Iteration steps of gsl_odeiv2_evolve_apply
    const double comp_tov_a_max                = hC-star_eos->h0;       // [1], Maximum value for h integration: stellar surface
    const double comp_tov_delta_a_max          = hC*comp_tov_params[1]; // [1], Maximal step size in h
    const double comp_tov_a1                   = comp_tov_params[0];    // [1], Inital point in a
    const double comp_tov_delta_a_min          = comp_tov_a1;           // [1], Minimal step size in h

    //region TOV equation {drr/da,dz/da} with a=hC-h, da=-dh, rr=r**2, z=m/r
    function<double (double, const double *, double *, int)> drrda_dzda_eq = [&](double a , const double * rr_z, double *drrda_dzda, int i)->double{

        // Auxiliary variables for {drr/da,dz/da} eqs.
        double h       = hC - a;
        double P        = star_eos->P(h);
        double rho      = star_eos->rho(h);
        double drrda    = 2.*rr_z[0]*(1. - 2.*rr_z[1])/(M_4PI*P*rr_z[0] + rr_z[1]);

        // {dr**2/da,d(m/r)/da} eqs.
        double df_i_dr;
        switch (i)
        {
            case 0: // dr**2/da
                df_i_dr = drrda;
                break;
            case 1: // d(m/r)/da
                df_i_dr = (M_2PI*rho - 0.5*rr_z[1]/rr_z[0])*drrda;
                break;
            default:
                GSL_ERROR_VAL ("comp_tov_eq:", GSL_FAILURE,GSL_FAILURE);
                df_i_dr = 0;
        }
        return df_i_dr;
    };

    // Wrapper struct for gsl_odeiv2_ode_struct; necessary since gsl_odeiv2 methods need the ODE system in a very specific form.
    gsl_odeiv2_ode_system tov_eq_drrda_dzda(2);
    tov_eq_drrda_dzda.df= &drrda_dzda_eq;
    struct gsl_odeiv2_ode_struct comp_tov_struct = {&tov_eq_drrda_dzda};
    //endregion

    // gsl_odeiv2_solver setup
    gsl_odeiv2_solver comp_tov_solver(2,mode,print);
    comp_tov_solver.set_parameters(comp_tov_delta_a_min,
                                   comp_tov_delta_a_max,
                                   comp_tov_a_max,
                                   comp_tov_nMax,
                                   (double *)comp_tov_odeiv2_control_params,
                                   (double *)comp_tov_error_scale
                                   );
    comp_tov_solver.set_system(&comp_tov_struct,comp_tov_step_type);

    // Inital Values
    double comp_tov_a       = comp_tov_a1;
    double comp_tov_delta_a = comp_tov_a1;
    unsigned int comp_tov_n = 0;

    if(print){
        printf("|-: comp_tov: 1-NLO(rr_1)/LO(rr_1)=%.8E | 1-NLO(z_1)/LO(z_1)=%.8E => r_1=%.8E km\n",
               1-tov_series(comp_tov_a,0,1)/tov_series(comp_tov_a,0,0),
               1-tov_series(comp_tov_a,1,1)/tov_series(comp_tov_a,1,0),
               sqrt(tov_series(comp_tov_a,0,1))
        );
    }


    double rr_z[2]     = {tov_series(comp_tov_a,0),tov_series(comp_tov_a,2)}; //
    double drr_z_da[2] = {tov_series(comp_tov_a,1),tov_series(comp_tov_a,3)}; //


    //region gsl_odeiv2_data data container setup and first explicit steps
    gsl_odeiv2_data comp_tov_data(4,comp_tov_nMax);

    comp_tov_data.put(comp_tov_n,{0,comp_tov_delta_a, // a, delta_a
                                  tov_series(0,0),tov_series(0,1), // rr, drr/da
                                  tov_series(0,2),tov_series(0,3) // z, dz/da
                                  });
    comp_tov_n++;

    comp_tov_data.put(comp_tov_n,{comp_tov_a,comp_tov_delta_a, // a, delta_a
                                  rr_z[0],drr_z_da[0], // rr, drr/da
                                  rr_z[1],drr_z_da[1] // z, dz/da
                                  });
    comp_tov_n++;
    //endregion

    //region Main loop
    comp_tov_solver.set_IC(comp_tov_n,
                           comp_tov_a,
                           comp_tov_delta_a,
                           rr_z,
                           drr_z_da,
                           &comp_tov_data);// Inital conditions for integration

    comp_tov_solver.comp();
    comp_tov_solver.get_nrh(comp_tov_n,comp_tov_a,comp_tov_delta_a);
    //endregion

    // Compute nuc offset via nu[Rs]==Log(1-2*Ms/Rs)
    n=comp_tov_n;
    count =(int)comp_tov_solver.odeiv2_evolve->count;
    fails =(int)comp_tov_solver.odeiv2_evolve->failed_steps;
    R=sqrt(rr_z[0]);
    M=rr_z[1]*R;
    Z=M/R;
    nuO = log(1-2*M/R)+2*(comp_tov_a-comp_tov_a_max);

    //region Paste computation Data to star if(mode==1)
    if(mode==1&&comp_tov_solver.solver_status==5){
        static_star_status=2;

        // Assign temporary method data to star
        comp_tov_data.copy(0,&a_i);
        comp_tov_data.copy(2,&rr_i);
        comp_tov_data.copy(3,&drrda_i);
        comp_tov_data.copy(4,&z_i);
        comp_tov_data.copy(5,&dzda_i);
    }
    //endregion
}

void Static_Star::comp_tov_print_params() {
    printf("|-: comp_tov_params = {h0, hmax, error_rel, error_abs, error_ay, error_day, error_scale[2], nMax}\n");
    printf("|-: {");
    for(int i=0;i<8;i++){
        printf("%.4E,",comp_tov_params[i]);
    }
    printf("%.4E}\n",comp_tov_params[8]);
}

void Static_Star::comp_tov_int() {
    if (static_star_status < 1) {
        GSL_ERROR_VAL ("comp_tov_int: Comp_tov was not successful or has not been run yet for this star.", GSL_FAILURE,);
    } else {
        if (static_star_status < 2) {
            GSL_ERROR_VAL ("comp_tov_int: No data supplied: run comp_tov with parameter mode=1.", GSL_FAILURE,);
        } else {
            static_star_status=3;
            rr_of_a=gsl_hermite_spline_obj(a_i,rr_i,drrda_i,n);
            z_of_a=gsl_hermite_spline_obj(a_i,z_i,dzda_i,n);

            r_i = vector<double>(n);
            h_i = vector<double>(n);
            dhdr_i = vector<double>(n);
            dzdr_i = vector<double>(n);

            for(int i=0;i<n;i++){
                r_i[i]=sqrt(rr_i[i]);
                h_i[i]=hC-a_i[i];
                dhdr_i[i]=-2.*r_i[i]/drrda_i[i];
                dzdr_i[i]=-dzda_i[i]*dhdr_i[i];
            }
            h_of_r = gsl_hermite_spline_obj(r_i,h_i,dhdr_i,n);
            z_of_r  = gsl_hermite_spline_obj(r_i,z_i,dzdr_i,n);
        }
    }
}