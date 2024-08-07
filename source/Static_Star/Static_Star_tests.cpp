//
// Static_Star test methods
// Created by M. J. Steil on 2017.03.01.
//

#include "../../include/Static_Star.hpp"

void Static_Star::comp_MGm1(int print) {
    //region Method for gravitational mass direct integration -> MG ...
    /* int print = 0: Do not print text to console
     *             1: Default - Print text to console
     * */
    //endregion

    if (static_star_status<3) {
        GSL_ERROR_VAL ("comp_MGm1: Background star with interpolations required.", GSL_FAILURE,);
    }

    // gsl_function via template wrapper gsl_function_pp
    auto comp_MGm1_dMG = [=](double r)->double{
        return 4 * M_PI * r*r * rho(r);
    };
    gsl_function_pp<decltype(comp_MGm1_dMG)> comp_MGm1_dMGp(comp_MGm1_dMG);
    gsl_function *comp_Mm1_dMG = static_cast<gsl_function*>(&comp_MGm1_dMGp);

    // Integration Method
    gsl_integration comp_Mm1_int(comp_Mm1_dMG,{0,R},1E-16,1E-16,"odeiv2",10000,6,print);

    if(print){
        printf("|=>I1: Gravitational mass MG: %.8E (error: %.8E) MS | 1-MG_I1/M = %.8E\n",comp_Mm1_int.F/MSkm,comp_Mm1_int.F_error/MSkm,1-comp_Mm1_int.F/M);
    }
}

double Static_Star::comp_MGm2(int print, double err_rel, double err_abs) {
    //region Method for gravitational mass direct integration -> MG ...
    /* int print = 0: Do not print text to console
     *             1: Default - Print text to console
     * */
    //endregion

    if (static_star_status<3) {
        GSL_ERROR_VAL ("comp_MGm2: Background star with interpolations required.", GSL_FAILURE,GSL_FAILURE);
    }

    gsl_function_pp_fkt comp_MGm2_dMG = gsl_function_pp_fkt([&](double r)->double{
        return r*r * (rho(r)+3.*P(r))*sqrt(expNu(r)*expLambda(r));
    });

    // Integration Method
    gsl_integration comp_Mm2_int(comp_MGm2_dMG.gsl_fkt(),{0,R},err_rel,err_abs,"odeiv2",(int)1E6,2,print);

    if(print){
        printf("|=>I2: Gravitational mass MG_I2: %.8E (error: %.8E) MS | 1-MG_I2/M = %.8E\n",M_4PI*comp_Mm2_int.F/MSkm,comp_Mm2_int.F_error/MSkm,1-M_4PI*comp_Mm2_int.F/M);
    }

    return M_4PI*comp_Mm2_int.F;
}

void Static_Star::comp_MR(double hC_min, double hC_max, double pth, double drho, int &type){
    vector<double> hCi, Ri, Mi,RiP,MiP;
    int nh=0;
    double h =hC_min;
    type=0;

    while (h < (1+hC_max)*hC_min) {
        hCi.push_back(h);
        comp_tov(h, 0, 0,{1E-16, 1E-3, 1.E-16, 1.E-16, 1., 1., 1, 1.E-3, 5E3});
        Mi.push_back(M / MSkm);
        Ri.push_back(R);

//        comp_tov_inP(h);
//        MiP.push_back(M / MSkm);
//        RiP.push_back(R);
        h = h*(1+hC_max*0.1);
        if(nh==0){
            printf("nh=%d: PC=%.16E MeVfm**-3 => (%.16E MS,%.16E km)\n",nh,PC/cMeVfm3km2,M/MSkm,R);
        }
        nh++;
    }
    printf("nh=%d: PC=%.16E MeVfm**-3 => (%.16E MS,%.16E km)\n",nh,PC/cMeVfm3km2,M/MSkm,R);
/*

    gsl_interp_accel *polyMRacc = gsl_interp_accel_alloc();
    gsl_spline *polyMR = gsl_spline_alloc (gsl_interp_polynomial, 4);

    double RiPolyData[4] = {Ri[nh-1],Ri[2*nh/3],Ri[1*nh/3],Ri[0]};
    double MiPolyData[4] = {Mi[nh-1],Mi[2*nh/3],Mi[1*nh/3],Mi[0]};

    gsl_spline_init (polyMR, RiPolyData, MiPolyData, 4);

    double dx0, dxm, ddx0, ddx3;
    dx0 = gsl_spline_eval_deriv (polyMR,RiPolyData[3],polyMRacc);
    dxm = gsl_spline_eval_deriv (polyMR,(RiPolyData[1]+RiPolyData[2])/2.,polyMRacc);
    ddx0 = gsl_spline_eval_deriv2 (polyMR,RiPolyData[3],polyMRacc);
    ddx3 = gsl_spline_eval_deriv2 (polyMR,RiPolyData[0],polyMRacc);

    if(dx0>0&&dxm>0){
        type=4; //D
    }
    if(dx0>0&&dxm<0){
        type=1; //A
    }
    if(dx0<0){
        type=3; //C
    }

    printf("%.4E, %.4E , %.4E ,%.4E => %d\n",dx0,dxm,ddx0,ddx3,type);
*/

/*    Plot test(Ri,Mi*//*,RiP,MiP*//*);
    stringstream name_tmp;
    name_tmp << "'P_th/rho_th="<<pth << ", deltaRho/rho_th="<<drho<<", type="<< type<<"'" ;

    test.setTitle(name_tmp.str());
    test.plot();*/
}

void Static_Star::test(double *array){

    double test[4] = {1.,1.,1.,1.};
    array=test;

    auto dfa = [=](double r)->double{
        return dh(r)*drhodh(r);
    };


    gsl_function_pp<decltype(dfa)> dfp(dfa);
    gsl_function *df = static_cast<gsl_function*>(&dfp);

    vector<double> range(r_i.begin(),r_i.begin() + 1300);


    gsl_integration comp_deltaMint0_I2(df,range,1E-16,1E-10,"trap",(int)1E6,4,1);

    printf("%.8E, %.8E \n",(rhoC+comp_deltaMint0_I2.F)/rho(range[1299])-1,range[1299]);

}

void Static_Star::comp_tov_inr(double hC_in){
    hC = hC_in;
    PC = star_eos->P(hC);
    rhoC = star_eos->rho(hC);
    nbarC = star_eos->nbar(hC);
    drhodhC = star_eos->drhodh(hC);

    //region TOV equation {dP/dr,dM/dr,dnu/dr}
    function<double (double, const double *, double *, int)> dPMeq = [&](double r , const double * PM, double *dPdM, int i)->double{
        double df_i_dr;
            switch (i)
            {
                case 0: // dP/dr
                    df_i_dr = -(PM[1]+4*M_PI *gsl_pow_3(r)*PM[0])*(star_eos->rhoofP(PM[0])+PM[0])/(r*(r-2*PM[1]));
                    break;
                case 1: // dM/dr
                    df_i_dr = 4*M_PI*star_eos->rhoofP(PM[0])*gsl_pow_2(r);
                    break;
                case 2: // dnu/dr
                    df_i_dr = -2*dPdM[0]/(PM[0]+star_eos->rhoofP(PM[0]));
                    break;
                default:
                    GSL_ERROR_VAL ("comp_tov_eq:", GSL_FAILURE,GSL_FAILURE);
                    df_i_dr = 0;
            }
            return df_i_dr;
    };

    gsl_odeiv2_ode_system tov_eq(3);
    tov_eq.df= &dPMeq;
    struct gsl_odeiv2_ode_struct comp_tov_params = {&tov_eq};
    //endregion

    //region Computational Parameters
    const double comp_tov_h0 = 1e-20;       // [km], Initial step size
    const double comp_tov_hmax = 1E5;        // [km], Maximal step size
    const double comp_tov_hmin = 1e-16;     // [km], Minimal step size
    const double comp_tov_hmin_n = 100;     // Minimal step size step threshold
    const double comp_tov_rMax = 50;        // [km], Maximum value for radial integration
    const int comp_tov_nMax = 10000;         // Maximum number of iteration steps of gsl_odeiv2_evolve_apply
    const size_t comp_tov_P0polyNodes = 10; // Nodes for the interpolation polynomial for the P<0 surface

    const double comp_tov_error_scale[3] = {1.E-35, 1.E-10, 1.E-10};  // Error scale for gsl_odeiv2_control_scaled_new
    const double comp_tov_error_rel = 1e-16;  // Relative error parameter for gsl_odeiv2_control_scaled_new
    const double comp_tov_error_abs = 1e-10;  // Absolute error parameter for gsl_odeiv2_control_scaled_new
    const double comp_tov_error_ay = 1;  // ay for gsl_odeiv2_control_scaled_new
    const double comp_tov_error_day = 1;  // day for gsl_odeiv2_control_scaled_new
    //endregion

    //region gsl_odeiv2 setup
    const gsl_odeiv2_step_type *comp_tov_type = comp_tov_step_type;    // Step type algorithm: "Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method."
    gsl_odeiv2_step *comp_tov_step = gsl_odeiv2_step_alloc(comp_tov_type, 3);
    gsl_odeiv2_control *comp_tov_control = gsl_odeiv2_control_scaled_new(comp_tov_error_rel, comp_tov_error_abs,
                                                                         comp_tov_error_ay, comp_tov_error_day,
                                                                         comp_tov_error_scale,
                                                                         3); // Adaptive Step-size Control
    gsl_odeiv2_evolve *comp_tov_evolve = gsl_odeiv2_evolve_alloc(3);
    gsl_odeiv2_system comp_tov_system = {gsl_odeiv2_ode, NULL, 3, &comp_tov_params};
    //endregion

    //region Initial steps 0 and 1
    unsigned int comp_tov_n = 0;

    // Temporary vector<double> storage container
    gsl_odeiv2_data comp_tov_data(6,comp_tov_nMax);

    // Explicit first step of iteration using TOV eq. r-> asymptotic's: Sec. r->0 Assymptotic: ISS - r0TOV-rules, Pertubation Formulas rotO1magO2.nb
    double r = comp_tov_h0;
    double h = comp_tov_h0;
    r0 = comp_tov_h0;
    // P_1, M_1, nu_1
    double PM[3] = {PC-2./3*M_PI*gsl_pow_2(r)*(PC+rhoC)*(3*PC+rhoC),
                    4./3*M_PI*gsl_pow_3(r)*rhoC,
                    4./3*M_PI*gsl_pow_2(r)*(3*PC+rhoC)};

    // dP_1, dM_1, dnu_1
    (comp_tov_evolve->dydt_in)[0] = -4./3*M_PI*r*(PC+rhoC)*(3*PC+rhoC);
    (comp_tov_evolve->dydt_in)[1] = 4*M_PI*gsl_pow_2(r)*(rhoC);
    (comp_tov_evolve->dydt_in)[2] = 8./3*M_PI*r*(3*PC+star_eos->rho(PC));

    double PM_prev[3];
    double r_pref, h_prev;
    double dydt_prev[3];
    // Central Values
    comp_tov_n++;
    comp_tov_n++;
    //endregion

    //region Main gsl_odeiv2_evolve_apply Loop
    int status=GSL_CONTINUE;

    while (1) {

        // Save initial conditions of the next step
        for(int i=0;i<3;i++){
            PM_prev[i] = PM[i];
            dydt_prev[i] = comp_tov_evolve->dydt_in[i];
        }
        r_pref = r;
        h_prev = h;

        //step
        status = gsl_odeiv2_evolve_apply(comp_tov_evolve, comp_tov_control, comp_tov_step, &comp_tov_system, &r,
                                         comp_tov_rMax, &h, PM);

        comp_tov_n++;

        // Negative Pressure Criterion and surface interpolation
        if (PM[0]<0) {
            //printf("negative Pressure@ node %d: r=%.16E km - p=%.8E MeVfm3 \n",comp_tov_n,r,PM[0]/cMeVfm3km2);
            comp_tov_n--;
            //comp_tov_evolve->failed_steps++;
            for(int i=0;i<3;i++) {
                comp_tov_evolve->dydt_in[i] = dydt_prev[i];
                PM[i] = PM_prev[i];
            }
            r= r_pref;
            h=h_prev*1E-6;
        };

        // h<hmin Criterion
        if (h<comp_tov_hmin&&comp_tov_n>comp_tov_hmin_n) {
            comp_tov_n--;
           
            printf("h Criterion@ node %d: r=%.16E km - M=%.8E km - h=%.8E km - p=%.8E MeVfm3 - rho=%.8E MeVfm3 - nu_arb=%.8E \n",comp_tov_n,r,PM[1], h, PM[0]/cMeVfm3km2,(star_eos)->rho(PM[0])/cMeVfm3km2,PM[2]);

            static_star_status=1;
            break;

        };

        // n>comp_tov_nMax: Integration failed to convert with the given maximum steps
        if(comp_tov_n>comp_tov_nMax){
            status = GSL_CONTINUE;
        }

        // r>comp_tov_rMax: Integration failed to convert with the given maximum steps
        if(r>comp_tov_rMax){
            status = GSL_CONTINUE;
        }

        // GSL Integration error
        if (status != GSL_SUCCESS) {
            comp_tov_n--;
            printf("GSL_ERROR: %s | Stopped at Step: %d - r=%.8E - p=%.8E MeVfm3 - M=%.8E km - nu_arb=%.8E \n", gsl_strerror(status), comp_tov_n, r, PM[0]/cMeVfm3km2, PM[1],PM[2]);
            break;
        }
    }
    //endregion

    n=comp_tov_n;
    count =(int)comp_tov_evolve->count;
    fails =(int)comp_tov_evolve->failed_steps;
    R=r;
    M=PM[1];
    Z=M/R;
    gsl_odeiv2_evolve_free (comp_tov_evolve);
    gsl_odeiv2_control_free (comp_tov_control);
    gsl_odeiv2_step_free (comp_tov_step);

}

void Static_Star::comp_tov_inP(double hC_in){

    // Initial values at r=0 <=> a=0
    hC = hC_in;
    PC = star_eos->P(hC);
    rhoC = star_eos->rho(hC);
    nbarC = star_eos->nbar(hC);
    drhodhC = star_eos->drhodh(hC);

    double comp_tov_odeiv2_control_params[4]   = {1E-16,1E-12,1.,1.};        // [error_rel, error_abs, error_ay, error_day]
    double comp_tov_error_scale[2]             = {1,1E-3};         // Error scale for gsl_odeiv2_solver
    unsigned long comp_tov_nMax                    = 10000;                     // Maximum number of Iteration steps of gsl_odeiv2_evolve_apply
    const double comp_tov_a_max                = PC-1E-18;//        // [1], Maximum value for h integration: stellar surface
    const double comp_tov_delta_a_max          = PC*1E+2;                    // [1], Maximal step size in h
    const double comp_tov_a1 = 1E-24;                                     // [1], Inital point in a
    const double comp_tov_delta_a_min = comp_tov_a1;                    // [1], Minimal step size in h

    //region TOV equation {drr/da,dz/da} with a=hC-h, da=-dh, rr=r**2, z=m/r
    function<double (double, const double *, double *, int)> drrda_dzda_eq = [&](double a , const double * rr_z, double *drrda_dzda, int i)->double{

        // Auxiliary variables for {drr/da,dz/da} eqs.
        double P       = PC - a;
        double rho      = star_eos->rhoofP(P);
        double drrda    = 2.*rr_z[0]*(1. - 2.*rr_z[1])/(M_4PI*P*rr_z[0] + rr_z[1])/(P+rho);

        // {dr**2/da,dm/da} eqs.
        double df_i_dr;
        switch (i)
        {
            case 0: // dr**2/da
                df_i_dr = drrda;
                break;
            case 1: // dm/r/da
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
    gsl_odeiv2_solver comp_tov_solver(2,0,1);
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

    r0 = comp_tov_a1/(2./3.*M_PI*(PC+rhoC)*(3*PC+rhoC));
    double rr_z[2]     = {r0,r0*M_4PI*rhoC/3.}; //
    double drr_z_da[2] = {0,0};

    //region gsl_odeiv2_data data container setup and first explicit steps
    gsl_odeiv2_data comp_tov_data(4,comp_tov_nMax);

/*    comp_tov_data.put(comp_tov_n,{0,comp_tov_delta_a, // a, delta_a
                                  0,0, // rr, drr/da
                                  0,0// z, dz/da
    });*/
    comp_tov_n++;

/*    comp_tov_data.put(comp_tov_n,{comp_tov_a,comp_tov_delta_a, // a, delta_a
                                  rr_z[0],drr_z_da[0], // rr, drr/da
                                  rr_z[1],drr_z_da[1] // z, dz/da
    });*/
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
    nuO = log(1-2*M/R)+2*(comp_tov_a-comp_tov_a_max);

}

void Static_Star::comp_tov_inh_rm(double hC_in) {

    // Initial values at r=0 <=> a=0
    hC = hC_in;
    PC = star_eos->P(hC);
    rhoC = star_eos->rho(hC);
    nbarC = star_eos->nbar(hC);
    drhodhC = star_eos->drhodh(hC);

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
        double drda  = rr_z[0]*(rr_z[0]-2.*rr_z[1])/(rr_z[1]+M_4PI*gsl_pow_3(rr_z[0])*P);

        // {dr/da,dm/da} eqs.
        double df_i_dr;
        switch (i)
        {
            case 0: // dr/da
                df_i_dr = drda;
                break;
            case 1: // dm/da
                df_i_dr = M_4PI*gsl_pow_2(rr_z[0])*rho*drda;
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
    gsl_odeiv2_solver comp_tov_solver(2,0,2);
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

    double rr_z[2]     = {sqrt(tov_series(comp_tov_a,0)),tov_series(comp_tov_a,2)*sqrt(tov_series(comp_tov_a,0))}; //
    double drr_z_da[2] = {0,0}; //

    r0 = rr_z[0];
    //region gsl_odeiv2_data data container setup and first explicit steps
    gsl_odeiv2_data comp_tov_data(4,comp_tov_nMax);

 /*   comp_tov_data.put(comp_tov_n,{0,comp_tov_delta_a, // a, delta_a
                                  tov_series(0,0),tov_series(0,1), // rr, drr/da
                                  tov_series(0,2),tov_series(0,3) // z, dz/da
    });*/
    comp_tov_n++;
/*
    comp_tov_data.put(comp_tov_n,{comp_tov_a,comp_tov_delta_a, // a, delta_a
                                  rr_z[0],drr_z_da[0], // rr, drr/da
                                  rr_z[1],drr_z_da[1] // z, dz/da
    });*/
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
    R=rr_z[0];
    M=rr_z[1];
    Z=M/R;
    nuO = log(1-2*M/R)+2*(comp_tov_a-comp_tov_a_max);

}