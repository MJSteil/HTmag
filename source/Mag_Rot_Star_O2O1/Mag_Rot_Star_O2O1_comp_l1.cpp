//
// Integration of the l=1 O[B**2Omega**1] structure equations
// Created by M. J. Steil on 2017.04.05.
//

#include "../../include/Mag_Rot_Star_O2O1.hpp"

void Mag_Rot_Star_O2O1::comp_l10(int mode, int print, vector<double> comp_params) {


    //region W10 ODEs
    function<double (double, const double *, double *, int)> dl1_eq = [&](double r , const double * uv1, double *duv1, int i)->double{
        double df_i_dr;

        double rr = r*r;
        double r3 = r*rr;
        double r4 = r*r3;

        double dNu = star->dNu(r);
        double dLambda = star->dLambda(r);

        // uv1 = {u1H, v1H, u1P, v1P}
        switch (i)
        {
            case 0: // du1H/dr
                df_i_dr = 2.*r3*( dNu + dLambda )*uv1[1];
                break;
            case 1: // dv1H/dr
                df_i_dr = -0.5*( dNu + dLambda )*uv1[1] + uv1[0]/r4;
                break;
            case 2: // du1P/dr
                df_i_dr = 2.*r3*( dNu + dLambda )*uv1[3] - ( S0(r) - S2(r) + S1(r) );
                break;
            case 3: // dv1P/dr
                df_i_dr = -0.5*( dNu + dLambda )*uv1[3] + uv1[2]/r4;
                break;
            default:
                GSL_ERROR_VAL ("dl1_eq:", GSL_FAILURE,GSL_FAILURE);
                df_i_dr = 0;
        }
        return df_i_dr;

    };

    gsl_odeiv2_ode_system uv1_eq(4);
    uv1_eq.df= &dl1_eq;
    struct gsl_odeiv2_ode_struct comp_l1_odes = {&uv1_eq};
    //endregion

    //region Computational Parameters
    if(comp_params.empty()-1){
        comp_l1_params = comp_params;
    }

    const double comp_l1_delta_r_0 = comp_l1_params[0];    // [km], Initial step size
    const double comp_l1_delta_r_max =  star->R*comp_l1_params[1]; // [km], Maximal step size
    const double comp_l1_rMax = star->R;            // [km], Maximum value for radial integration

    double comp_l1_odeiv2_control_params[4] = {comp_l1_params[2],comp_l1_params[3],comp_l1_params[4],comp_l1_params[5]};  // [error_rel, error_abs, error_ay, error_day]
    double comp_l1_error_scale[4]           = {comp_l1_params[6],comp_l1_params[7],comp_l1_params[8],comp_l1_params[9]};  // Error scale for *comp_l1_control

    unsigned long comp_l1_nMax = (unsigned long)comp_l1_params[10];   // Maximum number of Iteration steps of gsl_odeiv2_evolve_apply
    //endregion

    // gsl_odeiv2_solver setup
    gsl_odeiv2_solver comp_l1_solver(4,mode,print);
    comp_l1_solver.set_parameters(comp_l1_delta_r_0,
                                  comp_l1_delta_r_max,
                                  comp_l1_rMax,
                                  comp_l1_nMax,
                                  (double *)comp_l1_odeiv2_control_params,
                                  (double *)comp_l1_error_scale
    );
    comp_l1_solver.set_system(&comp_l1_odes,gsl_odeiv2_step_rk8pd);

    // gsl_odeiv2_data data container
    gsl_odeiv2_data comp_l1_data(8,comp_l1_nMax);

    double uv1Hscale = comp_l1_params[11];
    double Pc = star->PC;
    double rhoC = star->rhoC;
    double expNuC2 = sqrt(star->expNu(0));
    double omegabarC = rot_star->omegabar0C;

    double uv1[4] = {uv1Hscale*gsl_pow_5(comp_l1_delta_r_0),
                     uv1Hscale*( -gsl_pow_2(comp_l1_delta_r_0) + 5./(2*M_PI*(Pc+rhoC)))/8.,
                     gsl_pow_5(comp_l1_delta_r_0),
                     1/(expNuC2*16.)*( -2*gsl_pow_2(comp_l1_delta_r_0)*( expNuC2 + 4*omegabarC ) + ( 5*expNuC2 + 4* omegabarC )/( Pc + rhoC )/M_PI )
                    };
    double duv1[4] = {5*uv1Hscale*gsl_pow_4(comp_l1_delta_r_0),
                      -4*uv1Hscale*comp_l1_delta_r_0,
                      5*gsl_pow_4(comp_l1_delta_r_0),
                      -1/(expNuC2*4)*comp_l1_delta_r_0*( expNuC2 + 4*omegabarC )
    };
    double comp_l1_r = comp_l1_delta_r_0;
    double comp_l1_delta_r = comp_l1_delta_r_0;

    unsigned int comp_l1_n=0;

    // Central Values
    comp_l1_data.put(comp_l1_n,{0,comp_l1_delta_r_0,0,0,uv1Hscale*5./(16*M_PI*(Pc+rhoC)),0,0,0,( 5*expNuC2 + 4* omegabarC )/( Pc + rhoC )/(16.*M_PI)/expNuC2,0});
    comp_l1_n++;
    comp_l1_data.put(comp_l1_n,{comp_l1_r,comp_l1_delta_r,uv1[0],duv1[0],uv1[1],duv1[1],uv1[2],duv1[2],uv1[3],duv1[3]});
    comp_l1_n++;

    //region Main comp_l1_solver.comp loop
    int comp_l1_status;
    comp_l1_solver.set_IC(comp_l1_n,
                          comp_l1_r,
                          comp_l1_delta_r,
                          uv1,
                          duv1,
                          &comp_l1_data);// Inital conditions for integration

    comp_l1_solver.comp();// Integration
    comp_l1_solver.get_nrh(comp_l1_n,comp_l1_r,comp_l1_delta_r);
    nl1 = comp_l1_n;
    //endregion


    //region Solve LGS:
    double R = star->R;
    double RRRR = gsl_pow_4(R);

    double dW1HR = uv1[0]/RRRR;
    double W1HR = uv1[1];

    double dW1PR = uv1[2]/RRRR;
    double W1PR = uv1[3];

    double comp_l1_lgs_m_data[] = {W10_ext(R, 1, 0),-W1HR, dW10_ext(R, 1, 0),-dW1HR};
    double comp_l1_lgs_b_data[] = {W1PR- W10_ext(R, 0, 1),dW1PR- dW10_ext(R, 0, 1)};

    gsl_matrix_view comp_l1_lgs_m = gsl_matrix_view_array (comp_l1_lgs_m_data, 2, 2);
    gsl_vector_view comp_l1_lgs_b = gsl_vector_view_array (comp_l1_lgs_b_data, 2);

    gsl_vector *comp_l1_lgs_sol = gsl_vector_alloc (2);
    gsl_vector *comp_l1_lgs_norm = gsl_vector_alloc (2);
    gsl_vector *comp_l1_lgs_tau = gsl_vector_alloc (2);
    gsl_permutation * comp_l1_lgs_p = gsl_permutation_alloc (2);
    int comp_l1_lgs_s;

    gsl_linalg_QRPT_decomp  (&comp_l1_lgs_m.matrix, comp_l1_lgs_tau,comp_l1_lgs_p,&comp_l1_lgs_s,comp_l1_lgs_norm);
    gsl_linalg_QRPT_solve  (&comp_l1_lgs_m.matrix, comp_l1_lgs_tau,comp_l1_lgs_p, &comp_l1_lgs_b.vector, comp_l1_lgs_sol);

    deltaJ0 = gsl_vector_get(comp_l1_lgs_sol,0);
    cHl1 = gsl_vector_get(comp_l1_lgs_sol,1);

    if(print){
        printf("|-: comp_l1_lgs: {{%.4E,%.4E},{%.4E,%.4E}}.{deltaJ0,cHl1}=={%.4E,%.4E} \n",
               comp_l1_lgs_m_data[0], comp_l1_lgs_m_data[1], comp_l1_lgs_m_data[2],
               comp_l1_lgs_m_data[3], comp_l1_lgs_b_data[0], comp_l1_lgs_b_data[1]);
        printf("|-: comp_l1_lgs_sol = {deltaJ0=%.16E, cHl1=%.16E}\n",deltaJ0,cHl1);

        double comp_l1_lgs_c;
        gsl_vector *comp_l1_lgs_work = gsl_vector_alloc (6);
        gsl_linalg_QRPT_rcond (&comp_l1_lgs_m.matrix, &comp_l1_lgs_c, comp_l1_lgs_work);
        printf("|-: comp_l1_lgs_c = %.8E (=1/gsl_linalg_QRPT_rcond)\n",1/ comp_l1_lgs_c);

    };
    //endregion


    for (int i=2;i<=5;i++){
        comp_l1_data.scale(cHl1,i); // Scale homogeneous solution
        comp_l1_data.add(i+4,i);    // Add homogeneous and particular solution
    }

    //region Paste computation Data to star if(mode==1)
    if(mode==1&&comp_l1_solver.solver_status==5){
        // Assign temporary method data to star
        comp_l1_data.copy(0,&rl1_i);

        comp_l1_data.copy(2,&u1H0_i);
        comp_l1_data.copy(3,&du1H0_i);
        comp_l1_data.copy(4,&v1H0_i);
        comp_l1_data.copy(5,&dv1H0_i);

        comp_l1_data.copy(6,&u10_i);
        comp_l1_data.copy(7,&du10_i);
        comp_l1_data.copy(8,&v10_i);
        comp_l1_data.copy(9,&dv10_i);

        mag_rot_star_O2O1_l1_status = 1;
    }
    //endregion
}

void Mag_Rot_Star_O2O1::comp_l10_int() {
    if (mag_rot_star_O2O1_l1_status < 1) {
        GSL_ERROR_VAL ("comp_l10_int: no valid input data - comp_l10 has not been run successfully.", GSL_FAILURE,);
    } else {
        mag_rot_star_O2O1_l1_status = 2;

        u10_of_r = gsl_hermite_spline_obj(rl1_i,u10_i,du10_i,nl1);
        u1H0_of_r = gsl_hermite_spline_obj(rl1_i,u1H0_i,du1H0_i,nl1);

        v10_of_r = gsl_hermite_spline_obj(rl1_i,v10_i,dv10_i,nl1);
        v1H0_of_r = gsl_hermite_spline_obj(rl1_i,v1H0_i,dv1H0_i,nl1);
    }
}

void Mag_Rot_Star_O2O1::comp_l10_print_params() {
    printf("|-: comp_l1_params = {n0, delta_r_max, error_rel, error_abs, error_ay, error_day, error_scale[4], nMax, uv1Hscale}\n");
    printf("|-: {");
    for(int i=0;i<11;i++){
        printf("%.4E,",comp_l1_params[i]);
    }
    printf("%.4E}\n",comp_l1_params[11]);
}