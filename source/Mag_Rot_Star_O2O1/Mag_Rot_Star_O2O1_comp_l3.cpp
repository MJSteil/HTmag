//
// Integration of the l=3 O[B**2Omega**1] structure equations
// Created by M. J. Steil on 2017.04.05.
//

#include "../../include/Mag_Rot_Star_O2O1.hpp"

void Mag_Rot_Star_O2O1::comp_l30(int mode, int print, vector<double> comp_params) {


    //region W30 ODEs
    function<double (double, const double *, double *, int)> dl3_eq = [&](double r , const double * uv3, double *duv3, int i)->double{
        double df_i_dr;

        double rr = r*r;
        double r3 = r*rr;
        double r4 = r*r3;

        double dNu = star->dNu(r);
        double dLambda = star->dLambda(r);
        double expLambda = star->expLambda(r);

        // uv3 = {u3H, v3H, u3P, v3P}
        switch (i)
        {
            case 0: // du3H/dr
                df_i_dr = 2.*rr*( 5.*expLambda + r*(dNu + dLambda) )*uv3[1];
                break;
            case 1: // dv3H/dr
                df_i_dr = -0.5*( dNu + dLambda )*uv3[1] + uv3[0]/r4;
                break;
            case 2: // du3P/dr
                df_i_dr = 2.*rr*( 5.*expLambda + r*(dNu + dLambda) )*uv3[3] - ( S2(r) + S3(r) );
                break;
            case 3: // dv3P/dr
                df_i_dr = -0.5*( dNu + dLambda )*uv3[3] + uv3[2]/r4;
                break;
            default:
                GSL_ERROR_VAL ("dl3_eq:", GSL_FAILURE,GSL_FAILURE);
                df_i_dr = 0;
        }
        return df_i_dr;

    };

    gsl_odeiv2_ode_system uv3_eq(4);
    uv3_eq.df= &dl3_eq;
    struct gsl_odeiv2_ode_struct comp_l3_odes = {&uv3_eq};
    //endregion

    //region Computational Parameters
    if(comp_params.empty()-1){
        comp_l3_params = comp_params;
    }

    const double comp_l3_delta_r_0 = comp_l3_params[0];    // [km], Initial step size
    const double comp_l3_delta_r_max =  star->R*comp_l3_params[1]; // [km], Maximal step size
    const double comp_l3_rMax = star->R;            // [km], Maximum value for radial integration

    double comp_l3_odeiv2_control_params[4] = {comp_l3_params[2],comp_l3_params[3],comp_l3_params[4],comp_l3_params[5]};  // [error_rel, error_abs, error_ay, error_day]
    double comp_l3_error_scale[4]           = {comp_l3_params[6],comp_l3_params[7],comp_l3_params[8],comp_l3_params[9]};  // Error scale for *comp_l3_control

    unsigned long comp_l3_nMax = (unsigned long)comp_l3_params[10];   // Maximum number of Iteration steps of gsl_odeiv2_evolve_apply
    //endregion

    // gsl_odeiv2_solver setup
    gsl_odeiv2_solver comp_l3_solver(4,mode,print);
    comp_l3_solver.set_parameters(comp_l3_delta_r_0,
                                  comp_l3_delta_r_max,
                                  comp_l3_rMax,
                                  comp_l3_nMax,
                                  (double *)comp_l3_odeiv2_control_params,
                                  (double *)comp_l3_error_scale
    );
    comp_l3_solver.set_system(&comp_l3_odes,gsl_odeiv2_step_rk8pd);

    // gsl_odeiv2_data data container
    gsl_odeiv2_data comp_l3_data(8,comp_l3_nMax);

    double uv3Hscale = comp_l3_params[11];

    double uv3[4] = {uv3Hscale*gsl_pow_5(comp_l3_delta_r_0),
                     uv3Hscale*gsl_pow_2(comp_l3_delta_r_0),
                     gsl_pow_5(comp_l3_delta_r_0),
                     gsl_pow_2(comp_l3_delta_r_0)
                    };
    double duv3[4] = {5*uv3Hscale*gsl_pow_4(comp_l3_delta_r_0),
                      2*uv3Hscale*comp_l3_delta_r_0,
                      5*gsl_pow_4(comp_l3_delta_r_0),
                      2*comp_l3_delta_r_0
    };
    double comp_l3_r = comp_l3_delta_r_0;
    double comp_l3_delta_r = comp_l3_delta_r_0;

    unsigned int comp_l3_n=0;

    // Central Values
    comp_l3_data.put(comp_l3_n,{0,comp_l3_delta_r_0,0,0,0,0,0,0,0,0});
    comp_l3_n++;
    comp_l3_data.put(comp_l3_n,{comp_l3_r,comp_l3_delta_r,uv3[0],duv3[0],uv3[1],duv3[3],uv3[2],duv3[2],uv3[3],duv3[3]});
    comp_l3_n++;

    //region Main comp_l3_solver.comp loop
    int comp_l3_status;
    comp_l3_solver.set_IC(comp_l3_n,
                          comp_l3_r,
                          comp_l3_delta_r,
                          uv3,
                          duv3,
                          &comp_l3_data);// Inital conditions for integration

    comp_l3_solver.comp();// Integration
    comp_l3_solver.get_nrh(comp_l3_n,comp_l3_r,comp_l3_delta_r);
    nl3 = comp_l3_n;
    //endregion


    //region Solve LGS:
    double R = star->R;
    double RRRR = gsl_pow_4(R);

    double dW3HR = uv3[0]/RRRR;
    double W3HR = uv3[1];

    double dW3PR = uv3[2]/RRRR;
    double W3PR = uv3[3];

    double comp_l3_lgs_m_data[] = {W30_ext(R, 1, 0),-W3HR, dW30_ext(R, 1, 0),-dW3HR};
    double comp_l3_lgs_b_data[] = {W3PR- W30_ext(R, 0, 1),dW3PR- dW30_ext(R, 0, 1)};

    gsl_matrix_view comp_l3_lgs_m = gsl_matrix_view_array (comp_l3_lgs_m_data, 2, 2);
    gsl_vector_view comp_l3_lgs_b = gsl_vector_view_array (comp_l3_lgs_b_data, 2);

    gsl_vector *comp_l3_lgs_sol = gsl_vector_alloc (2);
    gsl_vector *comp_l3_lgs_norm = gsl_vector_alloc (2);
    gsl_vector *comp_l3_lgs_tau = gsl_vector_alloc (2);
    gsl_permutation * comp_l3_lgs_p = gsl_permutation_alloc (2);
    int comp_l3_lgs_s;

    gsl_linalg_QRPT_decomp  (&comp_l3_lgs_m.matrix, comp_l3_lgs_tau,comp_l3_lgs_p,&comp_l3_lgs_s,comp_l3_lgs_norm);
    gsl_linalg_QRPT_solve  (&comp_l3_lgs_m.matrix, comp_l3_lgs_tau,comp_l3_lgs_p, &comp_l3_lgs_b.vector, comp_l3_lgs_sol);

    W3asymp0 = gsl_vector_get(comp_l3_lgs_sol,0);
    cHl3 = gsl_vector_get(comp_l3_lgs_sol,1);

    if(print){
        printf("|-: comp_l3_lgs: {{%.4E,%.4E},{%.4E,%.4E}}.{W3asymp0,cHl3}=={%.4E,%.4E} \n",
               comp_l3_lgs_m_data[0], comp_l3_lgs_m_data[1], comp_l3_lgs_m_data[2],
               comp_l3_lgs_m_data[3], comp_l3_lgs_b_data[0], comp_l3_lgs_b_data[1]);
        printf("|-: comp_l3_lgs_sol = {W3asymp0=%.16E, cHl3=%.16E}\n",W3asymp0,cHl3);

        double comp_l3_lgs_c;
        gsl_vector *comp_l3_lgs_work = gsl_vector_alloc (6);
        gsl_linalg_QRPT_rcond (&comp_l3_lgs_m.matrix, &comp_l3_lgs_c, comp_l3_lgs_work);
        printf("|-: comp_l3_lgs_c = %.8E (=1/gsl_linalg_QRPT_rcond)\n",1/ comp_l3_lgs_c);

    };
    //endregion


    for (int i=2;i<=5;i++){
        comp_l3_data.scale(cHl3,i); // Scale homogeneous solution
        comp_l3_data.add(i+4,i);    // Add homogeneous and particular solution
    }

    //region Paste computation Data to star if(mode==1)
    if(mode==1&&comp_l3_solver.solver_status==5){
        // Assign temporary method data to star
        comp_l3_data.copy(0,&rl3_i);

        comp_l3_data.copy(2,&u3H0_i);
        comp_l3_data.copy(3,&du3H0_i);
        comp_l3_data.copy(4,&v3H0_i);
        comp_l3_data.copy(5,&dv3H0_i);

        comp_l3_data.copy(6,&u30_i);
        comp_l3_data.copy(7,&du30_i);
        comp_l3_data.copy(8,&v30_i);
        comp_l3_data.copy(9,&dv30_i);

        mag_rot_star_O2O1_l3_status = 1;
    }
    //endregion
}

void Mag_Rot_Star_O2O1::comp_l30_int() {
    if (mag_rot_star_O2O1_l3_status < 1) {
        GSL_ERROR_VAL ("comp_l30_int: no valid input data - comp_l30 has not been run successfully.", GSL_FAILURE,);
    } else {
        mag_rot_star_O2O1_l3_status = 2;

        u30_of_r = gsl_hermite_spline_obj(rl3_i,u30_i,du30_i,nl3);
        u3H0_of_r = gsl_hermite_spline_obj(rl3_i,u3H0_i,du3H0_i,nl3);

        v30_of_r = gsl_hermite_spline_obj(rl3_i,v30_i,dv30_i,nl3);
        v3H0_of_r = gsl_hermite_spline_obj(rl3_i,v3H0_i,dv3H0_i,nl3);
    }
}

void Mag_Rot_Star_O2O1::comp_l30_print_params() {
    printf("|-: comp_l1_params = {n0, delta_r_max, error_rel, error_abs, error_ay, error_day, error_scale[4], nMax, uv1Hscale}\n");
    printf("|-: {");
    for(int i=0;i<11;i++){
        printf("%.4E,",comp_l3_params[i]);
    }
    printf("%.4E}\n",comp_l3_params[11]);
}