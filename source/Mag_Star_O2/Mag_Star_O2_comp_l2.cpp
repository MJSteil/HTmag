//
// Integration of the l=0 O[B**2] structure equations
// Created by m. J. Steil on 2017.03.01.
//

#include "../../include/Mag_Star_O2.hpp"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix.h>

void Mag_Star_O2::comp_l20(int mode, int print, vector<double> comp_params) {

    //region quadrupole equations for n2H0[r], y2H0[r], n2P0[r] and y2P0[r] {dn2H0/dr,dy2H0/dr,dn2P0/dr,dy2P0/dr}
    function<double (double, const double *, double *, int)> dl2_eq = [&](double r , const double * n2y2, double *dn2y2, int i)->double{
        double df_i_dr;

        double rr = r*r;

        double aphi0 = mag_star->aphi0(r);
        double aphi02 = aphi0*aphi0;
        double daphi0 = mag_star->daphi0(r);
        double daphi02 = daphi0*daphi0;
        double cjphi0 = mag_star->cjphi0;

        double P = star->P(r);
        double rho = star->rho(r);

        double expLambda = star->expLambda(r);
        double dNu = star->dNu(r);

        switch (i)
        {
            case 0: // dn2H/dr
                df_i_dr = ( ( 2.*(1.-expLambda)/rr +8.*M_PI*expLambda*(P+rho))/dNu- dNu )*n2y2[0] -4.*expLambda/(rr*dNu)*n2y2[1];
                break;
            case 1: // dy2H/dr
                df_i_dr = -n2y2[0]*dNu;
                break;
            case 2: // dn2P/dr
                df_i_dr = ( (2.+ 2.*expLambda*(-1. + M_4PI*rr*(P+rho)))/(rr*dNu) - dNu )*n2y2[2] -4.*expLambda/(rr*dNu)*n2y2[3]+
                        4./3.*daphi0*aphi0/rr+
                        1./3.*daphi02*dNu/expLambda+
                        16./3.*cjphi0*M_PI*expLambda*aphi0*(P+rho)/dNu;
                break;
            case 3: // dy2P/dr
                df_i_dr = -n2y2[2]*dNu +
                        2./3.*( (-1. + 1/expLambda)/rr + M_4PI*(P+rho) )*aphi0*daphi0 +
                        0.5/expLambda*daphi02*dNu +
                        4./3.*r*M_PI*(P+rho)*cjphi0*( 2.*aphi0 + r*daphi0 );
                break;
            default:
                GSL_ERROR_VAL ("dl2_eq:", GSL_FAILURE,GSL_FAILURE);
                df_i_dr = 0;
        }
        return df_i_dr;

    };

    gsl_odeiv2_ode_system n2y2_eq(4);
    n2y2_eq.df= &dl2_eq;
    struct gsl_odeiv2_ode_struct comp_l2_odes = {&n2y2_eq};
    //endregion

    //region Computational Parameters
    if(comp_params.empty()-1){
        comp_l2_params = comp_params;
    }

    const double comp_l2_delta_r_0 = comp_l2_params[0];    // [km], Initial step size
    const double comp_l2_delta_r_max =  comp_l2_params[1]*star->R; // [km], Maximal step size
    const double comp_l2_rMax = star->R;            // [km], Maximum value for radial integration

    double comp_l2_odeiv2_control_params[4] = {comp_l2_params[2],comp_l2_params[3],comp_l2_params[4],comp_l2_params[5]};  // [error_rel, error_abs, error_ay, error_day]
    double comp_l2_error_scale[4]           = {comp_l2_params[6],comp_l2_params[7],comp_l2_params[8],comp_l2_params[9]};  // Error scale for *comp_l2_control

    unsigned long comp_l2_nMax = (unsigned long)comp_l2_params[10];   // Maximum number of Iteration steps of gsl_odeiv2_evolve_apply
    //endregion

    // gsl_odeiv2_solver setup
    gsl_odeiv2_solver comp_l2_solver(4,mode,print);
    comp_l2_solver.set_parameters(comp_l2_delta_r_0,
                                  comp_l2_delta_r_max,
                                  comp_l2_rMax,
                                  comp_l2_nMax,
                                  (double *)comp_l2_odeiv2_control_params,
                                  (double *)comp_l2_error_scale
    );
    comp_l2_solver.set_system(&comp_l2_odes);

    // gsl_odeiv2_data data container
    gsl_odeiv2_data comp_l2_data(8,comp_l2_nMax);

    double Hscale = 1E-5;
    double n2y2[4] = {gsl_pow_2(comp_l2_delta_r_0)*Hscale, -2./3.*M_PI*(3*star->PC+star->rhoC)*gsl_pow_4(comp_l2_delta_r_0)*Hscale,
                      gsl_pow_2(comp_l2_delta_r_0), -2./9.*M_PI*(3.*(3*star->PC+star->rhoC) + 3.*(-2.+mag_star->cjphi0)*star->PC + (-2.+ 3.*mag_star->cjphi0)*star->rhoC)*gsl_pow_4(comp_l2_delta_r_0)};
    double dn2y2[4] = {2.*comp_l2_delta_r_0*Hscale, -8./3.*M_PI*(3*star->PC+star->rhoC)*gsl_pow_3(comp_l2_delta_r_0)*Hscale,
                       2.*comp_l2_delta_r_0, -8./9.*M_PI*(3.*(3*star->PC+star->rhoC) + 3.*(-2.+mag_star->cjphi0)*star->PC + (-2.+ 3.*mag_star->cjphi0)*star->rhoC)*gsl_pow_3(comp_l2_delta_r_0)};
    double comp_l2_r = comp_l2_delta_r_0;
    double comp_l2_delta_r = comp_l2_delta_r_0;

    unsigned int comp_l2_n=0;

    // Central Values
    comp_l2_data.put(comp_l2_n,{0,comp_l2_delta_r_0,0,0,0,0,0,0,0,0});
    comp_l2_n++;
    comp_l2_data.put(comp_l2_n,{comp_l2_r,comp_l2_delta_r,n2y2[0],dn2y2[0],n2y2[1],dn2y2[1],n2y2[2],dn2y2[2],n2y2[3],dn2y2[3]});
    comp_l2_n++;

    //region Main comp_l2_solver.comp loop
    int comp_l2_status;
    comp_l2_solver.set_IC(comp_l2_n,
                          comp_l2_r,
                          comp_l2_delta_r,
                          n2y2,
                          dn2y2,
                          &comp_l2_data);// Inital conditions for integration

    comp_l2_solver.comp();// Integration
    comp_l2_solver.get_nrh(comp_l2_n,comp_l2_r,comp_l2_delta_r);
    //endregion

    //region Solve LGS: {{n20_ext(R,1,0),-n2y2[0]},{y20_ext(R,1,0),-n2y2[1]}}.{Qm0,cn2H00}=={n2y2[2]-n20_ext(R,0,1),n2y2[3]-y20_ext(R,0,1)};
    double R = star->R;

    double comp_l2_lgs_m_data[] = {n20_ext(R,1,0),-n2y2[0],y20_ext(R,1,0),-n2y2[1]};
    double comp_l2_lgs_b_data[] = {n2y2[2]-n20_ext(R,0,1),n2y2[3]-y20_ext(R,0,1)};

    gsl_matrix_view comp_l2_lgs_m = gsl_matrix_view_array (comp_l2_lgs_m_data, 2, 2);
    gsl_vector_view comp_l2_lgs_b = gsl_vector_view_array (comp_l2_lgs_b_data, 2);

    gsl_vector *comp_l2_lgs_sol = gsl_vector_alloc (2);
    gsl_vector *comp_l2_lgs_norm = gsl_vector_alloc (2);
    gsl_vector *comp_l2_lgs_tau = gsl_vector_alloc (2);
    gsl_permutation * comp_l2_lgs_p = gsl_permutation_alloc (2);
    int comp_l2_lgs_s;

    gsl_linalg_QRPT_decomp  (&comp_l2_lgs_m.matrix, comp_l2_lgs_tau,comp_l2_lgs_p,&comp_l2_lgs_s,comp_l2_lgs_norm);
    gsl_linalg_QRPT_solve  (&comp_l2_lgs_m.matrix, comp_l2_lgs_tau,comp_l2_lgs_p, &comp_l2_lgs_b.vector, comp_l2_lgs_sol);

    Qm0 = gsl_vector_get(comp_l2_lgs_sol,0);
    cn2H00 = gsl_vector_get(comp_l2_lgs_sol,1);
    nl2 = comp_l2_n;

    if(print){
        printf("|-: comp_l2_lgs: {{%.4E,%.4E},{%.4E,%.4E}}.{Qm0,cn2H00}=={%.4E,%.4E} \n",
               comp_l2_lgs_m_data[0], comp_l2_lgs_m_data[1], comp_l2_lgs_m_data[2],
               comp_l2_lgs_m_data[3], comp_l2_lgs_b_data[0], comp_l2_lgs_b_data[1]);
        printf("|-: comp_l2_lgs_sol = {Qm0=%.8E, cn2H00=%.8E}\n",Qm0,cn2H00);
        
        double comp_l2_lgs_c;
        gsl_vector *comp_l2_lgs_work = gsl_vector_alloc (6);
        gsl_linalg_QRPT_rcond (&comp_l2_lgs_m.matrix, &comp_l2_lgs_c, comp_l2_lgs_work);
        printf("|-: comp_l2_lgs_c = %.8E (=1/gsl_linalg_QRPT_rcond)\n",1/ comp_l2_lgs_c);

    };
    //endregion


    for (int i=2;i<=5;i++){
        comp_l2_data.scale(cn2H00,i); // Scale homogeneous solution
        comp_l2_data.add(i+4,i);    // Add homogeneous and particular solution
    }

    //region Paste computation Data to star if(mode==1)
    if(mode==1&&comp_l2_solver.solver_status==5){
        // Assign temporary method data to star
        comp_l2_data.copy(0,&rl2_i);

        comp_l2_data.copy(2,&n2H0_i);
        comp_l2_data.copy(3,&dn2H0_i);
        comp_l2_data.copy(4,&y2H0_i);
        comp_l2_data.copy(5,&dy2H0_i);

        comp_l2_data.copy(6,&n20_i);
        comp_l2_data.copy(7,&dn20_i);
        comp_l2_data.copy(8,&y20_i);
        comp_l2_data.copy(9,&dy20_i);

        mag_star_O2_l2_status = 1;
    }
    //endregion
}

void Mag_Star_O2::comp_l20_print_params() {
    printf("|-: comp_l2_params = {h0, hmax, error_rel, error_abs, error_ay, error_day, error_scale[4], nMax}\n");
    printf("|-: {");
    for(int i=0;i<10;i++){
        printf("%.4E,",comp_l2_params[i]);
    }
    printf("%.4E}\n",comp_l2_params[10]);
}

void Mag_Star_O2::comp_l20_int() {
    if (mag_star_O2_l2_status < 1) {
        GSL_ERROR_VAL ("comp_l20_int: no valid input data - comp_l20 has not been run successfully.", GSL_FAILURE,);
    } else {
        mag_star_O2_l2_status = 2;

        n20_of_r = gsl_hermite_spline_obj(rl2_i,n20_i,dn20_i,nl2);
        n2H0_of_r = gsl_hermite_spline_obj(rl2_i,n2H0_i,dn2H0_i,nl2);

        y20_of_r = gsl_hermite_spline_obj(rl2_i,y20_i,dy20_i,nl2);
        y2H0_of_r = gsl_hermite_spline_obj(rl2_i,y2H0_i,dy2H0_i,nl2);
    }
}
