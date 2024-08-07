//
// Integration of the l=0 O[B**2] structure equations
// Created by M. J. Steil on 2017.03.01.
//

#include "../../include/Mag_Star_O2.hpp"

void Mag_Star_O2::comp_l00_odeiv(int mode, int print, vector<double> comp_params) {
    //region Method for m0 and h0 Integration => m0[r], h0[r]
    /*      int mode = 0: Do not copy data to Rot_Star_O1
     *                 1: Default - Copy data to Rot_Star_O1
     *      int print = 0: Do not print text to console
     *                  1: Default - Print text to console
     * */
    //endregion

    //region monopole equations for m0[r], h0[r], MB[r] {dm0/dr,dh0/dr, dMB/dr}
    function<double (double, const double *, double *, int)> dl0_eq = [&](double r , const double * m0h0, double *dm0h0, int i)->double{
        double df_i_dr;

        double rr = r*r;

        double aphi0 = mag_star->aphi0(r);
        double aphi02 = aphi0*aphi0;
        double daphi0 = mag_star->daphi0(r);
        double daphi02 = daphi0*daphi0;
        double cjphi0 = mag_star->cjphi0;

        double P = star->P(r);
        double rho = star->rho(r);
        double drhodh = star->drhodh(r);

        double expLambda = star->expLambda(r);
        double dNu = star->dNu(r);


        switch (i)
        {
            case 0: // dm0/dr
                df_i_dr = ( 2.*aphi02/rr + 1./(expLambda)*daphi02 )/3. + M_4PI*rr*drhodh*m0h0[1];
                break;
            case 1: // dh0/dr
                df_i_dr = expLambda*( 2./3./r*aphi02/rr - M_4PI*r*(P+rho)*m0h0[1] - m0h0[0]*(1.+r*dNu)/rr ) + 2./3.*cjphi0*daphi0 - daphi02/(3.*r) ;
                break;
            case 2: // dMB/dr
                if(r>0){
                    df_i_dr = M_4PI*star->star_eos->mB*sqrt(expLambda)*r*r*( expLambda*m0h0[0]*star->nbar(r)/r + star->dnbardh(r)*m0h0[1]);
                }else{
                    df_i_dr = 0;
                }
                break;
            default:
                GSL_ERROR_VAL ("dl0_eq:", GSL_FAILURE,GSL_FAILURE);
                df_i_dr = 0;
        }
        return df_i_dr;
    };

    gsl_odeiv2_ode_system m0h0_eq(3);
    m0h0_eq.df= &dl0_eq;
    struct gsl_odeiv2_ode_struct comp_l0_ode = {&m0h0_eq};
    //endregion

    // Computational Parameters
    if(comp_params.empty()-1){
        comp_l0_params = comp_params;
    }

    const double comp_l0_delta_r_0 = comp_l0_params[0];   // [km], Initial step size
    const double comp_l0_delta_r_max = comp_l0_params[1]*star->R;   // [km], Maximal step size
    const double comp_l0_rMax = star->R;    // [km], Maximum value for radial integration
    double comp_l0_odeiv2_control_params[4] = {comp_l0_params[2],comp_l0_params[3],comp_l0_params[4],comp_l0_params[5]};  // [error_rel, error_abs, error_ay, error_day]
    double comp_l0_error_scale[3]           = {comp_l0_params[6], comp_l0_params[7], comp_l0_params[8]};       // Error scale for *comp_l0_control

    unsigned long comp_l0_nMax = (unsigned long)comp_l0_params[9];   // Maximum number of Iteration steps of gsl_odeiv2_evolve_apply

    // gsl_odeiv2_solver setup
    gsl_odeiv2_solver comp_l0_solver(3,mode,print);
    comp_l0_solver.set_parameters(comp_l0_delta_r_0,
                                  comp_l0_delta_r_max,
                                  comp_l0_rMax,
                                  comp_l0_nMax,
                                  (double *)comp_l0_odeiv2_control_params,
                                  (double *)comp_l0_error_scale
    );
    comp_l0_solver.set_system(&comp_l0_ode);

    // gsl_odeiv2_data data container
    gsl_odeiv2_data comp_l0_data(6,comp_l0_nMax);

    double Pc = star->PC;
    double rhoC = star->rhoC;
    double drhodhC = star->drhodhC;
    double dnbardhC = star->dnbardhC;
    double cjphi0 = mag_star->cjphi0;

    ch0r0 = comp_l0_params[10];

    double m0h0[3] = {( 1 + 8*M_PI*drhodhC*ch0r0 )*gsl_pow_3(comp_l0_delta_r_0)/6.,
                      ( -4*M_PI*(3*Pc+3*rhoC+drhodhC)*ch0r0 - (1+2*mag_star->cjphi0) )*gsl_pow_2(comp_l0_delta_r_0)/6. + ch0r0,
                      M_4PI*star->star_eos->mB*dnbardhC*ch0r0*gsl_pow_2(comp_l0_delta_r_0)};
    double dm0h0[3] = { m0h0[1]*3/comp_l0_delta_r_0, m0h0[1]*2/comp_l0_delta_r_0, m0h0[2]*2/comp_l0_delta_r_0 };
    double comp_l0_r = comp_l0_delta_r_0;
    double comp_l0_delta_r = comp_l0_delta_r_0;

    unsigned int comp_l0_n=0;

    // Central Values
    comp_l0_data.put(comp_l0_n,{0,comp_l0_delta_r_0,0,0,ch0r0,0,0,0});
    comp_l0_n++;
    comp_l0_data.put(comp_l0_n,{comp_l0_r,comp_l0_delta_r,m0h0[0],dm0h0[0],m0h0[1],dm0h0[1],m0h0[2],dm0h0[2]});
    comp_l0_n++;

    //region Main comp_l0_solver.comp loop
    int comp_l0_status;
    comp_l0_solver.set_IC(comp_l0_n,
                          comp_l0_r,
                          comp_l0_delta_r,
                          m0h0,
                          dm0h0,
                          &comp_l0_data);// Initial conditions for integration
    comp_l0_solver.comp();// Integration
    comp_l0_solver.get_nrh(comp_l0_n,comp_l0_r,comp_l0_delta_r);
    //endregion

    nl0=comp_l0_n;
    m0R = m0h0[0];
    h0R = m0h0[1];
    MBR = m0h0[2];

    if(print){printf("|-: l0 integration successfull@ n=%.d - r=%.16E - h=%.8E - m0=%.8E - h0=%.8E - MB=%.8E\n",comp_l0_n,comp_l0_r,comp_l0_delta_r, m0R,h0R,MBR);};

    if(print){printf("|-: m0(Rs)=%.8E | h0(Rs)=%.8E\n",m0R,h0R);};

    //region Paste computation Data to star if(mode==1)
    if(mode==1&&comp_l0_solver.solver_status==5){
        comp_l0_data.copy(0,&rl0_i);

        //comp_l0_data.scale(1/1,2);
        //comp_l0_data.scale(1/1,3);
        //comp_l0_data.scale(1/1,4);
        //comp_l0_data.scale(1/1,5);

        comp_l0_data.copy(2,&m00_i);
        comp_l0_data.copy(3,&dm00_i);
        comp_l0_data.copy(4,&h00_i);
        comp_l0_data.copy(5,&dh00_i);
        comp_l0_data.copy(6,&MB0_i);
        comp_l0_data.copy(7,&dMB0_i);
    }
    //endregion
}

void Mag_Star_O2::comp_l00_match(int print) {
    deltaM0=0;
    double R = star->R;
    double M = star->M;
    double QQ0 = gsl_pow_2(mag_star->Qe0);

    deltaM0=m0R-m00_ext(R,0) + deltaM0_C*4*M_PI*gsl_pow_3(R)*(2*M-R)/(M-R)*star->rho(R)*(h0R) +QQ0*(R-2*M)/(R-M)/(2*R);

    ch00=0;
    double n00R= -h0R+2./3*mag_star->cjphi0*mag_star->aphi0(R);
    ch00= n00_ext(R,deltaM0)-n00R;

    mag_star_O2_l0_status=1;

    if(print){
        printf("|-: deltaM = %.8E MS | ch0 = %.8E\n",deltaM()/MSkm,ch00*BB());
        printf("|-: deltaMB = %.8E MS | ch0r0 = %.8E\n",MBR*BB()/MSkm,ch0r0*BB());
    }

}

void Mag_Star_O2::comp_l00_MB0(int print, int mode) {
    double MBref = star->MB;
    double hCref = star->hC;
    comp_l0_params[10] = 0;
    comp_l00_odeiv(0,0);

    gsl_function_pp_fkt comp_MB0_fkt = gsl_function_pp_fkt([&](double f0)->double{
        if(mode==0){
            comp_l0_params[10] = f0;
            comp_l00_odeiv(0,0);
            return MBR;
        }else{
            comp_l0_params[10] = 0;

            star->comp(hCref-f0,0);
            mag_star->comp(mag_star->B0,1,0);
            comp_l00_odeiv(0,0);

            return star->MB+MBR*BB()-MBref;
        }
    });

    double Z = star->M/star->R;
    vector<double> comp_cjphi0_range_Z;
    if(mode==0){
        comp_cjphi0_range_Z = {-5E2,0,0};
    }else{
        comp_cjphi0_range_Z = {0,0,hCref*0.8};
    }

    gsl_rootfinder star_target(comp_MB0_fkt.gsl_fkt(), 200, 1e-16,0, comp_cjphi0_range_Z,print);
}

void Mag_Star_O2::comp_l00(int mode, int print, vector<double> comp_params) {
    comp_l00_odeiv(mode,print,comp_params);
    comp_l00_match(print);
}

void Mag_Star_O2::comp_l00_print_params() {
    printf("|-: comp_l0_params = {h0, hmax, error_rel, error_abs, error_ay, error_day, error_scale[3], nMax, ch0r0}\n");
    printf("|-: {");
    for(int i=0;i<9;i++){
        printf("%.4E,",comp_l0_params[i]);
    }
    printf("%.4E}\n",comp_l0_params[10]);
}


void Mag_Star_O2::comp_l00_int() {
    if (mag_star_O2_l0_status < 1) {
        GSL_ERROR_VAL ("comp_l00_int: no valid input data - comp_l00 has not been run successfully.", GSL_FAILURE,);
    } else {
        mag_star_O2_l0_status = 2;

        m00_of_r = gsl_hermite_spline_obj(rl0_i,m00_i,dm00_i,nl0);
        h00_of_r = gsl_hermite_spline_obj(rl0_i,h00_i,dh00_i,nl0);
    }
}

void Mag_Star_O2::comp_deltaB0(int print) {
    if (mag_star_O2_l0_status<2) {
        GSL_ERROR_VAL ("comp_deltaB0: Background star with interpolations required.", GSL_FAILURE,);
    }

    gsl_function_pp_fkt comp_deltaB0_dB = gsl_function_pp_fkt([&](double r)->double{
        double expLambda = star->expLambda(r);

        double P = star->P(r);
        double rho = star->rho(r);
        double nbar = star->nbar(r);
        double dnbardh = star->dnbardh(r);

        if(r>0){
            return sqrt(expLambda)*r*r*( expLambda*m00_of_r.f(r)*nbar/r + dnbardh*h00_of_r.f(r) );
        }else{
            return 0;
        }

    });

    // Integration Method
    gsl_integration comp_deltaB0_int(comp_deltaB0_dB.gsl_fkt(),{0,star->R},1E-14,1E-10,"odeiv2",(int)1E5,6,print==2);

    deltaB0 = M_4PI*comp_deltaB0_int.F;
    deltaMB0 = deltaB0*(star->star_eos->mB);

    if(print==1||print==2){
        printf("|=> Baryon number perturbation deltaB: %.8E -> B = %.8E\n",BB()*deltaB0,BB()*deltaB0+star->B);
        printf("|=> Baryonic mass perturbation deltaMB: %.8E MS -> MB = %.8E MS\n",BB()*deltaMB0/MSkm,(star->MB+BB()*deltaMB0)/MSkm);

        /*vector<double> ri, fi;
        for(double r=0; r<=1.0*star->R; r=r+star->R/1000.){
            ri.push_back(r/star->R);
            fi.push_back(comp_deltaB0_dBf(r));
        }
        Plot test(ri,fi);
        test.plot();*/
    }
}

double Mag_Star_O2::deltaM_EM(int print) {

    gsl_function_pp_fkt test_deltaM0_int_dM_Tem = gsl_function_pp_fkt([&](double r)->double{
        double rr = r*r;

        double R = star->R;
        double expLambda = star->expLambda(r);
        double expNu = star->expNu(r);

        double aphi0 = mag_star->aphi0(r);
        double daphi0 = mag_star->daphi0(r);
        double QQ0 = gsl_pow_2(mag_star->Qe0);

        if(r>0){
            return 2./3./rr*sqrt(expNu/expLambda)*( 2.*expLambda*aphi0*aphi0 + rr*daphi0*daphi0 ) + QQ0/(rr)*(r>R);
        }else{
            return 0;
        }

    });

    //region Integration
    gsl_integration test_deltaM0_int_M_Tem_i(test_deltaM0_int_dM_Tem.gsl_fkt(),{0,star->R},1E-16,1E-10,"odeiv2",(int)1E6,6,print==2);
    double Mem0_i = test_deltaM0_int_M_Tem_i.F;

    double Mem0_e_int = 0;
    if(print==2){
        gsl_integration test_deltaM0_int_M_Tem_e(test_deltaM0_int_dM_Tem.gsl_fkt(),{star->R,INFINITY},1E-16,1E-12,"qagiu",1000,0,print==2);
        Mem0_e_int = test_deltaM0_int_M_Tem_e.F;
    }


    double y = 1-2*star->M/star->R;
    double Mem0_Q = gsl_pow_2(mag_star->Qe0)/star->R;
    double Mem0_e = -3*gsl_pow_2(mag_star->mu0)*(3-4*y+y*y+2*log(y))*(-1+y*y-2*y*log(y))/(8*gsl_pow_3(star->M*(y-1))) + Mem0_Q;
    double Mem0 = Mem0_i+Mem0_e;
    //endregion

    //region Printing
    if(print==1||print==2){
        printf("|=> deltaM_EM = %.8E MS = %.2f deltaM = [%.2E | %.2E (%.2E|%.2E)] = [int, ext(Q,B)]\n",Mem0*BB()/MSkm,Mem0/deltaM0,
               Mem0_i/Mem0,Mem0_e/Mem0,Mem0_Q/Mem0_e,1-Mem0_Q/Mem0_e);
        if(print==2)printf("|-: abs(1-Mem0_e_int/Mem0_e) = %.4E \n",abs(1-Mem0_e_int/Mem0_e));

        /*vector<double> ri, fi;
        for(double r=0; r<=1.0*star->R; r=r+star->R/1000.){
            ri.push_back(r/star->R);
            fi.push_back(test_deltaM0_int_dM_Tem_f(r));
        }
        Plot test(ri,fi);
        test.plot();*/
    }
    //endregion
    return Mem0*BB();
}

double Mag_Star_O2::deltaM_F(int print) {

    gsl_function_pp_fkt test_deltaM0_int_dM_TF = gsl_function_pp_fkt([&](double r)->double{
        double expLambda = star->expLambda(r);
        double expNu = star->expNu(r);
        double detg = M_4PI*r*r*sqrt(expLambda*expNu);

        double P = star->P(r);
        double rho = star->rho(r);
        double drhodh = star->drhodh(r);

        double cjphi0 = mag_star->cjphi0;
        double aphi0 = mag_star->aphi0(r);

        if(r>0) {
            return detg *(  expLambda/r*(3*P+rho)*m00(r) + n00(r)*(3*P+rho) + (drhodh+3*(P+rho))*h00(r) );
//            return detg * (expLambda / r * (3. * P + rho) * m00(r) + (drhodh + 2. * rho) * h00(r) +
//                           (P + rho / 3.) * (3. * ch00 + 2. * cjphi0 * aphi0));
        }else{
            return 0;
        }
    });

    gsl_integration test_deltaM0_int_M_Tem_i(test_deltaM0_int_dM_TF.gsl_fkt(),{0,star->R},1E-16,1E-10,"odeiv2",(int)1E6,6,print==2);
    double Mf0_i = test_deltaM0_int_M_Tem_i.F;

    if(print==1||print==2){
        printf("|=> deltaM_F = %.8E MS = %.2f deltaM \n",Mf0_i*BB()/MSkm,Mf0_i/deltaM0);
    }

    return Mf0_i*BB();
}

double Mag_Star_O2::deltaM_error(int print){
    double dMF  = deltaM_F(print);
    double dMEM = deltaM_EM(print);
    double dM = dMF + dMEM;

    if(print==1||print==2){
        printf("|=> deltaM: %.8E [%.4E | %.4E] = %.4E M_0\n",dM/MSkm,dMF/dM,dMEM/dM,dM/star->M);
        printf("|=> |1-deltaM_int/deltaM| = %.8E \n",abs(1-dM/deltaM()));
    }

    return abs(1-dM/deltaM());
}

double Mag_Star_O2::Bmax_hc(){
    return sqrt(-star->hC/ch0r0);
}