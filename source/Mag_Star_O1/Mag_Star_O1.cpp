//
// Created by M. J. Steil on 2016.07.04.
//

#include "../../include/Mag_Star_O1.hpp"

Mag_Star_O1::Mag_Star_O1() {

}

Mag_Star_O1::Mag_Star_O1(Static_Star * star_in) {
    mag_star_O1_status=0;
    star = star_in;
}

Mag_Star_O1::Mag_Star_O1(Static_Star * star_in, Rot_Star_O1 * rotstar_in) {
    mag_star_O1_status=0;
    star = star_in;
    set_rot(rotstar_in);
}


Mag_Star_O1::Mag_Star_O1(Static_Star * star_in,double Br00, double Qe,int mode, int print, vector<double> comp_cjphi0_range_in, vector<double> comp_params) {

    mag_star_O1_status=0;
    star = star_in;

    set_B0Q0(Br00,Qe);
    if(magnetic_l1){
        comp(Br00, mode, print,comp_cjphi0_range_in, comp_params);
    }
}


void Mag_Star_O1::set_B0Q0(double Br00, double Qe) {
    if(Br00!=0.){
        B0 = Br00;
        eA = B0;
        magnetic_l1 = 1;
    }

    if(Qe!=0.){
        if(magnetic_l1 == 0){
            eA = Qe;
        }
        Qe0 = Qe/eA;
        electric_l0 = 1;
    }
}


void Mag_Star_O1::set_rot(Rot_Star_O1 *rot_star_in) {

    double M = star->M;
    double Z = star->Z;
    double M3 = M*M*M;

    double R = star->R;
    double Y = 1.-2.*M/R;
    double logY = log(Y);

    rot_star = rot_star_in;

    if(magnetic_l1){
        electric_l2  = 1;
//        if(Z<0.01){
            QEr00 = ( M3*mu0*(0.6 + 0.4*logY + (-0.8 + 0.2*Y)*Y)
                      + J0()*mu0*(-0.6 - 0.1*logY + Y*(-4.1 - 2.7*logY + Y*(3.9 - 3.3*logY + (0.9 + 0.1*logY - 0.1*Y)*Y)))
                    )/(M*(-1. + Y*(-9. - 6.*logY + Y*(9. - 6.*logY + Y))));
     /*   }else {
            QEr00 = (J0() * mu0 * (113190 + Z * (16170 + Z * (6468 + Z * (-3696 + (-20944 - 54080 * Z) * Z)))) +
                     mu0 * gsl_pow_3(R) *
                     (113190 + Z * (-56595 + Z * (-3234 + Z * (6468 + Z * (19712 + 43120 * Z)))))) / (339570. * R);
        }*/
    }

    Catoi = (at00_ext(R)-at00(R))*Omega();
}

double Mag_Star_O1::Omega() {
    if(electric_l2){
        return rot_star->Omega;
    }else{
        return 0;
    }
}
double Mag_Star_O1::J0() {
    if(electric_l2){
        return rot_star->J0;
    }else{
        return 0;
    }
}

void Mag_Star_O1::info() {
    double R = star->R;

    printf("Mag_Star_O1:\n");
    printf("|---------------- Magnetic ---------------\n");
    printf("|-> Brtet(r=0,theta=0)= %.8E T = %.8E Gs\n",Btet(0,0,1)/cTkm,Btet(0,0,1)/cGskm);
    printf("|=> mu0 = %.8E km**3 | cjphi0=%.8E| in %d steps \n",mu0, cjphi0,n);
    printf("|=> mu = %.8E km**2 = %.8E Msol**2 = %8E  Am**2 \n",mu(),mu()/gsl_pow_2(MSkm),mu()/cAm2km2);
    printf("|=> cjphi = %.8E A/m = %.8E km**-1 = %.8E \n",cjphi()/cA*1e-3,cjphi(),-cjphi()/(R/star->rhoC/star->rhoC)/cA*1e-6);
    printf("|=> Series values: mu =  %.8E *(-cjphi*R**3) | B0 = %.8E*(-cjphi()) \n",mu()/(-cjphi()*R*R*R), B0/-cjphi());
    printf("|=> Brtet(r=Rs,theta=0)= %.8E T = %.8E Gs\n",Btet(star->R,0,1)/cTkm,Btet(star->R,0,1)/cGskm);
    printf("|=> Bthtet(r=Rs,theta=Pi/2)= %.8E T = %.8E Gs\n",Btet(star->R,M_PI/2,2)/cTkm,Btet(star->R,M_PI/2,2)/cGskm);
    printf("|=> Jphitet(r=Rs/2,theta=Pi/2)= %.8E A m**-2 \n",Jtet(star->R*0.5,M_PI/2,3)/cA*1e-6);
    printf("\n");

    if(electric_l0){
        printf("|---------------- Electric l=0 ---------------\n");
        printf("|-> Q= %.8E C\n",Qe()/cCkm);
        printf("\n");
    }

    if(electric_l2) {
        printf("|---------------- Electric l=2 ---------------\n");
        printf("|-> QeInf= %.8E Cm**2 (QeInf00= %.8E km**5)\n", QEr() / (cCkm) * 1E6,QEr00);
        printf("|-> Q= %.8E C = %.8E + %.8E C \n", Qe() / cCkm, QeI() / cCkm, QeS() / cCkm);
        printf("\n");
        printf("|-> Ertet(r=0,theta=0)= %.8E V/m\n", Etet(0, 0, 1) / (cV) * 1E-3);
        printf("|-> Ertet(r=R,theta=0)= %.8E V/m\n", Etet(R, 0, 1) / (cV) * 1E-3);
        printf("|-> Ertet(r=R+epsilon,theta=0)= %.8E V/m\n", Etet(R, 0, 1, 1) / (cV) * 1E-3);
        printf("|-> Ertet(r=2R,theta=0)= %.8E V/m\n", Etet(2 * R, 0, 1) / (cV) * 1E-3);

        printf("\n");
        printf("|-> Ertet(r=0,theta=M_PI/2)= %.8E V/m\n", Etet(0, M_PI / 2., 1) / (cV) * 1E-3);
        printf("|-> Ertet(r=R,theta=M_PI/2)= %.8E V/m\n", Etet(R, M_PI / 2., 1) / (cV) * 1E-3);
        printf("|-> Ertet(r=R+epsilon,theta=M_PI/4)= %.8E V/m\n", Etet(R, M_PI / 2., 1, 1) / (cV) * 1E-3);
        printf("|-> Ertet(r=2R,theta=M_PI/2)= %.8E V/m\n", Etet(2 * R, M_PI / 2., 1) / (cV) * 1E-3);

        printf("\n");
        printf("|-> Ertet(r=0,theta=M_PI/4)= %.8E V/m\n", Etet(0, M_PI / 4., 1) / (cV) * 1E-3);
        printf("|-> Ertet(r=R,theta=M_PI/4)= %.8E V/m\n", Etet(R, M_PI / 4., 1) / (cV) * 1E-3);
        printf("|-> Ertet(r=R+epsilon,theta=M_PI/4)= %.8E V/m\n", Etet(R, M_PI / 4., 1, 1) / (cV) * 1E-3);
        printf("|-> Ertet(r=2R,theta=M_PI/4)= %.8E V/m\n", Etet(2 * R, M_PI / 4., 1) / (cV) * 1E-3);

        printf("\n");

        printf("|-> Ethtet(r=0,theta=M_PI/4)= %.8E V/m\n", Etet(0, M_PI / 4., 2) / (cV) * 1E-3);
        printf("|-> Ethtet(r=R,theta=M_PI/4)= %.8E V/m\n", Etet(R, M_PI / 4., 2) / (cV) * 1E-3);
        printf("|-> Ethtet(r=R+epsilon,theta=M_PI/4)= %.8E V/m\n", Etet(R, M_PI / 4., 2, 1) / (cV) * 1E-3);
        printf("|-> Ethtet(r=2R,theta=M_PI/4)= %.8E V/m\n", Etet(2 * R, M_PI / 4., 2) / (cV) * 1E-3);
        printf("\n");
    }
}

void Mag_Star_O1::info_E() {
    double R = star->R;
    vector <double> ri, Er0, Er1, Er2, Et0, Et1, Et2, at0i;

    for(double r=0; r<2*R; r=r+R/1.E2){
        ri.push_back(r/R);

        Er0.push_back(Etet(r,0,1)/(cV)*1E-3);
        Er1.push_back(Etet(r,M_PI/4,1)/(cV)*1E-3);
        Er2.push_back(Etet(r,M_PI/2,1)/(cV)*1E-3);

        Et0.push_back(Etet(r,0,2)/(cV)*1E-3);
        Et1.push_back(Etet(r,M_PI/4,2)/(cV)*1E-3);
        Et2.push_back(Etet(r,M_PI/2,2)/(cV)*1E-3);
        at0i.push_back(Atet_2D_plot(r,0,0));

    }

 /*   Plot p1({ri,Er0,ri,Er1,ri,Er2},3);
    p1.plot();

    Plot p2({ri,Et0,ri,Et1,ri,Et2},3);
    p2.plot();

    Plot p3(ri,at0i);
    p3.plot();*/
}


