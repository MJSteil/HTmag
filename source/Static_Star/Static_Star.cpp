//
// Created by martin on 6/23/16.
//

#include "../../include/Static_Star.hpp"

using namespace std;
using namespace gsl_wrapper;
using namespace units;

Static_Star::Static_Star() {
}

Static_Star::Static_Star(eos &eos_in, string name) {
    eos *typeEoS = &eos_in;
    star_eos = typeEoS;
    static_star_status = 0;
    this->name=name;
}

Static_Star::Static_Star(eos &eos_in, double hC_in, int print, vector<double> comp_params) {
    eos *typeEoS = &eos_in;
    star_eos = typeEoS;
    static_star_status = 0;
    comp(hC_in,print,comp_params);
}


//region Minimization method on Static_Star: for calculating stars of specific target Mass, Radius or Compactness as well as maximum Mass finder
double Static_Star::target_fkt(double hC, double target,int type) {

    comp_tov(hC,1,0);

    double f;

    switch ( type )
    {
        case 0: // Search for maximum Mass
            f = -M;
            break;
        case 1: // Search for specific Mass
            f = target-M;
            break;
        case 2: // Search for specific Radius
            f = target-R;
            break;
        case 3: // Search for specific Compactness
            f = target-M/R;
            break;
        case 4: // Search for specific baryon rest mass
            comp_tov_int();
            comp_B(0);
            f = target-MB;
            break;
        case 5: // Search  for minimum Mass
            f = M;
            break;
        case 6: // Search  for minimum Radius
            f = R;
            break;
        case 7: // Search  for maximum Radius
            f = -R;
            break;
        default:
            GSL_ERROR("Static_Star::target_fkt: Invalid type!",GSL_EDOM);
            f = 0;
    }
    return f;
}

void Static_Star::target(int type, double target, vector<double> range, int nMax, double abs_error_threshold, double rel_error_threshold, int print){
     if (type==3&&target>0.5){
        cout << target << endl;
        GSL_ERROR_VAL ("star_target: Target compactness bigger then 4/9 requested, which exceeds Buchdahl’s limit.", GSL_ESANITY,);
    }
    if (type==2&&target<=8){
        cerr << "star_target: Very small target radius of "<< target << " km requested." <<endl;
    }

    // gsl_function via template wrapper gsl_function_pp
    auto star_target_f = [&](double hC_in)->double{
        return  target_fkt(hC_in, target, type);
    };
    gsl_function_pp<decltype(star_target_f)> star_target_fp(star_target_f);
    gsl_function *star_target_ff = static_cast<gsl_function*>(&star_target_fp);

    double x_min;

    if(type==0||type>=5){
        gsl_minimizer star_min(star_target_ff, nMax, abs_error_threshold,rel_error_threshold, range, print==2);
        x_min = star_min.x_min;
    }else{
        gsl_rootfinder star_target(star_target_ff, nMax, abs_error_threshold,rel_error_threshold, range, print==2);
        x_min = star_target.x_min;
    }

    comp(x_min,print==2||print==1);
}

double Static_Star::target_hC(int type, double target_int, vector<double> range, int nMax, double abs_error_threshold, double rel_error_threshold, int print){
    target(type,target_int,range,nMax,abs_error_threshold,rel_error_threshold,print);
    return hC;
}
//endregion


void Static_Star::comp(double hC_in, int print, vector<double> comp_params) {
    //region Method for full BG star computation
    /*      double Pc = Central Pressure in MeVfm**-3
     *      int print = 0: Do not print final BG star info to console
     *                  1: Default - Print final BG star info to console
     *                  2: Print computational information
     * */
    //endregion

    //if(hC_in!=hC){
        comp_tov(hC_in,1,print==2,comp_params);
        comp_tov_int();
        comp_B(print==2);
    //}

    if(print==2||print==1){
        info();
    }
}

void Static_Star::info() {
    printf("Static star: %s\n",name.c_str());
    printf("|-> EoS: %s \n",star_eos->name.c_str());
    printf("|-> P0=%.12E MeV/fm**3 <=> rho0=%.12E MeV/fm**3 <=> nbar0=%.8E fm**-3 <=> logh0=%.17g <=> cs=%.4E \n",PC/cMeVfm3km2,rhoC/cMeVfm3km2,nbarC/cfm3km3,hC,star_eos->cs(hC));
    printf("|-: %d nodes (%d steps, f/s=%.2g)\n",n,count,1.*fails/count);
    printf("|=> RS=%.17g km | MsG/Rs=%.17g | MsG=%.17g Msol | MsB=%.17g Msol\n",R,M/R,M/MSkm,MB/MSkm);
    printf("\n");
}

void Static_Star::save(string filename) {
    FILE *out;
    out=fopen(filename.c_str(),"w+");

    fprintf(out,"#Static_Star \n");
    fprintf(out,"#|-> EoS:%s \n",star_eos->name.c_str());
    fprintf(out,"#|-> P0=%.17g MeV/fm**3 <=> rho0=%.17g MeV/fm**3 <=> logh0=%.17g \n",PC/cMeVfm3km2,rhoC/cMeVfm3km2,hC);
    fprintf(out,"#|=> RS=%.17g km | MsG/Rs=%.17g | MsG=%.17g Msol | MsB=%.17g Msol | %.d points\n",R,M/R,M/MSkm,MB/MSkm,n);
    fprintf(out,"# logh[1], r[km], dlogh/dr [km**-1], m[km], dm/dr[1], P[km**-2], dP/dr[km**-3] \n");

    for(int j=0;j<n;j++){
        fprintf(out,"%.16E, %.16E, %.16E, %.16E, %.16E, %.16E, %.16E\n",h_i[j],r_i[j],dhdr_i[j],m(r_i[j]),dm(r_i[j]),P(r_i[j]),dP(r_i[j]));
    }
    int j=n;
    fprintf(out,"%.16E, %.16E, %.16E, %.16E, %.16E, %.16E, %.16E",h_i[j],r_i[j],dhdr_i[j],m(r_i[j]),dm(r_i[j]),P(r_i[j]),dP(r_i[j]));
    fclose(out);
}