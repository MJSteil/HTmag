//
// Created by M. J. Steil on 2017.05.12.
//

#ifndef HTMAG_ISS_O11_SERIES_HPP
#define HTMAG_ISS_O11_SERIES_HPP

#include <gsl/gsl_poly.h>

#include "Static_Star.hpp"
#include "Static_Star_ISS.hpp"


class ISS_O11_Series {
public:
    ISS_O11_Series(Static_Star* iss_in, int print = 1);
    ISS_O11_Series(Static_Star_ISS* iss_in, int print = 1);


    Static_Star *iss;
    Static_Star_ISS *iss_ana;

    double *Z, *M, *R;

    double J0(int order);
    double omegabar0(double r, int d,int order);

    double mu0(int order);
    double B0(int order);

    double aphi0(double r, int d,int order);
    double Btet0(double r, double theta, int i, int order);

    double Qe00(int order);

private:

    double wI[21] = {0.,1.,8.57142857142857E-1,1.00952380952381,1.367965367965368,2.00868845440274,3.107779839208411,4.988344958496219,8.2281714784026,1.386022335083852E1,2.373917569650862E1,4.121151809358846E1,7.234385639676935E1,1.281807689943514E2,2.289060239186798E2,4.115294351727705E2,7.441150622479728E2,1.352162461548796E3,2.467591308474046E3,4.519799724118707E3,8.30507932491463E3};

    double wi[21][41] = {
            {1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-2.,0.,1.2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-3.E-1,0.,-1.8,0.,1.414285714285714,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-2.857142857142857E-1,0.,-3.6E-1,0.,-2.057142857142857,0.,1.895238095238095,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-3.41547619047619E-1,0.,-3.728571428571428E-1,0.,-3.921428571428571E-1,0.,-2.75,0.,2.762175324675325,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-4.553246753246753E-1,0.,-4.912857142857143E-1,0.,-3.783673469387755E-1,0.,-5.328571428571428E-1,0.,-4.002597402597402,0.,4.253481518481519,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-6.451596379810666E-1,0.,-7.341753246753246E-1,0.,-4.527397959183673E-1,0.,-5.274489795918368E-1,0.,-7.757305194805195E-1,0.,-6.160077422577422,0.,6.809107808857809,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-9.47038818324533E-1,0.,-1.188314942200656,0.,-5.956530612244898E-1,0.,-6.599501133786848E-1,0.,-7.650371057513915E-1,0.,-1.195162837162837,0.,-9.85808158508159,0.,1.12185624963272E1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-1.416352535763723,0.,-2.034717689453404,0.,-8.15588482435931E-1,0.,-9.33350597814884E-1,0.,-9.4851557668522E-1,0.,-1.18002294134437,0.,-1.913500907425908,0.,-1.623884229426036E1,0.,1.889835384246172E1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-2.129575635848088,0.,-3.634687250137818,0.,-1.117558800892985,0.,-1.430036060538101,0.,-1.317318346966886,0.,-1.46664483016983,0.,-1.889364593739594,0.,-3.152975556376397,0.,-2.73518535043434E1,0.,3.240183589834229E1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-3.176539830884248,0.,-6.71205610565476,0.,-1.464532733015034,0.,-2.328212131349376,0.,-1.954852898357826,0.,-2.047354089903603,0.,-2.34776422924298,0.,-3.113608664654673,0.,-5.311727348752806,0.,-4.689154091647152E1,0.,5.635684839107994E1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-4.621854850193863,0.,-1.273122124771019E1,0.,-1.68012700620415,0.,-3.990119835401471,0.,-3.021699062652715,0.,-3.068447726496067,0.,-3.274539265970321,0.,-3.869745013953633,0.,-5.24575902972387,0.,-9.10754737699736,0.,-8.15538353412232E1,0.,9.91956812816561E1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-6.382807925060304,0.,-2.468434560766691E1,0.,-1.198269390990551,0.,-7.163960349273301,0.,-4.782776799017114,0.,-4.829573021615355,0.,-4.897093715652909,0.,-5.399056795919848,0.,-6.52006183879652,0.,-8.99487127852501,0.,-1.584138030442162E1,0.,-1.435391216301188E2,0.,1.763582335396428E2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-7.89779662822539,0.,-4.874037276213144E1,0.,1.561475917574373,0.,-1.343525269165858E1,0.,-7.641263107869459,0.,-7.888886574456563,0.,-7.672378782694945,0.,-8.07948022339275,0.,-9.09704166333395,0.,-1.118046687011353E1,0.,-1.564592568201857E1,0.,-2.788367067142283E1,0.,-2.551866116269776E2,0.,3.1624305617124E2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {-7.291646265790306,0.,-9.77188191010702E1,0.,1.071546208230642E1,0.,-2.624729542554509E1,0.,-1.213429978382535E1,0.,-1.328982503117874E1,0.,-1.242142973292298E1,0.,-1.267474139818235E1,0.,-1.361280551895957E1,0.,-1.560027621878028E1,0.,-1.944826122077808E1,0.,-2.754032518308143E1,0.,-4.957483287942578E1,0.,-4.575843446949453E2,0.,5.712986212372352E2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {6.343533609599749E-1,0.,-1.984442578019848E2,0.,3.655897478221529E1,0.,-5.323631859872981E1,0.,-1.875044570967779E1,0.,-2.302876846659568E1,0.,-2.059157058370308E1,0.,-2.057322744426698E1,0.,-2.135082024145026E1,0.,-2.334586069132159E1,0.,-2.713733551265366E1,0.,-3.423413706320704E1,0.,-4.896532232993282E1,0.,-8.88982188643302E1,0.,-8.26616005798712E2,0.,1.038755412825174E3,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {2.910573004249907E1,0.,-4.073852297693708E2,0.,1.041811097771345E2,0.,-1.116258576613659E2,0.,-2.715204084203302E1,0.,-4.104673336444683E1,0.,-3.471097912227812E1,0.,-3.427437225991082E1,0.,-3.463700749157261E1,0.,-3.662036364999298E1,0.,-4.06122930677183E1,0.,-4.777015521126612E1,0.,-6.086772764637196E1,0.,-8.78065297276755E1,0.,-1.605981229413583E2,0.,-1.502955908621476E3,0.,1.899484431758825E3,0.,0.,0.,0.,0.,0.,0.,0.},
            {1.101655636597088E2,0.,-8.44037791928959E2,0.,2.736504802914801E2,0.,-2.40774063755725E2,0.,-3.363417925107924E1,0.,-7.541070966763625E1,0.,-5.913408332056123E1,0.,-5.830564914401429E1,0.,-5.7633612410269E1,0.,-5.941968344874764E1,0.,-6.370586179673208E1,0.,-7.149208611524662E1,0.,-8.49361280706807E1,0.,-1.091519753237292E2,0.,-1.58627792590549E2,0.,-2.920080552150657E2,0.,-2.74828884142286E3,0.,3.491014500271629E3,0.,0.,0.,0.,0.,0.},
            {3.190836935565262E2,0.,-1.76242878046851E3,0.,6.869835668640713E2,0.,-5.314976054735777E2,0.,-2.367055286508275E1,0.,-1.432845665862976E2,0.,-1.011597314441518E2,0.,-1.009597025999027E2,0.,-9.77986510115404E1,0.,-9.89064903395325E1,0.,-1.033683526260625E2,0.,-1.121479767229753E2,0.,-1.271164574764684E2,0.,-1.523152038977315E2,0.,-1.971918933769873E2,0.,-2.884282025336116E2,0.,-5.339747995420792E2,0.,-5.05094987852669E3,0.,6.445058538291366E3,0.,0.,0.,0.},
            {8.28995898431942E2,0.,-3.704722311734869E3,0.,1.676840910085849E3,0.,-1.194816211952608E3,0.,4.872805409927331E1,0.,-2.826099639157569E2,0.,-1.723884849214054E2,0.,-1.776509593447438E2,0.,-1.68535692546115E2,0.,-1.679563756969864E2,0.,-1.720533870107035E2,0.,-1.819754925195075E2,0.,-1.994082120549339E2,0.,-2.27960095255598E2,0.,-2.751729091997684E2,0.,-3.585514733868661E2,0.,-5.274326164508901E2,0.,-9.81385376483359E2,0.,-9.3248936192837E3,0.,1.194710853984579E4,0.,0.},
            {2.032244411837068E3,0.,-7.832069504562398E3,0.,4.017003738113203E3,0.,-2.72341071268953E3,0.,3.161836093336329E2,0.,-5.802369857274802E2,0.,-2.8928842830835E2,0.,-3.175852680244064E2,0.,-2.939743867046293E2,0.,-2.89847921553776E2,0.,-2.921291590296861E2,0.,-3.029047199350916E2,0.,-3.235726423086523E2,0.,-3.576072407252659E2,0.,-4.118377717015768E2,0.,-5.003479407018101E2,0.,-6.556680848493793E2,0.,-9.69367734114809E2,0.,-1.811829146245683E3,0.,-1.728526785293723E4,0.,2.222745028090414E4}
    };


    double aMu[21] = {0.,1.,-8.57142857142857E-1,-3.904761904761905E-1,-3.175757575757576E-1,-3.226316540602255E-1,-3.505530659816374E-1,-3.628401462216588E-1,-2.839488421437243E-1,6.538777728122702E-2,1.106935412373765,3.809714536257503,1.03588909441164E1,2.561635687915638E1,6.029564267689408E1,1.378390294013194E2,3.092738585128257E2,6.852403592558777E2,1.504932037792311E3,3.284320162614748E3};

    double ai[21][43] = {
        {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,5.E-1,0.,-3.E-1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,2.25E-1,0.,1.5E-1,0.,-2.464285714285714E-1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,4.192857142857143E-1,0.,-3.15E-1,0.,3.75E-1,0.,-3.345238095238095E-1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,8.42738095238095E-1,0.,-7.534285714285715E-1,0.,1.671428571428572E-1,0.,4.285714285714286E-1,0.,-4.771103896103896E-1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,1.754861471861472,0.,-1.780142857142857,0.,4.866734693877551E-1,0.,-4.023809523809524E-2,0.,6.461038961038961E-1,0.,-7.287712287712288E-1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,3.726015104538319,0.,-4.182100974025974,0.,1.378366071428571,0.,-2.049030612244898E-1,0.,4.596996753246754E-2,0.,9.89089035964036E-1,0.,-1.159663773726274,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,8.00515197005376,0.,-9.78983152169259,0.,3.783535969387755,0.,-7.399386337868481E-1,0.,1.07994086270872E-1,0.,3.791052697302697E-2,0.,1.585306880619381,0.,-1.90277342143151,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,1.733138683660994E1,0.,-2.285517041101755E1,0.,1.012304438352464E1,0.,-2.39350226242012,0.,3.565953927025356E-1,0.,4.767375481661196E-3,0.,7.452703546453547E-2,0.,2.612012803373097,0.,-3.195491942751521,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,3.772244511296129E1,0.,-5.323834682125367E1,0.,2.654503928698852E1,0.,-7.300629613023031,0.,1.244654701467364,0.,-1.193275177203749E-1,0.,7.176997110032825E-2,0.,1.226217348827643E-1,0.,4.400957225340449,0.,-5.465797512762006,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,8.24193078961191E1,0.,-1.23772167641869E2,0.,6.850690497985275E1,0.,-2.141395487474544E1,0.,4.218222617895787,0.,-5.875951560614709E-1,0.,9.46577596510632E-2,0.,1.016489650055826E-1,0.,2.118196508164823E-1,0.,7.54665681768889,0.,-9.4888178308128,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,1.805946047388233E2,0.,-2.872580488397057E2,0.,1.745447717052429E2,0.,-6.10150305256272E1,0.,1.371911805396405E1,0.,-2.270518425336705,0.,2.382511604331815E-1,0.,5.232492170819376E-2,0.,1.842388258751697E-1,0.,3.68660108636526E-1,0.,1.31276230272695E1,0.,-1.667601061864358E1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,3.965935915776434E2,0.,-6.656492372036294E2,0.,4.400457128850396E2,0.,-1.699363669735181E2,0.,4.299124902453651E1,0.,-8.05346510109996,0.,9.02885056291201E-1,0.,-1.054244049065907E-1,0.,1.364050843055294E-1,0.,3.220445191694639E-1,0.,6.494615340417845E-1,0.,2.310886584864649E1,0.,-2.961029057897748E1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,8.72477945793537E2,0.,-1.540294666360098E3,0.,1.099664984512042E3,0.,-4.646172569212653E2,0.,1.305828328956677E2,0.,-2.719021092820728E1,0.,3.603342118840186,0.,-6.400234480851426E-1,0.,3.569106801750441E-2,0.,2.314376533907982E-1,0.,5.728973940025726E-1,0.,1.155002877359193,0.,4.108864453821876E1,0.,-5.303984833252715E1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,1.922156608413855E3,0.,-3.559588852323329E3,0.,2.727564598321399E3,0.,-1.2508664210235E3,0.,3.86470363565329E2,0.,-8.85356865469985E1,0.,1.372147204104778E1,0.,-2.49434938429458,0.,-1.336881760971547E-1,0.,-9.75447988178117E-4,0.,4.244617626633679E-1,0.,1.025966247926282,0.,2.071220691477107,0.,7.36855098772736E1,0.,-9.57296905436191E1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,4.239780337860013E3,0.,-8.21635925115682E3,0.,6.721985464181E3,0.,-3.323955603521819E3,0.,1.119144561847916E3,0.,-2.799254176560451E2,0.,4.953020019860208E1,0.,-8.98279935123767,0.,-2.118153064565842E-1,0.,-5.962827058298746E-1,0.,5.01677859666157E-2,0.,7.712383240727565E-1,0.,1.850986958893996,0.,3.741194486230681,0.,1.33123584994524E2,0.,-1.739214009335561E2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,9.36135255681745E3,0.,-1.8944498512476E4,0.,1.647374200614495E4,0.,-8.7342683176088E3,0.,3.181452958217075E3,0.,-8.63242450475536E2,0.,1.709242238529431E2,0.,-3.159007898425207E1,0.,8.68604764982026E-1,0.,-2.113603469987957,0.,-8.55655210313521E-1,0.,1.141484829468358E-1,0.,1.4109928389745,0.,3.360542099621568,0.,6.800730763330733,0.,2.42065161430144E2,0.,-3.178157720722704E2,0.,0.,0.,0.,0.,0.,0.,0.},
        {0.,0.,2.068774191726166E4,0.,-4.363620874511033E4,0.,4.017490345388215E4,0.,-2.272805709909474E4,0.,8.90156359249803E3,0.,-2.605324931141741E3,0.,5.68692955863066E2,0.,-1.093161016718608E2,0.,8.00449755204644,0.,-6.103488361675224,0.,-2.900560007218337,0.,-1.516567545376303,0.,2.623667796232514E-1,0.,2.591198453073327,0.,6.135776124229666,0.,1.243209823695649E1,0.,4.426687731584063E2,0.,-5.837504595157506E2,0.,0.,0.,0.,0.,0.},
        {0.,0.,4.575289747802215E4,0.,-1.00416112602326E5,0.,9.75503291282309E4,0.,-5.863818267843718E4,0.,2.456497622121566E4,0.,-7.716324159006982E3,0.,1.836223528052196E3,0.,-3.717551498410959E2,0.,4.076973453722323E1,0.,-1.725788651536651E1,0.,-7.194170054924813,0.,-5.352596656432817,0.,-2.599667246011841,0.,5.576011742177508E-1,0.,4.777827678619766,0.,1.125958714572772E1,0.,2.284072507190477E1,0.,8.13610264176749E2,0.,-1.077130051078365E3,0.,0.,0.,0.},
        {0.,0.,1.012550444193055E5,0.,-2.30877340442027E5,0.,2.359482845798257E5,0.,-1.501440413383696E5,0.,6.697511153625731E4,0.,-2.247769742932245E4,0.,5.781732757404288E3,0.,-1.240504774983364E3,0.,1.717830440113628E2,0.,-5.073429426178545E1,0.,-1.535043072794446E1,0.,-1.408363762648639E1,0.,-9.30829667161525,0.,-4.547581643198817,0.,1.148881421028046,0.,8.84228792297518,0.,2.075616807730496E1,0.,4.215271928392089E1,0.,1.502142320462474E3,0.,-1.995697989426721E3,0.,0.},
        {0.,0.,2.242209230892802E5,0.,-5.304041276897598E5,0.,5.687089283724897E5,0.,-3.818576693146662E5,0.,1.806612822815443E5,0.,-6.452083956049345E4,0.,1.78184762885547E4,0.,-4.060090544001445E3,0.,6.576733169057277E2,0.,-1.573625076527357E2,0.,-2.791160133292036E1,0.,-3.372471650171839E1,0.,-2.435753879517773E1,0.,-1.658342654465128E1,0.,-8.01484118781483,0.,2.31727600898699,0.,1.642094242622291E1,0.,3.841956470787614E1,0.,7.810848379345432E1,0.,2.784611494954315E3,0.,-3.711364623940379E3}
    };

    double aEi[21][20] = {
            {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {2.E-1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {3.E-1,-1.714285714285714E-1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {4.8E-1,-2.571428571428571E-1,-7.809523809523809E-2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {8.E-1,-4.114285714285714E-1,-1.171428571428571E-1,-6.351515151515152E-2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {1.371428571428571,-6.857142857142858E-1,-1.874285714285714E-1,-9.52727272727273E-2,-6.45263308120451E-2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {2.4,-1.175510204081633,-3.123809523809524E-1,-1.524363636363637E-1,-9.67894962180677E-2,-7.011061319632748E-2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {4.266666666666667,-2.057142857142857,-5.355102040816327E-1,-2.54060606060606E-1,-1.548631939489082E-1,-1.051659197944912E-1,-7.256802924433177E-2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {7.68,-3.657142857142857,-9.37142857142857E-1,-4.355324675324676E-1,-2.581053232481804E-1,-1.68265471671186E-1,-1.088520438664976E-1,-5.678976842874485E-2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {1.396363636363636E1,-6.582857142857143,-1.666031746031746,-7.621818181818182E-1,-4.424662684254521E-1,-2.804424527853099E-1,-1.741632701863962E-1,-8.51846526431173E-2,1.30775554562454E-2,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {2.56E1,-1.196883116883117E1,-2.998857142857143,-1.354989898989899,-7.743159697445412E-1,-4.807584904891027E-1,-2.902721169773271E-1,-1.362954442289877E-1,1.961633318436811E-2,2.213870824747529E-1,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {4.726153846153847E1,-2.194285714285714E1,-5.452467532467533,-2.438981818181818,-1.376561723990295,-8.4132735835593E-1,-4.976093433897035E-1,-2.271590737149794E-1,3.138613309498897E-2,3.320806237121294E-1,7.619429072515006E-1,0.,0.,0.,0.,0.,0.,0.,0.,0.},
            {8.77714285714286E1,-4.050989010989011E1,-9.99619047619048,-4.434512396694215,-2.477811103182532,-1.495693081521653,-8.70816350931981E-1,-3.894155549399647E-1,5.231022182498162E-2,5.313289979394071E-1,1.142914360877251,2.071778188823281,0.,0.,0.,0.,0.,0.,0.,0.},
            {1.6384E2,-7.52326530612245E1,-1.845450549450549E1,-8.12993939393939,-4.505111096695512,-2.692247546738975,-1.548117957212411,-6.814772211449384E-1,8.96746659856828E-2,8.85548329899012E-1,1.828662977403601,3.107667283234921,5.123271375831275,0.,0.,0.,0.,0.,0.,0.},
            {3.072E2,-1.404342857142857E2,-3.427265306122449E1,-1.500911888111888E1,-8.25937034394177,-4.89499553952541,-2.78661232298234,-1.211515059813224,1.569306654749449E-1,1.518082851255449,3.047771629006002,4.972267653175874,7.684907063746913,1.205912853537882E1,0.,0.,0.,0.,0.,0.},
            {5.782588235294118E2,-2.633142857142857E2,-6.397561904761905E1,-2.787407792207792E1,-1.524806832727712E1,-8.97415848912992,-5.066567859967891,-2.180727107663803,2.789878497332353E-1,2.656644989697035,5.22475136401029,8.28711275529312,1.229585130199506E1,1.808869280306822E1,2.756780588026387E1,0.,0.,0.,0.,0.},
            {1.092266666666667E3,-4.956504201680672E2,-1.199542857142857E2,-5.203161212121212E1,-2.831784117922894E1,-1.656767721070139E1,-9.28870774327447,-3.96495837757055,5.021781295198235E-1,4.722924426128063,9.14331488701801,1.420647900907392E1,2.04930855033251E1,2.894190848490916E1,4.13517088203958E1,6.185477170256513E1,0.,0.,0.,0.},
            {2.069557894736842E3,-9.36228571428572E2,-2.257963025210084E2,-9.75592727272727E1,-5.285997020122734E1,-3.076854339130258E1,-1.714838352604517E1,-7.269090358879343,9.13051144581498E-1,8.50126396703051,1.625478202136535E1,2.486133826587937E1,3.513100371998589E1,4.823651414151526E1,6.616273411263329E1,9.27821575538477E1,1.370480718511755E2,0.,0.,0.},
            {3.93216E3,-1.773906766917293E3,-4.265041269841269E2,-1.836409839572192E2,-9.91124441273013E1,-5.743461433043147E1,-3.184699797694102E1,-1.341985912408494E1,1.673927098399412,1.545684357641912E1,2.925860763845762E1,4.419793469489665E1,6.14792565099753E1,8.26911670997405E1,1.102712235210555E2,1.484514520861563E2,2.055721077767633E2,3.009864075584623E2,0.,0.},
            {7.489828571428572E3,-3.370422857142857E3,-8.08113082706767E2,-3.468774141414142E2,-1.865646007102142E2,-1.07689901869559E2,-5.944772955695658E1,-2.492259551615775E1,3.090326950891222,2.833754655676838E1,5.319746843355931E1,7.955628245081398E1,1.092964560177339E2,1.447095424245458E2,1.890363831789523E2,2.474190868102605E2,3.289153724428213E2,4.514796113376934E2,6.568640325229497E2,0.},
            {1.429876363636363E4,-6.419853061224489E3,-1.535414857142857E3,-6.572414162679425E2,-3.523998013415157E2,-2.027104035191699E2,-1.114644929192936E2,-4.652217829682779E1,5.739178623083697,5.231547056634162E1,9.75286921281921E1,1.446477862742072E2,1.96733620831921E2,2.57261408754748E2,3.308136705631664E2,4.241470059604466E2,5.481922874047021E2,7.223673781403095E2,9.85296048784424E2,1.426914772106178E3}
    };

    double aB0[21] = {0.,-1.,-4.5E-1,-8.38571428571429E-1,-1.685476190476191,-3.509722943722944,-7.452030209076637,-1.601030394010751E1,-3.466277367321987E1,-7.544489022592257E1,-1.648386157922382E2,-3.611892094776466E2,-7.931871831552868E2,-1.744955891587074E3,-3.844313216827711E3,-8.47956067572003E3,-1.87227051136349E4,-4.137548383452333E4,-9.15057949560443E4,-2.02510088838611E5,-4.484418461785604E5};

    double aBthetaI[21][41] = {
            {0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {1.,0,-1.2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {4.5E-1,0,-4.E-1,0,-2.785714285714286E-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {8.38571428571429E-1,0,-1.71,0,1.15,0,-5.976190476190476E-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {1.68547619047619,0,-3.852285714285714,0,2.037857142857143,0,3.785714285714286E-1,0,-7.556277056277056E-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {3.509722943722943,0,-8.80604761904762,0,5.514469387755102,0,-9.19761904761905E-1,0,9.82467532467532E-1,0,-1.146769896769897,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {7.452030209076637,0,-2.023812683982684E1,0,1.454802976190476E1,0,-3.47169387755102,0,6.289258658008658E-1,0,1.318743756243756,0,-1.792283757908758,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {1.601030394010751E1,0,-4.661135629584701E1,0,3.767475824056896E1,0,-1.147215787981859E1,0,2.241894944341373,0,4.850045787545787E-2,0,2.136672702297702,0,-2.909533542437954,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {3.466277367321987E1,0,-1.074309855841777E2,0,9.61715772833799E1,0,-3.523989343949701E1,0,7.857227878272521,0,-5.133190916226631E-1,0,4.046440018315018E-1,0,3.471247686137392,0,-4.844974836180332,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {7.544489022592257E1,0,-2.476161609582346E2,0,2.426857653959476E2,0,-1.032896552664852E2,0,2.641458231409644E1,0,-3.1976623426177,0,9.38562392964179E-1,0,5.10571738310709E-1,0,5.813005971389291,0,-8.23283280187305,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {1.648386157922382E2,0,-5.705335607933985E2,0,6.066634303275213E2,0,-2.928766858679132E2,0,8.51402740952805E1,0,-1.389004288658672E1,0,2.783356151536955,0,7.484551808485632E-1,0,8.70304791071287E-1,0,9.91357884568991,0,-1.421690496987238E1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {3.611892094776466E2,0,-1.313870811151061E3,0,1.504634855685972E3,0,-8.10016367278237E2,0,2.645716025361063E2,0,-5.26410096101907E1,0,9.57997582139632,0,7.77704906102884E-1,0,1.490987353735661,0,1.434623180370326,0,1.717276862526817E1,0,-2.487675614371453E1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {7.931871831552868E2,0,-3.023786158292164E3,0,3.70688716477294E3,0,-2.196937675848825E3,0,7.973244796075694E2,0,-1.84683150034874E2,0,3.452022149302405E1,0,-8.9883528051531E-1,0,2.804869755336038,0,2.37364858451892,0,2.441902163454604,0,3.012627726416763E1,0,-4.401115588743693E1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {1.744955891587074E3,0,-6.954365848595678E3,0,9.07999225114795E3,0,-5.865115542897067E3,0,2.34207622852163E3,0,-6.148896741318351E2,0,1.23584275297087E2,0,-1.096480026062145E1,0,6.106397942141186,0,4.029155597515737,0,4.061512012561742,0,4.215552495180105,0,5.341456018838135E1,0,-7.859240398350761E1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {3.844313216827711E3,0,-1.598331130088039E4,0,2.212997266379114E4,0,-1.545421740559182E4,0,6.732996515177554E3,0,-1.968729274585325E3,0,4.314254464468213E2,0,-5.59544685772885E1,0,1.594860962300682E1,0,6.801638276084532,0,7.068063874239184,0,6.98949106252792,0,7.377905444990323,0,9.55621492539151E1,0,-1.414725651867717E2,0,0,0,0,0,0,0,0,0,0,0,0},
            {8.47956067572003E3,0,-3.670975022145501E4,0,5.369779024858579E4,0,-4.027303667696039E4,0,1.900493725161476E4,0,-6.11156542521483E3,0,1.462320619190501E3,0,-2.342171261534358E2,0,4.920094389974961E1,0,1.04024560595329E1,0,1.291773863001914E1,0,1.208601859963437E1,0,1.221728952255708E1,0,1.305551044926965E1,0,1.722967001422791E2,0,-2.564352103712587E2,0,0,0,0,0,0,0,0,0,0},
            {1.87227051136349E4,0,-8.4257554725624E4,0,1.297857324330832E5,0,-1.039593595671033E5,0,5.280832795862922E4,0,-1.849780175891158E4,0,4.81704242006492E3,0,-8.93043534280898E2,0,1.684105761172769E2,0,9.9519478106029,0,2.486931155739502E1,0,2.169871117333618E1,0,2.113346187745226E1,0,2.158804235288944E1,0,2.331990966920128E1,0,3.127470307168432E2,0,-4.67652024092682E2,0,0,0,0,0,0,0,0},
            {4.137548383452333E4,0,-1.932675400940763E5,0,3.125876344353369E5,0,-2.661563469357278E5,0,1.447524063457124E5,0,-5.480860225610919E4,0,1.546032146349691E4,0,-3.213824669284519E3,0,6.011588417621611E2,0,-2.115203833316579E1,0,5.145057034827451E1,0,3.979037309202081E1,0,3.809237945887453E1,0,3.730608265927345E1,0,3.851340387665994E1,0,4.198973154420738E1,0,5.71052849436636E2,0,-8.57431545669761E2,0,0,0,0,0,0},
            {9.15057949560443E4,0,-4.430399342438275E5,0,7.504854571930092E5,0,-6.765056654636984E5,0,3.92083015728276E5,0,-1.594682732295023E5,0,4.848119210009733E4,0,-1.109525519872613E4,0,2.154928804267723E3,0,-1.995811101729676E2,0,1.185174935852696E2,0,7.310112876838348E1,0,7.074558835777228E1,0,6.717448057989662E1,0,6.651071669618541E1,0,6.926969701861935E1,0,7.613387613758225E1,0,1.048166458316252E3,0,-1.579612229804849E3,0,0,0,0},
            {2.025100888386109E5,0,-1.015015156724152E6,0,1.796666415970996E6,0,-1.708543240542939E6,0,1.050921138031051E6,0,-4.56714004922914E5,0,1.489302647044337E5,0,-3.707985868714501E4,0,7.62138888109653E3,0,-9.89466396701462E2,0,3.135305902053473E2,0,1.303615205772837E2,0,1.345427227477245E2,0,1.244468383261387E2,0,1.197603780007387E2,0,1.195462053386126E2,0,1.254704592303654E2,0,1.388844440203358E2,0,1.93289217655383E3,0,-2.922531705355494E3,0,0},
            {4.484418461785604E5,0,-2.32412659959765E6,0,4.290010034525023E6,0,-4.290406578708892E6,0,2.790684892931906E6,0,-1.290031539046457E6,0,4.492134482408548E5,0,-1.206457924667385E5,0,2.64176325974639E4,0,-4.13232800131984E3,0,9.49677824563752E2,0,2.095537711094498E2,0,2.62108635393612E2,0,2.34931605112424E2,0,2.220464589290684E2,0,2.152064946932382E2,0,2.164076899608008E2,0,2.286766101290344E2,0,2.547179291708982E2,0,3.579311158799398E3,0,-5.428018949510553E3}
    };

    double aBthetaE[21][20] = {
            {0.,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {-2.E-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {-4.E-1,1.714285714285714E-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {-7.4E-1,3.428571428571429E-1,7.809523809523809E-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {-1.36,6.342857142857143E-1,1.561904761904762E-1,6.351515151515152E-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {-2.512142857142857,1.165714285714286,2.889523809523809E-1,1.27030303030303E-1,6.45263308120451E-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {-4.672857142857143,2.153265306122449,5.310476190476191E-1,2.350060606060606E-1,1.290526616240902E-1,7.011061319632748E-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {-8.75059523809524,4.00530612244898,9.80931972789116E-1,4.319030303030303E-1,2.387474240045669E-1,1.40221226392655E-1,7.256802924433177E-2,0,0,0,0,0,0,0,0,0,0,0,0,0},
            {-1.64847619047619E1,7.500510204081634,1.824639455782313,7.977956709956709E-1,4.387790495219066E-1,2.594092688264117E-1,1.451360584886635E-1,5.678976842874485E-2,0,0,0,0,0,0,0,0,0,0,0,0},
            {-3.121586715367966E1,1.412979591836735E1,3.416899092970521,1.483986147186147,8.10496805235581E-1,4.767521697350269E-1,2.685017082040275E-1,1.135795368574897E-1,-1.30775554562454E-2,0,0,0,0,0,0,0,0,0,0,0},
            {-5.93758143939394E1,2.675645756029685E1,6.436907029478458,2.778976911976912,1.50761162918714,8.80639380755299E-1,4.934625988614561E-1,2.10122143186356E-1,-2.615511091249081E-2,-2.21387082474753E-1,0,0,0,0,0,0,0,0,0,0},
            {-1.133768615238928E2,5.089355519480519E1,1.218905288857968E1,5.23516075036075,2.823219015678199,1.638084398322766,9.1150628161541E-1,3.86170425315465E-1,-4.838695518810799E-2,-4.427741649495058E-1,-7.619429072515005E-1,0,0,0,0,0,0,0,0,0},
            {-2.172207906676657E2,9.71801670204795E1,2.318484181096681E1,9.91340265971402,5.318506000122325,3.067547989878602,1.695500168987208,7.133200555853416E-1,-8.89273771024687E-2,-8.19132205156586E-1,-1.523885814503001,-2.071778188823281,0,0,0,0,0,0,0,0},
            {-4.174051966392982E2,1.861892491437134E2,4.427096497599623E1,1.885631923783287E1,1.007122685271593E1,5.778783827691583,3.175067255717028,1.326852375217317,-1.642634376414825E-1,-1.50543216082832,-2.819188756830552,-4.143556377646561,-5.123271375831276,0,0,0,0,0,0,0},
            {-8.04157087365759E2,3.577758828336842E2,8.48195468321361E1,3.600574269001201E1,1.915651720908961E1,1.09428179379979E1,5.981333419953041,2.484721385925531,-3.055477421241338E-1,-2.780779889513237,-5.181211769310204,-7.665579298646139,-1.024654275166256E1,-1.205912853537882E1,0,0,0,0,0,0},
            {-1.552822557521752E3,6.89277503456365E2,1.629867910686783E2,6.898405715748898E1,3.657896436561066E1,2.081437378095209E1,1.1326369802477E1,4.680829055872118,-5.721819725067371E-1,-5.172551048392263,-9.5705471600126,-1.408809168399831E1,-1.895610409057572E1,-2.411825707075763E1,-2.756780588026387E1,0,0,0,0,0},
            {-3.004635056628421E3,1.330990763590073E3,3.140041960190108E2,1.325577715387832E2,7.008230298937898E1,3.974460641857615E1,2.154392917672705E1,8.86370933479965,-1.077901939962627,-9.68634374839686,-1.780225178299756E1,-2.602301389318385E1,-3.483824535565267E1,-4.461877558090162E1,-5.513561176052774E1,-6.185477170256513E1,0,0,0,0},
            {-5.82448854265073E3,2.575401477110075E3,6.063402367465891E2,2.553807962300957E2,1.346681290050705E2,7.614741416350564E1,4.113767701348203E1,1.685969374849977E1,-2.041136169085175,-1.824756671693095E1,-3.333726987947711E1,-4.840561754029252E1,-6.435194795999498E1,-8.2002074040576E1,-1.020008817569763E2,-1.237095434051303E2,-1.370480718511756E2,0,0,0},
            {-1.130942387989262E4,4.992418750843483E3,1.173238450683479E3,4.931388000856959E2,2.594465312210682E2,1.463226714385742E2,7.881642344824016E1,3.219322855559872E1,-3.882452527482383,-3.455394878036305E1,-6.280223705531536E1,-9.06464617675329E1,-1.197015762167436E2,-1.514712680675975E2,-1.874610799857944E2,-2.28862655299491E2,-2.740961437023511E2,-3.009864075584623E2,0,0},
            {-2.199232702877449E4,9.69379189705082E3,2.274324097606475E3,9.5419925434745E2,5.009897101952726E2,2.818997325069303E2,1.514513625822832E2,6.167959199962799E1,-7.41346097016882,-6.572519159118341E1,-1.189235428572565E2,-1.707638508111532E2,-2.24158370524094E2,-2.817529245658865E2,-3.462713331460288E2,-4.206124475774429E2,-5.070778658493495E2,-6.019728151169245E2,-6.568640325229497E2,0},
            {-4.282431147867957E4,1.885056602466385E4,4.416060753100927E3,1.849716361423626E3,9.69390378167367E2,5.443467084646977E2,2.917804751649753E2,1.185217222905023E2,-1.420358468102936E1,-1.255008629645935E2,-2.262049031987184E2,-3.233617635709929E2,-4.222795440203026E2,-5.276227636863213E2,-6.441020931024509E2,-7.769401145640056E2,-9.31926888587994E2,-1.113649707966311E3,-1.313728065045899E3,-1.426914772106178E3}
    };

    //double QeInfi[21] = {0.,-6.666666666666666E-2,6.380952380952381E-2,-4.444444444444445E-3,4.273757988043702E-3,7.13539884152129E-3,7.81043597370128E-3,4.179852863207348E-3,-1.055972955055271E-2,-5.283514508961481E-2,-1.604526587720938E-1,-4.185939439768732E-1,-1.016630272550158,-2.371625206178944,-5.39581751814005,-1.207433852729258E1,-2.671035302838972E1,-5.860436429715177E1,-1.278125368228011E2,-2.775111453193934E2,-6.005283661346088E2};
};


#endif //HTMAG_ISS_O11_SERIES_HPP