#include "magneto/wmm.h"

#include <stddef.h>
#include <math.h>

#include "magneto/model.h"
#include "common_private.h"

#define EPOCH           REAL(2000.0)
#define N_MAX           (12U)
#define M_MAX           (12U)
#define TOTAL_COEFFS    (90U)
#define NUM_MODELS      (1U)
#define INTERVAL        REAL(1.0)

STATIC_ASSERT(N_MAX == M_MAX, n_and_m_max_must_be_same);
STATIC_ASSERT(MAGNETO_CALC_INDEX(1, 0) == 0, first_index_must_be_0);
STATIC_ASSERT(MAGNETO_CALC_INDEX(N_MAX, M_MAX) == (TOTAL_COEFFS - 1), last_index_must_match_total);

static const magneto_SphericalHarmonicCoeff COEFFS_WMM2020[TOTAL_COEFFS] = {
    // Auto-generated table by `tools/gen_coeffs.py`
    { .g = REAL(-2.9404500000000000e+04), .h = REAL( 0.0000000000000000e+00) },  // (n =   1, m =   0)
    { .g = REAL(-1.4507000000000000e+03), .h = REAL( 4.6528999999999996e+03) },  // (n =   1, m =   1)
    { .g = REAL(-1.2500000000000000e+03), .h = REAL( 0.0000000000000000e+00) },  // (n =   2, m =   0)
    { .g = REAL( 1.7216585027234639e+03), .h = REAL(-1.7272010653076843e+03) },  // (n =   2, m =   1)
    { .g = REAL( 4.8405046568858222e+02), .h = REAL(-2.1211848890026849e+02) },  // (n =   2, m =   2)
    { .g = REAL( 6.8195000000000005e+02), .h = REAL( 0.0000000000000000e+00) },  // (n =   3, m =   0)
    { .g = REAL(-1.4580587693916866e+03), .h = REAL(-5.0337014214194305e+01) },  // (n =   3, m =   1)
    { .g = REAL( 4.7877820125816083e+02), .h = REAL( 9.3648737311295335e+01) },  // (n =   3, m =   2)
    { .g = REAL( 8.3120468297525861e+01), .h = REAL(-8.5840027085270663e+01) },  // (n =   3, m =   3)
    { .g = REAL( 1.1288750000000000e+02), .h = REAL( 0.0000000000000000e+00) },  // (n =   4, m =   0)
    { .g = REAL( 1.2797737690701432e+02), .h = REAL( 4.4588115008374153e+01) },  // (n =   4, m =   1)
    { .g = REAL( 9.6374529830240956e+00), .h = REAL(-1.7709658381798338e+01) },  // (n =   4, m =   2)
    { .g = REAL(-1.8490186586403070e+01), .h = REAL( 1.1940333807250592e+01) },  // (n =   4, m =   3)
    { .g = REAL( 1.0120722200373986e+00), .h = REAL(-7.3972126145113428e+00) },  // (n =   4, m =   4)
    { .g = REAL(-2.9300000000000001e+01), .h = REAL( 0.0000000000000000e+00) },  // (n =   5, m =   0)
    { .g = REAL( 5.8595010541996380e+01), .h = REAL( 7.6975544005872409e+00) },  // (n =   5, m =   1)
    { .g = REAL( 2.2909254212466820e+01), .h = REAL( 2.5422196900309295e+01) },  // (n =   5, m =   2)
    { .g = REAL(-1.0510541583334323e+01), .h = REAL(-9.0613268945163714e+00) },  // (n =   5, m =   3)
    { .g = REAL(-5.3244718047896544e+00), .h = REAL( 1.1339152917607600e+00) },  // (n =   5, m =   4)
    { .g = REAL( 1.5256162559167558e-01), .h = REAL( 1.1035662113967188e+00) },  // (n =   5, m =   5)
    { .g = REAL( 1.3729166666666668e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =   6, m =   0)
    { .g = REAL( 1.7893866999351369e+00), .h = REAL(-5.2099521293843176e-01) },  // (n =   6, m =   1)
    { .g = REAL( 1.5742120572497456e+00), .h = REAL( 5.3911371823621423e-01) },  // (n =   6, m =   2)
    { .g = REAL(-1.7467284470853341e+00), .h = REAL( 7.5763447869462641e-01) },  // (n =   6, m =   3)
    { .g = REAL(-2.8504820672779224e-01), .h = REAL(-5.0710233462071330e-01) },  // (n =   6, m =   4)
    { .g = REAL( 4.5327541726427639e-02), .h = REAL( 3.0218361150951761e-02) },  // (n =   6, m =   5)
    { .g = REAL(-6.2710758763313612e-02), .h = REAL( 6.6006223675141515e-02) },  // (n =   6, m =   6)
    { .g = REAL( 1.6791666666666665e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =   7, m =   0)
    { .g = REAL(-2.1166010488516722e+00), .h = REAL(-1.4165793477991662e+00) },  // (n =   7, m =   1)
    { .g = REAL(-1.8677134651661542e-01), .h = REAL(-3.7804320740712521e-01) },  // (n =   7, m =   2)
    { .g = REAL( 8.9901224571182281e-01), .h = REAL( 3.6596958674994552e-02) },  // (n =   7, m =   3)
    { .g = REAL( 1.5160303638031344e-01), .h = REAL( 2.2548552879350417e-01) },  // (n =   7, m =   4)
    { .g = REAL( 3.0704412431455885e-02), .h = REAL(-1.0554641773312961e-02) },  // (n =   7, m =   5)
    { .g = REAL(-1.3548669069928364e-02), .h = REAL(-5.1183860930840486e-02) },  // (n =   7, m =   6)
    { .g = REAL( 4.9286297770162728e-03), .h = REAL(-9.5555067105417522e-04) },  // (n =   7, m =   7)
    { .g = REAL( 6.1458333333333337e-02), .h = REAL( 0.0000000000000000e+00) },  // (n =   8, m =   0)
    { .g = REAL( 3.4027777777777782e-02), .h = REAL( 2.9166666666666667e-02) },  // (n =   8, m =   1)
    { .g = REAL(-5.0838716890091395e-02), .h = REAL(-4.4447563909622768e-02) },  // (n =   8, m =   2)
    { .g = REAL(-8.5821441757406218e-04), .h = REAL( 2.7462861362369990e-02) },  // (n =   8, m =   3)
    { .g = REAL(-2.9222182540084251e-02), .h = REAL(-1.6342263221468919e-02) },  // (n =   8, m =   4)
    { .g = REAL( 1.1753844594950104e-02), .h = REAL( 1.1446554540180167e-02) },  // (n =   8, m =   5)
    { .g = REAL( 4.8719821694714520e-03), .h = REAL( 1.2802288912479728e-03) },  // (n =   8, m =   6)
    { .g = REAL(-2.1425868521027145e-03), .h = REAL(-8.9599086542477150e-04) },  // (n =   8, m =   7)
    { .g = REAL(-9.7390311459214279e-06), .h = REAL( 9.0897624028600001e-05) },  // (n =   8, m =   8)
    { .g = REAL( 1.3020833333333332e-02), .h = REAL( 0.0000000000000000e+00) },  // (n =   9, m =   0)
    { .g = REAL( 2.8649620961716051e-02), .h = REAL(-8.1406849805851722e-02) },  // (n =   9, m =   1)
    { .g = REAL( 8.6407547150381802e-03), .h = REAL( 3.3073233564456482e-02) },  // (n =   9, m =   2)
    { .g = REAL(-3.1859584804880834e-03), .h = REAL( 2.2301709363416588e-02) },  // (n =   9, m =   3)
    { .g = REAL(-1.7006255919061533e-03), .h = REAL(-7.8847186533830740e-03) },  // (n =   9, m =   4)
    { .g = REAL(-1.2288210729148123e-02), .h = REAL(-5.7283388361442375e-03) },  // (n =   9, m =   5)
    { .g = REAL( 5.2482445172476969e-04), .h = REAL( 3.7214824758665483e-03) },  // (n =   9, m =   6)
    { .g = REAL( 1.8387048357799655e-03), .h = REAL( 8.2638419585616418e-05) },  // (n =   9, m =   7)
    { .g = REAL(-6.5901529973659871e-04), .h = REAL(-1.0629279028009655e-04) },  // (n =   9, m =   8)
    { .g = REAL(-1.9875737739993781e-04), .h = REAL( 1.6201231603188207e-04) },  // (n =   9, m =   9)
    { .g = REAL(-4.9479166666666660e-04), .h = REAL( 0.0000000000000000e+00) },  // (n =  10, m =   0)
    { .g = REAL(-2.1771037225375529e-03), .h = REAL( 1.1938955897786579e-03) },  // (n =  10, m =   1)
    { .g = REAL(-3.0410115006309491e-05), .h = REAL(-6.0820230012618981e-05) },  // (n =  10, m =   2)
    { .g = REAL( 4.0554616724186750e-04), .h = REAL( 8.3494799138031546e-04) },  // (n =  10, m =   3)
    { .g = REAL(-1.5181647085108494e-04), .h = REAL( 8.0968784453911964e-04) },  // (n =  10, m =   4)
    { .g = REAL( 6.4011444562398634e-05), .h = REAL(-9.1749737206104709e-04) },  // (n =  10, m =   5)
    { .g = REAL(-5.3675228017305988e-05), .h = REAL(-5.9639142241451100e-06) },  // (n =  10, m =   6)
    { .g = REAL( 5.4965543232612220e-05), .h = REAL(-1.2150277977735333e-04) },  // (n =  10, m =   7)
    { .g = REAL( 1.6534434043570574e-05), .h = REAL(-4.0155054105814248e-05) },  // (n =  10, m =   8)
    { .g = REAL(-9.1962493307687283e-06), .h = REAL(-3.8317705544869705e-07) },  // (n =  10, m =   9)
    { .g = REAL(-3.3415587792658218e-06), .h = REAL(-7.5399275019331374e-06) },  // (n =  10, m =  10)
    { .g = REAL( 7.8125000000000004e-04), .h = REAL( 0.0000000000000000e+00) },  // (n =  11, m =   0)
    { .g = REAL(-4.9364816694836568e-04), .h = REAL( 0.0000000000000000e+00) },  // (n =  11, m =   1)
    { .g = REAL(-7.7313943488978439e-04), .h = REAL( 8.0406501228537581e-04) },  // (n =  11, m =   2)
    { .g = REAL( 5.9509499195552760e-04), .h = REAL(-1.2397812332406825e-04) },  // (n =  11, m =   3)
    { .g = REAL(-1.6297347547619890e-04), .h = REAL(-7.2432655767199522e-05) },  // (n =  11, m =   4)
    { .g = REAL( 3.5932273867492175e-05), .h = REAL( 7.1864547734984351e-05) },  // (n =  11, m =   5)
    { .g = REAL(-4.9809553756847286e-05), .h = REAL(-1.4231301073384939e-05) },  // (n =  11, m =   6)
    { .g = REAL(-3.7502771216246222e-06), .h = REAL(-6.3754711067618577e-05) },  // (n =  11, m =   7)
    { .g = REAL( 2.4090432186113738e-05), .h = REAL(-2.7531922498415701e-05) },  // (n =  11, m =   8)
    { .g = REAL(-3.9986503997038993e-06), .h = REAL(-1.9993251998519501e-05) },  // (n =  11, m =   9)
    { .g = REAL( 4.1133676800104636e-07), .h = REAL(-4.1133676800104633e-06) },  // (n =  11, m =  10)
    { .g = REAL( 1.3593080508189816e-06), .h = REAL(-1.1400648168159201e-06) },  // (n =  11, m =  11)
    { .g = REAL(-4.3402777777777779e-05), .h = REAL( 0.0000000000000000e+00) },  // (n =  12, m =   0)
    { .g = REAL(-2.9486381097515517e-06), .h = REAL(-3.5383657317018615e-05) },  // (n =  12, m =   1)
    { .g = REAL( 1.3068441657910044e-05), .h = REAL( 1.3068441657910044e-05) },  // (n =  12, m =   2)
    { .g = REAL( 2.7742878622516237e-05), .h = REAL( 2.7742878622516237e-05) },  // (n =  12, m =   3)
    { .g = REAL(-1.9206608277126624e-05), .h = REAL(-2.8809912415689940e-05) },  // (n =  12, m =   4)
    { .g = REAL( 7.6857810047984364e-06), .h = REAL( 1.0979687149712054e-06) },  // (n =  12, m =   5)
    { .g = REAL( 2.0541113764093508e-06), .h = REAL( 4.7929265449551517e-06) },  // (n =  12, m =   6)
    { .g = REAL( 1.9238515705535426e-06), .h = REAL(-3.8477031411070853e-07) },  // (n =  12, m =   7)
    { .g = REAL(-3.8477031411070847e-07), .h = REAL( 1.1543109423321252e-06) },  // (n =  12, m =   8)
    { .g = REAL(-4.1981883085339448e-07), .h = REAL( 1.6792753234135781e-07) },  // (n =  12, m =   9)
    { .g = REAL( 3.1005675498566781e-08), .h = REAL(-2.7905107948710101e-07) },  // (n =  12, m =  10)
    { .g = REAL(-1.0057382384307833e-07), .h = REAL( 0.0000000000000000e+00) },  // (n =  12, m =  11)
    { .g = REAL(-5.5989670430932392e-09), .h = REAL( 9.3316117384887328e-09) },  // (n =  12, m =  12)
};

static const magneto_ModelCoeffs SUBMODELS_WMM2020[NUM_MODELS] = {
    { .coeffs = COEFFS_WMM2020 }
};

static const magneto_SphericalHarmonicCoeff SECULAR_WMM2020[TOTAL_COEFFS] = {
    // Auto-generated table by `tools/gen_coeffs.py`
    { .g = REAL( 6.7000000000000002e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =   1, m =   0)
    { .g = REAL( 7.7000000000000002e+00), .h = REAL(-2.5100000000000001e+01) },  // (n =   1, m =   1)
    { .g = REAL(-5.7500000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =   2, m =   0)
    { .g = REAL(-4.0991869112463428e+00), .h = REAL(-1.7435978129526696e+01) },  // (n =   2, m =   1)
    { .g = REAL(-6.3508529610858833e-01), .h = REAL(-6.8993357168160268e+00) },  // (n =   2, m =   2)
    { .g = REAL( 1.3999999999999999e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =   3, m =   0)
    { .g = REAL(-3.7967091013139260e+00), .h = REAL( 3.4905228834660287e+00) },  // (n =   3, m =   1)
    { .g = REAL( 1.3168143377105215e+00), .h = REAL(-3.8729833462074165e-01) },  // (n =   3, m =   2)
    { .g = REAL(-1.9289893727027114e+00), .h = REAL( 1.7392527130926089e-01) },  // (n =   3, m =   3)
    { .g = REAL(-1.3750000000000001e-01), .h = REAL( 0.0000000000000000e+00) },  // (n =   4, m =   0)
    { .g = REAL(-2.5298221281347039e-01), .h = REAL( 3.1622776601683798e-02) },  // (n =   4, m =   1)
    { .g = REAL(-6.7082039324993703e-01), .h = REAL( 7.7144345223742761e-01) },  // (n =   4, m =   2)
    { .g = REAL( 3.2271172452028629e-01), .h = REAL( 2.2111729272686284e-01) },  // (n =   4, m =   3)
    { .g = REAL(-1.1620871002517104e-01), .h = REAL(-1.1832159566199232e-01) },  // (n =   4, m =   4)
    { .g = REAL(-3.7499999999999999e-02), .h = REAL( 0.0000000000000000e+00) },  // (n =   5, m =   0)
    { .g = REAL( 9.6824583655185412e-02), .h = REAL( 1.6137430609197572e-02) },  // (n =   5, m =   1)
    { .g = REAL(-8.5391256382996661e-02), .h = REAL( 3.0496877279641665e-01) },  // (n =   5, m =   2)
    { .g = REAL( 7.4701788083399601e-03), .h = REAL(-6.7231609275059639e-02) },  // (n =   5, m =   3)
    { .g = REAL( 4.2257712736425833e-02), .h = REAL( 1.0564428184106459e-01) },  // (n =   5, m =   4)
    { .g = REAL( 1.1135885079684349e-02), .h = REAL( 5.5679425398421746e-03) },  // (n =   5, m =   5)
    { .g = REAL(-1.2499999999999999e-02), .h = REAL( 0.0000000000000000e+00) },  // (n =   6, m =   0)
    { .g = REAL(-1.0910894511799617e-02), .h = REAL( 2.7277236279499044e-03) },  // (n =   6, m =   1)
    { .g = REAL( 1.0782274364724285e-02), .h = REAL(-3.8816187713007426e-02) },  // (n =   6, m =   2)
    { .g = REAL( 2.0126912147485330e-02), .h = REAL(-2.0126912147485330e-02) },  // (n =   6, m =   3)
    { .g = REAL(-1.1023963796102461e-02), .h = REAL( 7.0868338689230115e-03) },  // (n =   6, m =   4)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 3.3575956834390848e-04) },  // (n =   6, m =   5)
    { .g = REAL( 7.7540350866539239e-04), .h = REAL( 9.6925438583174038e-04) },  // (n =   6, m =   6)
    { .g = REAL(-2.0833333333333333e-03), .h = REAL( 0.0000000000000000e+00) },  // (n =   7, m =   0)
    { .g = REAL(-8.2679728470768446e-03), .h = REAL( 1.3779954745128076e-02) },  // (n =   7, m =   1)
    { .g = REAL(-2.2502571869471738e-03), .h = REAL( 1.3501543121683042e-02) },  // (n =   7, m =   2)
    { .g = REAL( 1.1138204814128777e-02), .h = REAL(-1.1138204814128777e-02) },  // (n =   7, m =   3)
    { .g = REAL( 1.9190257769659928e-03), .h = REAL(-1.9190257769659928e-03) },  // (n =   7, m =   4)
    { .g = REAL(-2.3987822212074910e-03), .h = REAL(-5.7570773308979785e-03) },  // (n =   7, m =   5)
    { .g = REAL(-1.5054076744364850e-03), .h = REAL( 3.7635191860912125e-04) },  // (n =   7, m =   6)
    { .g = REAL( 5.0292140581798699e-04), .h = REAL( 1.5087642174539610e-04) },  // (n =   7, m =   7)
    { .g = REAL(-2.6041666666666666e-04), .h = REAL( 0.0000000000000000e+00) },  // (n =   8, m =   0)
    { .g = REAL( 3.4722222222222224e-04), .h = REAL(-1.0416666666666667e-03) },  // (n =   8, m =   1)
    { .g = REAL(-2.9050695365766515e-04), .h = REAL( 2.0335486756036560e-03) },  // (n =   8, m =   2)
    { .g = REAL( 1.0727680219675777e-03), .h = REAL(-4.2910720878703109e-04) },  // (n =   8, m =   3)
    { .g = REAL(-1.3849375611414336e-04), .h = REAL( 6.9246878057071680e-04) },  // (n =   8, m =   4)
    { .g = REAL( 3.0729005476993739e-04), .h = REAL(-2.3046754107745301e-04) },  // (n =   8, m =   5)
    { .g = REAL( 1.7780956822888511e-04), .h = REAL(-1.7780956822888511e-04) },  // (n =   8, m =   6)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 5.1941499444914289e-05) },  // (n =   8, m =   7)
    { .g = REAL( 1.2985374861228572e-05), .h = REAL( 3.2463437153071431e-06) },  // (n =   8, m =   8)
    { .g = REAL(-2.6041666666666666e-04), .h = REAL( 0.0000000000000000e+00) },  // (n =   9, m =   0)
    { .g = REAL(-6.9877124296868428e-04), .h = REAL(-1.0481568644530263e-03) },  // (n =   9, m =   1)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 5.9591411827849519e-04) },  // (n =   9, m =   2)
    { .g = REAL( 9.1027385156802396e-04), .h = REAL(-9.1027385156802396e-04) },  // (n =   9, m =   3)
    { .g = REAL(-4.6380697961076901e-04), .h = REAL( 6.1840930614769215e-04) },  // (n =   9, m =   4)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 9.2392561873294154e-05) },  // (n =   9, m =   5)
    { .g = REAL( 1.4313394137948263e-04), .h = REAL( 0.0000000000000000e+00) },  // (n =   9, m =   6)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL(-4.1319209792808209e-05) },  // (n =   9, m =   7)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 3.5430930093365516e-05) },  // (n =   9, m =   8)
    { .g = REAL(-6.6809202487374053e-06), .h = REAL( 3.3404601243687027e-06) },  // (n =   9, m =   9)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  10, m =   0)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  10, m =   1)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 3.0410115006309491e-05) },  // (n =  10, m =   2)
    { .g = REAL( 4.7711313793160887e-05), .h = REAL(-7.1566970689741327e-05) },  // (n =  10, m =   3)
    { .g = REAL(-1.6868496761231663e-05), .h = REAL( 1.6868496761231663e-05) },  // (n =  10, m =   4)
    { .g = REAL(-2.1337148187466214e-05), .h = REAL(-2.1337148187466214e-05) },  // (n =  10, m =   5)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 5.9639142241451100e-06) },  // (n =  10, m =   6)
    { .g = REAL(-2.8929233280322222e-06), .h = REAL( 0.0000000000000000e+00) },  // (n =  10, m =   7)
    { .g = REAL(-2.3620620062243678e-06), .h = REAL(-1.1810310031121839e-06) },  // (n =  10, m =   8)
    { .g = REAL(-3.8317705544869705e-07), .h = REAL( 7.6635411089739410e-07) },  // (n =  10, m =   9)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  10, m =  10)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  11, m =   0)
    { .g = REAL(-3.5260583353454692e-05), .h = REAL( 0.0000000000000000e+00) },  // (n =  11, m =   1)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 3.0925577395591381e-05) },  // (n =  11, m =   2)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  11, m =   3)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 3.6216327883599761e-05) },  // (n =  11, m =   4)
    { .g = REAL(-1.1977424622497393e-05), .h = REAL( 0.0000000000000000e+00) },  // (n =  11, m =   5)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  11, m =   6)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 3.7502771216246222e-06) },  // (n =  11, m =   7)
    { .g = REAL(-1.7207451561509813e-06), .h = REAL( 0.0000000000000000e+00) },  // (n =  11, m =   8)
    { .g = REAL(-6.6644173328398332e-07), .h = REAL(-6.6644173328398332e-07) },  // (n =  11, m =   9)
    { .g = REAL(-2.0566838400052318e-07), .h = REAL( 0.0000000000000000e+00) },  // (n =  11, m =  10)
    { .g = REAL(-4.3848646800612310e-08), .h = REAL( 0.0000000000000000e+00) },  // (n =  11, m =  11)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  12, m =   0)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  12, m =   1)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  12, m =   2)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL(-2.1340675863474029e-06) },  // (n =  12, m =   3)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 1.6005506897605523e-06) },  // (n =  12, m =   4)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  12, m =   5)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  12, m =   6)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  12, m =   7)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 1.9238515705535424e-07) },  // (n =  12, m =   8)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  12, m =   9)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  12, m =  10)
    { .g = REAL( 0.0000000000000000e+00), .h = REAL( 0.0000000000000000e+00) },  // (n =  12, m =  11)
    { .g = REAL(-1.8663223476977468e-09), .h = REAL(-1.8663223476977468e-09) },  // (n =  12, m =  12)
};

const magneto_Model MODEL_WMM2020 = {
    .epoch = {
        .year = EPOCH
    },
    .nm_max = N_MAX,
    .num_model_coeffs = TOTAL_COEFFS,
    .num_models = 1U,
    .model_interval = {                 // Unused
        .year = 1
    },
    .models = SUBMODELS_WMM2020,
    .last_secular = {
        .coeffs = SECULAR_WMM2020
    }
};
