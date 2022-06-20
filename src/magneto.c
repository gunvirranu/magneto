#include "magneto/magneto.h"

#define C(x) MAGNETO_REAL(x)

typedef magneto_real real;

const real magneto_pi = C(3.141592653589793238462643383279);
const real magneto_rad_per_deg = C(0.017453292519943295769);
const real magneto_deg_per_rad = C(57.29577951308232087680);

const real magneto_wgs84_a = C(6378137.0);
const real magneto_wgs84_f = (real) (1000000000.0 / 298257223563LL);
const real magneto_wgs84_f_inv = (real) (298257223563LL / 1000000000.0);
const real magneto_wgs84_e_sq = (magneto_wgs84_f * (2 - magneto_wgs84_f));

real magneto_rad_to_deg(const real rad) {
    return rad * magneto_deg_per_rad;
}

real magneto_deg_to_rad(const real deg) {
    return deg * magneto_rad_per_deg;
}
