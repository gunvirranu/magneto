#include "magneto/magneto.h"

#include <stddef.h>
#include <math.h>

#define REAL(x) MAGNETO_REAL(x)

typedef magneto_real real;

const real magneto_PI = REAL(3.141592653589793238462643383279);
static const real RAD_PER_DEG = REAL(0.017453292519943295769);
static const real DEG_PER_RAD = REAL(57.29577951308232087680);

const real magneto_WGS84_A = REAL(6378137.0);
const real magneto_WGS84_F = (real) (1000000000.0 / 298257223563LL);
const real magneto_WGS84_F_INV = (real) (298257223563LL / 1000000000.0);
const real magneto_WGS84_E_SQ = (magneto_WGS84_F * (2 - magneto_WGS84_F));

real magneto_rad_to_deg(const real rad) {
    return rad * DEG_PER_RAD;
}

real magneto_deg_to_rad(const real deg) {
    return deg * RAD_PER_DEG;
}

magneto_FieldState magneto_FieldState_from_ned(const magneto_real *B_ned) {
    magneto_FieldState field = { 0 };
    if (B_ned != NULL) {
        field.B_ned[0] = B_ned[0];
        field.B_ned[1] = B_ned[1];
        field.B_ned[2] = B_ned[2];
        field.H = REAL(hypot)(B_ned[0], B_ned[1]);
        field.F = REAL(hypot)(field.H, B_ned[2]);
        field.D = magneto_rad_to_deg(REAL(atan2)(B_ned[1], B_ned[0]));
        field.I = magneto_rad_to_deg(REAL(atan2)(B_ned[2], field.H));
    }
    return field;
}
