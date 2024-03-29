#include "magneto/magneto.h"

#include <stdbool.h>
#include <stddef.h>
#include <math.h>

#include "common_private.h"

const real magneto_PI = REAL(3.141592653589793238462643383279);
const real magneto_WGS84_A = REAL(6378137.0);
const real magneto_WGS84_B = REAL(6356752.314245);
const real magneto_WGS84_F = (real) (1000000000.0 / 298257223563LL);
const real magneto_WGS84_F_INV = (real) (298257223563LL / 1000000000.0);
const real magneto_WGS84_E_SQ = (magneto_WGS84_F * (2 - magneto_WGS84_F));

static const real RAD_PER_DEG = REAL(0.017453292519943295769);
static const real DEG_PER_RAD = REAL(57.29577951308232087680);
static const int DAYS_IN_YEAR = 365;

static bool is_leap_year(const int year) {
    return ((year % 4) == 0) && ((year % 100) != 0) && ((year % 400) == 0);
}

static real calc_day_of_year(
    const int year, const int month, const int day,
    const int hour, const int minute, const int sec
) {
    // Cumulative number of days in year at start of month
    const int month_days[] = { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 };
    // Check for invalid month only b/c it's used for indexing
    if ((month < 1) || (month > 12)) {
        return 0;
    }
    // Month is now guaranteed to be 1 - 12
    real doy = month_days[month - 1];
    // Add an extra day for Feb. in a leap year
    if ((month > 2) && is_leap_year(year)) {
        doy += 1;
    }
    // Add day and time simply with no checks
    doy += day;
    doy += (hour / (real) 24);
    doy += (minute / (real) 1440);
    doy += (sec / (real) 86400);
    return doy;
}

real magneto_rad_to_deg(const real rad) {
    return rad * DEG_PER_RAD;
}

real magneto_deg_to_rad(const real deg) {
    return deg * RAD_PER_DEG;
}

DecYear magneto_DecYear_from_datetime(
    int_fast16_t year, int_fast8_t month, int_fast8_t day,
    int_fast8_t hour, int_fast8_t minute, int_fast8_t sec
) {
    const real doy = calc_day_of_year(year, month, day, hour, minute, sec);
    const int days_in_year = (DAYS_IN_YEAR + (is_leap_year(year) ? 1 : 0));
    const real year_frac = (doy - 1) / days_in_year;
    return (DecYear) { .year = ((real) year) + year_frac };
}

DecYear magneto_DecYear_from_date(
    int_fast16_t year, int_fast8_t month, int_fast8_t day
) {
    return magneto_DecYear_from_datetime(year, month, day, 0, 0, 0);
}

Coords magneto_Coords_from_spherical(const SphericalCoords pos) {
    const EcefPosition ecef = magneto_EcefPosition_from_spherical(pos);
    return magneto_Coords_from_ecef(ecef);
}

Coords magneto_Coords_from_ecef(const EcefPosition pos) {
    Coords coords = { 0 };

    // Following algorithm adapted from equations found in
    // "Fundamentals of Spacecraft Attitude Determination and Control"
    // by Markley & Crassidis
    const real eps_sq = sq(magneto_WGS84_A / magneto_WGS84_B) - 1;
    if (eps_sq == REAL(0.0)) {
        return coords;
    }
    const real rho_sq = sq(pos.x) + sq(pos.y);
    const real p = REAL(fabs)(pos.z) / eps_sq;
    const real s = (rho_sq / (magneto_WGS84_E_SQ * eps_sq));
    const real q = (sq(p) - sq(magneto_WGS84_B) + s);
    if (q == REAL(0.0)) {
        return coords;
    }
    const real u = (p / SQRT(q));
    const real v = sq(magneto_WGS84_B * u) / q;
    const real P = (27 * v * s / q);
    // TODO: Use `sq(cbrt(...))` for better precision
    const real Q = REAL(pow)(SQRT(P + 1) + SQRT(P), (REAL(2.0) / 3));
    if (Q == REAL(0.0)) {
        return coords;
    }
    const real t = (1 + Q + (1/Q)) / 6;
    const real c = SQRT(sq(u) - 1 + (2 * t));
    const real w = (c - u) / 2;
    real d = SQRT(q) * (w + SQRT(SQRT(sq(t) + v) - ((u*w) + (t/2) + REAL(0.25))));
    d = (pos.z < 0) ? -d : d;
    const real N = (magneto_WGS84_A * SQRT(1 + (eps_sq * sq(d / magneto_WGS84_B))));
    if (N == REAL(0.0)) {
        return coords;
    }
    const real phi = REAL(asin)((eps_sq + 1) * d / N);
    const real h = (SQRT(rho_sq) * COS(phi)) + (pos.z * SIN(phi)) - (sq(magneto_WGS84_A) / N);
    const real lambda = REAL(atan2)(pos.y, pos.x);

    coords.latitude = rad_to_deg(phi);
    coords.longitude = rad_to_deg(lambda);
    coords.height = h;
    return coords;
}

SphericalCoords magneto_SphericalCoords_from_coords(const Coords pos) {
    SphericalCoords spherical = { 0 };

    const Coords pos_with_no_lon = {
        .latitude = pos.latitude,
        .longitude = 0,
        .height = pos.height
    };
    const EcefPosition ecef_no_lon = magneto_EcefPosition_from_coords(pos_with_no_lon);
    const SphericalCoords partial_no_azi = magneto_SphericalCoords_from_ecef(ecef_no_lon);

    spherical.radius = partial_no_azi.radius;
    spherical.polar = partial_no_azi.polar;
    spherical.azimuth = pos.longitude;
    return spherical;
}

SphericalCoords magneto_SphericalCoords_from_ecef(const EcefPosition pos) {
    SphericalCoords spherical = { 0 };

    const real r = HYPOT(HYPOT(pos.x, pos.y), pos.z);
    if (r == REAL(0.0)) {
        return spherical;
    }
    spherical.radius = r;
    spherical.polar = rad_to_deg(REAL(asin)(pos.z / r));
    spherical.azimuth = rad_to_deg(REAL(atan2)(pos.y, pos.x));
    return spherical;
}

EcefPosition magneto_EcefPosition_from_coords(const Coords pos) {
    EcefPosition ecef = { 0 };

    const real lat = deg_to_rad(pos.latitude);
    const real lon = deg_to_rad(pos.longitude);
    const real sin_lat = SIN(lat);
    const real cos_lat = COS(lat);
    const real sin_lon = SIN(lon);
    const real cos_lon = COS(lon);

    const real denom = SQRT(1 - magneto_WGS84_E_SQ * sq(sin_lat));
    if (denom == REAL(0.0)) {
        return ecef;
    }
    const real R_c = (magneto_WGS84_A / denom);  // Radius of curvature

    ecef.x = (R_c + pos.height) * cos_lat * cos_lon;
    ecef.y = (R_c + pos.height) * cos_lat * sin_lon;
    ecef.z = (R_c * (1 - magneto_WGS84_E_SQ) + pos.height) * sin_lat;
    return ecef;
}

EcefPosition magneto_EcefPosition_from_spherical(const SphericalCoords pos) {
    EcefPosition ecef = { 0 };

    const real pol = deg_to_rad(pos.polar);     // Geocentric latitude
    const real azi = deg_to_rad(pos.azimuth);   // Longitude
    const real sin_pol = SIN(pol);
    const real cos_pol = COS(pol);
    const real sin_azi = SIN(azi);
    const real cos_azi = COS(azi);

    ecef.x = (pos.radius * cos_pol * cos_azi);
    ecef.y = (pos.radius * cos_pol * sin_azi);
    ecef.z = (pos.radius * sin_pol);
    return ecef;
}

FieldState magneto_FieldState_from_ned(const real *const B_ned) {
    FieldState field = { 0 };
    if (B_ned == NULL) {
        return field;
    }
    field.B_ned[0] = B_ned[0];
    field.B_ned[1] = B_ned[1];
    field.B_ned[2] = B_ned[2];
    field.H = HYPOT(B_ned[0], B_ned[1]);
    field.F = HYPOT(field.H, B_ned[2]);
    field.D = rad_to_deg(REAL(atan2)(B_ned[1], B_ned[0]));
    field.I = rad_to_deg(REAL(atan2)(B_ned[2], field.H));
    return field;
}

void magneto_convert_vector_ned_to_ecef(
    const magneto_Coords pos,
    const magneto_real *const ned,
    magneto_real *const ecef
) {
    if ((ned == NULL) || (ecef == NULL)) {
        return;
    }

    magneto_real A[9] = { 0 };
    magneto_matrix_ned_to_ecef(pos, A);
    // Standard matrix multiplication
    ecef[0] = (A[0] * ned[0]) + (A[1] * ned[1]) + (A[2] * ned[2]);
    ecef[1] = (A[3] * ned[0]) + (A[4] * ned[1]) + (A[5] * ned[2]);
    ecef[2] = (A[6] * ned[0]) + (A[7] * ned[1]) + (A[8] * ned[2]);
}

void magneto_convert_vector_ecef_to_ned(
    const magneto_Coords pos,
    const magneto_real *const ecef,
    magneto_real *const ned
) {
    if ((ecef == NULL) || (ned == NULL)) {
        return;
    }

    magneto_real A_T[9] = { 0 };
    magneto_matrix_ned_to_ecef(pos, A_T);
    // Standard matrix multiplication
    ned[0] = (A_T[0] * ecef[0]) + (A_T[3] * ecef[1]) + (A_T[6] * ecef[2]);
    ned[1] = (A_T[1] * ecef[0]) + (A_T[4] * ecef[1]) + (A_T[7] * ecef[2]);
    ned[2] = (A_T[2] * ecef[0]) + (A_T[5] * ecef[1]) + (A_T[8] * ecef[2]);
}

void magneto_matrix_ned_to_ecef(const Coords pos, real *const matrix) {
    if (matrix == NULL) {
        return;
    }

    const real lat = deg_to_rad(pos.latitude);
    const real lon = deg_to_rad(pos.longitude);
    const real sin_lat = SIN(lat);
    const real cos_lat = COS(lat);
    const real sin_lon = SIN(lon);
    const real cos_lon = COS(lon);

    matrix[0] = (-sin_lat * cos_lon);
    matrix[1] = (-sin_lon);
    matrix[2] = (-cos_lat * cos_lon);
    matrix[3] = (-sin_lat * sin_lon);
    matrix[4] = cos_lon;
    matrix[5] = (-cos_lat * sin_lon);
    matrix[6] = cos_lat;
    matrix[7] = 0;
    matrix[8] = (-sin_lat);
}
