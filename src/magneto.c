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

static const uint16_t YEAR_MIN = 1583U;     ///< Beginning of valid Gregorian calendar
static const uint16_t YEAR_MAX = 9999U;
static const uint8_t MONTH_MIN = 1U;
static const uint8_t MONTH_MAX = 12U;
static const uint8_t DAY_MIN = 1U;
/// Number of days in each month, indexed by month (0 is unused)
static const uint8_t DAY_MAX[13U] = { 0U, 31U, 28U, 31U, 30U, 31U, 30U, 31U, 31U, 30U, 31U, 30U, 31U };
static const uint8_t HOUR_MIN = 0U;
static const uint8_t HOUR_MAX = 24U;
static const uint8_t MINUTE_MIN = 0U;
static const uint8_t MINUTE_MAX = 59U;
static const uint8_t SEC_MIN = 0U;
static const uint8_t SEC_MAX = 60U;

static const uint16_t DAYS_IN_YEAR = 365U;  ///< Not in a leap year
/// Cumulative number of days in year at start of month, index by month (0 is unused)
static const uint16_t DAY_PER_MONTH[13U] = { 0U, 0U, 31U, 59U, 90U, 120U, 151U, 181U, 212U, 243U, 273U, 304U, 334U };

real magneto_rad_to_deg(const real rad) {
    return rad * DEG_PER_RAD;
}

real magneto_deg_to_rad(const real deg) {
    return deg * RAD_PER_DEG;
}

/// Based on Gregorian calendar, assumes valid `year`
static bool is_leap_year(const uint16_t year) {
    const bool div_4   = (year % 4U)   == 0U;
    const bool div_100 = (year % 100U) == 0U;
    const bool div_400 = (year % 400U) == 0U;
    return (div_4 && !div_100) || div_400;
}

/// Assumes valid `year`
static uint16_t days_in_year(const uint16_t year) {
    // Add extra day for February 29th in leap years
    return (DAYS_IN_YEAR + (is_leap_year(year) ? 1U : 0U));
}

/// Assumes valid `year` and `month`
static uint8_t days_in_month(const uint16_t year, const uint8_t month)
{
    // Check month b/c it's used for indexing
    uint8_t day = DAY_MAX[(month < ARRAY_SIZE(DAY_MAX)) ? month : 0U];
    // Add extra day for February 29th in leap years
    if ((month == 2U) && is_leap_year(year))
    {
        day += 1U;
    }
    return day;
}

/// Mostly assumes a valid `Datetime`
static real day_of_year(const DateTime t) {
    // Check for invalid month only b/c it's used for indexing
    real doy = DAY_PER_MONTH[(t.month < ARRAY_SIZE(DAY_PER_MONTH)) ? t.month : 0U];
    // Add an extra day after February 29th in a leap year
    if ((t.month > 2U) && is_leap_year(t.year)) {
        doy += REAL(1.0);
    }
    // Add day and time simply with no checks
    doy += (real) t.day;
    doy += (t.hour / REAL(24.0));
    doy += (t.minute / REAL(1440.0));
    doy += (t.sec / REAL(86400.0));
    return doy;
}

bool magneto_DateTime_is_valid(const DateTime t) {
    bool valid = true;
    valid &= (YEAR_MIN <= t.year) && (t.year <= YEAR_MAX);
    valid &= (MONTH_MIN <= t.month) && (t.month <= MONTH_MAX);
    valid &= (DAY_MIN <= t.day) && (t.day <= days_in_month(t.year, t.month));
    valid &= (HOUR_MIN <= t.hour) && (t.hour <= HOUR_MAX);
    valid &= (MINUTE_MIN <= t.minute) && (t.minute <= MINUTE_MAX);
    valid &= (SEC_MIN <= t.sec) && (t.sec <= SEC_MAX);
    return valid;
}

bool magneto_DecYear_is_valid(const DecYear t) {
    return (t.year >= (real) YEAR_MIN) && (t.year <= (real) YEAR_MAX);
}

DecYear magneto_DecYear_from_date_time(const DateTime t) {
    if (!magneto_DateTime_is_valid(t)) {
        return (DecYear) { .year = REAL(0.0) };
    }
    const real doy = day_of_year(t);
    const real frac_year = (doy - REAL(1.0)) / days_in_year(t.year);
    return (DecYear) {
        .year = (real) t.year + frac_year
    };
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

SphericalCoords  magneto_SphericalCoords_from_ecef(const EcefPosition pos) {
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
    // Standard matrix multiplication, transposed
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
