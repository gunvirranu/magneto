#ifndef MAGNETO_MAGNETO_H
#define MAGNETO_MAGNETO_H

#include <stdint.h>

#ifdef MAGNETO_SINGLE_PRECISION
typedef float magneto_real;
#define MAGNETO_REAL(x) x ## f
#else
typedef double magneto_real;
#define MAGNETO_REAL(x) x
#endif

/// Pi
extern const magneto_real magneto_PI;
/// Equatorial radius or semi-major axis of the WGS84 ellipse in [m]
extern const magneto_real magneto_WGS84_A;
/// Flattening of the WGS84 ellipsoid
extern const magneto_real magneto_WGS84_F;
/// Reciprocal flattening (1/f) of the WGS ellipsoid
extern const magneto_real magneto_WGS84_F_INV;
/// Eccentricity squared (e^2) of the WGS84 ellipsoid
extern const magneto_real magneto_WGS84_E_SQ;

/// Convert radians to degrees
magneto_real magneto_rad_to_deg(magneto_real rad);
/// Convert degrees to radians
magneto_real magneto_deg_to_rad(magneto_real deg);

typedef struct {
    magneto_real year;
} magneto_DecYear;

/// Earth-centered geographic coordinates
typedef struct {
    magneto_real latitude;      ///< Geodetic latitude [deg]
    magneto_real longitude;     ///< Longitude [deg]
    magneto_real height;        ///< Height above WGS84 reference ellipsoid [m]
} magneto_Coords;

typedef struct {
    magneto_real B_ned[3];
    magneto_real F;
    magneto_real H;
    magneto_real D;
    magneto_real I;
} magneto_FieldState;

// Date-time conversions

magneto_DecYear magneto_DecYear_from_datetime(
    int_fast16_t year, int_fast8_t month, int_fast8_t day,
    int_fast8_t hour, int_fast8_t minute, int_fast8_t sec
);

magneto_DecYear magneto_DecYear_from_date(
    int_fast16_t year, int_fast8_t month, int_fast8_t day
);

// Magnetic field conversions

magneto_FieldState magneto_FieldState_from_ned(const magneto_real *B_ned);

#endif  // MAGNETO_MAGNETO_H
