#ifndef MAGNETO_MAGNETO_H
#define MAGNETO_MAGNETO_H

#include <stdint.h>

#ifdef MAGNETO_SINGLE_PRECISION
typedef float magneto_real;
#else
typedef double magneto_real;
#endif

// Constants

/// Pi
extern const magneto_real magneto_PI;
/// Equatorial radius or semi-major axis of the WGS84 ellipse in [m]
extern const magneto_real magneto_WGS84_A;
/// Semi-minor axis of the WGS84 ellipsoid in [m]
extern const magneto_real magneto_WGS84_B;
/// Flattening of the WGS84 ellipsoid
extern const magneto_real magneto_WGS84_F;
/// Reciprocal flattening (1/f) of the WGS ellipsoid
extern const magneto_real magneto_WGS84_F_INV;
/// Eccentricity squared (e^2) of the WGS84 ellipsoid
extern const magneto_real magneto_WGS84_E_SQ;

// Time

typedef struct {
    magneto_real year;
} magneto_DecYear;

// Position

/// Earth-centered geographic coordinates
typedef struct {
    magneto_real latitude;      ///< Geodetic latitude [deg]
    magneto_real longitude;     ///< Longitude [deg]
    magneto_real height;        ///< Height above WGS84 reference ellipsoid [m]
} magneto_Coords;

typedef struct {
    magneto_real radius;        ///< Distance from centre of WGS84 ellipsoid [m]
    magneto_real polar;         ///< Polar angle (a.k.a. geoentric latitude) [deg]
    magneto_real azimuth;       ///< Azimuthal angle (a.k.a. longitude) [deg]
} magneto_SphericalCoords;

typedef struct {
    magneto_real x;
    magneto_real y;
    magneto_real z;
} magneto_EcefPosition;

// Magnetic field

typedef struct {
    magneto_real B_ned[3];
    magneto_real F;
    magneto_real H;
    magneto_real D;
    magneto_real I;
} magneto_FieldState;

// Math

/// Convert radians to degrees
magneto_real magneto_rad_to_deg(magneto_real rad);
/// Convert degrees to radians
magneto_real magneto_deg_to_rad(magneto_real deg);

// Date-time conversions

magneto_DecYear magneto_DecYear_from_datetime(
    int_fast16_t year, int_fast8_t month, int_fast8_t day,
    int_fast8_t hour, int_fast8_t minute, int_fast8_t sec
);
magneto_DecYear magneto_DecYear_from_date(
    int_fast16_t year, int_fast8_t month, int_fast8_t day
);

// Position conversions

magneto_Coords magneto_Coords_from_ecef(magneto_EcefPosition pos);
magneto_SphericalCoords magneto_SphericalCoords_from_coords(magneto_Coords pos);
magneto_SphericalCoords magneto_SphericalCoords_from_ecef(magneto_EcefPosition pos);
magneto_EcefPosition magneto_EcefPosition_from_coords(magneto_Coords pos);
magneto_EcefPosition magneto_EcefPosition_from_spherical(magneto_SphericalCoords pos);

// Magnetic field conversions

magneto_FieldState magneto_FieldState_from_ned(const magneto_real *B_ned);

#endif  // MAGNETO_MAGNETO_H
