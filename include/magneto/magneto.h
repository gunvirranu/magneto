#ifndef MAGNETO_MAGNETO_H
#define MAGNETO_MAGNETO_H

#include <stdint.h>
#include <stdbool.h>

// Default is double-precision if not defined
#ifdef MAGNETO_SINGLE_PRECISION
typedef float magneto_real;
#else
typedef double magneto_real;
#endif

// Constants

/// [ ] Pi
extern const magneto_real magneto_PI;
/// [m] Equatorial radius or semi-major axis of the WGS84 ellipse
extern const magneto_real magneto_WGS84_A;
/// [m] Semi-minor axis of the WGS84 ellipsoid
extern const magneto_real magneto_WGS84_B;
/// [ ] Flattening of the WGS84 ellipsoid
extern const magneto_real magneto_WGS84_F;
/// [ ] Reciprocal flattening (1/f) of the WGS ellipsoid
extern const magneto_real magneto_WGS84_F_INV;
/// [ ] Eccentricity squared (e^2) of the WGS84 ellipsoid
extern const magneto_real magneto_WGS84_E_SQ;

// Time

typedef struct {
    magneto_real year;  ///< [year] Decimal year with fractional day & time in [1583, 9999]
} magneto_DecYear;

/// Follows ISO8601, based on Gregorian calendar
typedef struct {
    uint16_t year;  ///< [year]     Integer Gregorian year in [1583, 9999]
    uint8_t month;  ///< [month]    Month of year in [1, 12], starting at January
    uint8_t day;    ///< [day]      Day of month in [1, 31]
    uint8_t hour;   ///< [hour]     Hour of day in [0, 24)
    uint8_t minute; ///< [min]      Minute of hour in [0, 60)
    uint8_t sec;    ///< [sec]      Second in [0, 60)
} magneto_DateTime;

// Position

/// Earth-centered geographic coordinates (geodetic)
typedef struct {
    magneto_real latitude;      ///< [deg]  Geodetic latitude
    magneto_real longitude;     ///< [deg]  Longitude
    magneto_real height;        ///< [m]    Height above WGS84 reference ellipsoid
} magneto_Coords;

/// Earth-centered spherical coordinates (geocentric)
typedef struct {
    magneto_real polar;         ///< [deg]  Polar angle (a.k.a. geoentric latitude)
    magneto_real azimuth;       ///< [deg]  Azimuthal angle (a.k.a. longitude)
    magneto_real radius;        ///< [m]    Distance from centre of WGS84 ellipsoid
} magneto_SphericalCoords;

/// Earth-centered cartesian coordinates (TODO: Define axes)
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

bool magneto_DateTime_is_valid(magneto_DateTime t);
bool magneto_DecYear_is_valid(magneto_DecYear t);

magneto_DecYear magneto_DecYear_from_date_time(magneto_DateTime t);
// magneto_DateTime magneto_DateTime_from_dec_year(magneto_DecYear t);  TODO: Do we need this?

// Position conversions

magneto_Coords magneto_Coords_from_spherical(magneto_SphericalCoords pos);
magneto_Coords magneto_Coords_from_ecef(magneto_EcefPosition pos);

magneto_SphericalCoords magneto_SphericalCoords_from_coords(magneto_Coords pos);
magneto_SphericalCoords magneto_SphericalCoords_from_ecef(magneto_EcefPosition pos);

magneto_EcefPosition magneto_EcefPosition_from_coords(magneto_Coords pos);
magneto_EcefPosition magneto_EcefPosition_from_spherical(magneto_SphericalCoords pos);

// Magnetic field conversions

magneto_FieldState magneto_FieldState_from_ned(const magneto_real *B_ned);

void magneto_convert_vector_ned_to_ecef(
    magneto_Coords pos,
    const magneto_real *ned,
    magneto_real *ecef
);
void magneto_convert_vector_ecef_to_ned(
    magneto_Coords pos,
    const magneto_real *ecef,
    magneto_real *ned
);
void magneto_matrix_ned_to_ecef(magneto_Coords pos, magneto_real *matrix);

#endif  // MAGNETO_MAGNETO_H
