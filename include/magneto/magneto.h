#ifndef MAGNETO_MAGNETO_H
#define MAGNETO_MAGNETO_H

#ifdef MAGNETO_SINGLE_PRECISION
typedef float magneto_real;
#define MAGNETO_REAL(x) x ## f
#else
typedef double magneto_real;
#define MAGNETO_REAL(x) x
#endif

/// Pi
extern const magneto_real magneto_pi;
/// Equatorial radius or semi-major axis of the WGS84 ellipse in [m]
extern const magneto_real magneto_wgs84_a;
/// Flattening of the WGS84 ellipsoid
extern const magneto_real magneto_wgs84_f;
/// Reciprocal flattening (1/f) of the WGS ellipsoid
extern const magneto_real magneto_wgs84_f_inv;
/// Eccentricity squared (e^2) of the WGS84 ellipsoid
extern const magneto_real magneto_wgs84_e_sq;

/// Convert radians to degrees
magneto_real magneto_rad_to_deg(magneto_real rad);
/// Convert degrees to radians
magneto_real magneto_deg_to_rad(magneto_real deg);

#endif  // MAGNETO_MAGNETO_H
