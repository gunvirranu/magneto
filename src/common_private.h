#ifndef MAGNETO_COMMON_PRIVATE_H
#define MAGNETO_COMMON_PRIVATE_H

// General macros

#define STATIC_ASSERT(COND, MSG) typedef char static_assertion_##MSG[(COND) ? 1 : -1]

// Support single and double precision floating point

#ifdef MAGNETO_SINGLE_PRECISION
#define REAL(x) x ## f
#else
#define REAL(x) x
#endif

#define SIN         REAL(sin)
#define COS         REAL(cos)
#define SQRT        REAL(sqrt)
#define HYPOT       REAL(hypot)

// Avoid some typing in source files

#define NS(x)       magneto_ ## x

#define rad_to_deg  NS(rad_to_deg)
#define deg_to_rad  NS(deg_to_rad)

typedef NS(real)            real;
typedef NS(DecYear)         DecYear;
typedef NS(Coords)          Coords;
typedef NS(SphericalCoords) SphericalCoords;
typedef NS(EcefPosition)    EcefPosition;
typedef NS(FieldState)      FieldState;

// General utility functions

static inline real sq(const real x) {
    return x * x;
}

#endif  // MAGNETO_COMMON_PRIVATE_H
