#ifndef MAGNETO_COMMON_PRIVATE_H
#define MAGNETO_COMMON_PRIVATE_H

#define STATIC_ASSERT(COND, MSG) typedef char static_assertion_##MSG[(COND) ? 1 : -1]

#define NS(x)       magneto_ ## x
#define REAL        MAGNETO_REAL

#define SIN         REAL(sin)
#define COS         REAL(cos)
#define SQRT        REAL(sqrt)
#define HYPOT       REAL(hypot)

#define rad_to_deg  NS(rad_to_deg)
#define deg_to_rad  NS(deg_to_rad)

typedef NS(real)            real;
typedef NS(DecYear)         DecYear;
typedef NS(Coords)          Coords;
typedef NS(SphericalCoords) SphericalCoords;
typedef NS(EcefPosition)    EcefPosition;
typedef NS(FieldState)      FieldState;

static inline real sq(const real x) {
    return x * x;
}

#endif  // MAGNETO_COMMON_PRIVATE_H
