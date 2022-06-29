#ifndef MAGNETO_COMMON_PRIVATE_H
#define MAGNETO_COMMON_PRIVATE_H

#define NS(x) magneto_ ## x
#define REAL(x) MAGNETO_REAL(x)

typedef NS(real) real;
typedef NS(FieldState) FieldState;
typedef NS(Coords) Coords;
typedef NS(DecYear) DecYear;

#endif  // MAGNETO_COMMON_PRIVATE_H
