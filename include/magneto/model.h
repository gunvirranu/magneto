#ifndef MAGNETO_MODEL_H
#define MAGNETO_MODEL_H

#include <stddef.h>

#include "magneto.h"

#define MAGNETO_CALC_INDEX(n, m) (((n) * ((n) + 1) / 2) + (m) - 1)

typedef struct {
    const magneto_real g;
    const magneto_real h;
} magneto_SphericalHarmonicCoeff;

typedef struct {
    // Length is `magneto_Model.num_model_coeffs`
    const magneto_SphericalHarmonicCoeff *const coeffs;
} magneto_ModelCoeffs;

typedef struct {
    const magneto_DecYear epoch;
    const size_t nm_max;
    const size_t num_model_coeffs;
    const size_t num_models;
    const magneto_DecYear model_interval;
    // Length is `num_models`
    const magneto_ModelCoeffs *const models;
    const magneto_ModelCoeffs last_secular;
} magneto_Model;

magneto_FieldState eval_field(
    const magneto_Model *model,
    magneto_DecYear t,
    magneto_Coords coords
);

#endif  // MAGNETO_MODEL_H
