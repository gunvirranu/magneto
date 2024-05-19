#include "magneto/model.h"

#include <math.h>
#include <stdbool.h>
#include <stddef.h>

#include <stdio.h>

#include "common_private.h"

static inline real calc_K(const real n, const real m) {
    if (n <= 1U) {
        return 0;
    }
    return (sq(n - 1) - sq(m)) / (((2 * n) - 1) * ((2 * n) - 3));
}

static void rotate_vector_spherical_to_ned(
    const Coords pos,
    const SphericalCoords pos_sph,
    const SphericalCoords B_spherical,
    real *const B_ned
) {
    const real B_r = B_spherical.radius;
    const real B_theta = B_spherical.azimuth;
    const real B_phi = B_spherical.polar;
    const real eps = deg_to_rad(pos_sph.polar - pos.latitude);
    const real sin_eps = SIN(eps);
    const real cos_eps = COS(eps);
    B_ned[0] = (-B_theta * cos_eps) - (B_r * sin_eps);
    B_ned[1] = B_phi;
    B_ned[2] = (B_theta * sin_eps) - (B_r * cos_eps);
}

static inline void calc_g_and_h(
    const magneto_Model *const model,
    const size_t i_model,
    const magneto_real t,
    const size_t n,
    const size_t m,
    real *const g_n_m,
    real *const h_n_m
) {
    const magneto_SphericalHarmonicCoeff *const coeffs_i = model->models[i_model].coeffs;
    const bool is_last_submodel = ((i_model + 1U) >= model->num_models);
    const size_t idx_coeff = MAGNETO_CALC_INDEX(n, m);

    printf("(%d, %d) -> %d\n", (int) n, (int) m, (int) idx_coeff);

    real g_dot = 0;
    real h_dot = 0;
    if (is_last_submodel) {
        g_dot = model->last_secular.coeffs[idx_coeff].g;
        h_dot = model->last_secular.coeffs[idx_coeff].h;
    } else {
        const magneto_SphericalHarmonicCoeff *const coeffs_next = model->models[i_model + 1U].coeffs;
        g_dot = (coeffs_next[idx_coeff].g - coeffs_i[idx_coeff].g) / model->model_interval.year;
        h_dot = (coeffs_next[idx_coeff].h - coeffs_i[idx_coeff].h) / model->model_interval.year;
    }

    *g_n_m = (coeffs_i->g + (t * g_dot));
    *h_n_m = (coeffs_i->h + (t * h_dot));
}

/// Compute vector as gradient of spherical harmonic potential expansion
///
/// @todo There is definitely some additional optimization to be had here.
///       Need to get some basic benchmarking set up and test out some stuff.
///
/// @param[in]  model           Spherical harmonic model and coefficients
/// @param[in]  i_model         Index of which sub-model to use
/// @param[in]  t               Time delta into i-th sub-model in years
/// @param[in]  pos             Geocentric spherical coordinates
/// @param[out] B_spherical     Output vector in spherical reference frame
static void eval_spherical_expansion(
    const magneto_Model *const model,
    const size_t i_model,
    const real t,
    const SphericalCoords pos,
    SphericalCoords *const B_spherical
) {
    const real theta = deg_to_rad(pos.azimuth);
    const real phi = deg_to_rad(pos.polar);
    const real sin_theta = SIN(theta);
    const real cos_theta = COS(theta);
    const real normed_r = (magneto_WGS84_A / pos.radius);

    // Output vector
    real B_r = 0;       ///< Bz
    real B_theta = 0;   ///< Bx
    real B_phi = 0;     ///< By

    // Recursive values for P_{n,m} and dP_{n,m}/dtheta
    real P_n_n = 1;
    real dP_n_n = 0;

    // Compute Gaussian normalized associated Legendre polynomials recursively
    for (size_t m = 0U; m <= model->nm_max; ++m) {
        const real P_nprev_nprev = P_n_n;
        const real dP_nprev_nprev = dP_n_n;

        // Compute P_{n,n}, but skip first iter
        P_n_n = (m > 0U) ? (sin_theta * P_nprev_nprev) : 1;
        real P_nprev_m = 1;
        real P_nprevprev_m = 0;

        // Compute dP_{n,n}, but skip first iter
        dP_n_n = (m > 0U) ? ((sin_theta * dP_nprev_nprev) + (cos_theta * P_nprev_nprev)) : 0;
        real dP_nprev_m = 0;
        real dP_nprevprev_m = 0;

        const real sin_mphi = SIN((real) m * phi);
        const real cos_mphi = COS((real) m * phi);

        // Condition is enforced by loop bounds: (m <= n)
        for (size_t n = MAX_OF(m, 1U); n <= model->nm_max; ++n) {

            // Default values for first iter, computed pre-loop
            real P_n_m = P_n_n;     // From  P_{n,n}
            real dP_n_m = dP_n_n;   // From dP_{n,n}

            // Skipping first iteration, otherwise compute P_{n,m} and dP_{n,m}
            if (n != m) {
                const real K_n_m = calc_K((real) n, (real) m);

                P_n_m = (cos_theta * P_nprev_m) - (K_n_m * P_nprevprev_m);
                dP_n_m = (cos_theta * dP_nprev_m) - (sin_theta * P_nprev_m) - (K_n_m * dP_nprevprev_m);
            }

            // Save recursive values for next iter
            P_nprevprev_m = P_nprev_m;
            P_nprev_m = P_n_m;
            dP_nprevprev_m = dP_nprev_m;
            dP_nprev_m = dP_n_m;

            printf("P_{n=%zu}{m=%zu} = %.10f, dP = %f\n", n, m, P_n_m, dP_n_m);

            real g_n_m = 0;
            real h_n_m = 0;
            calc_g_and_h(model, i_model, t, n, m, &g_n_m, &h_n_m);

            const real r_scalar = REAL(pow)(normed_r, (real) (n + 2U));

            B_r += (r_scalar * ((g_n_m * cos_mphi) + (h_n_m * sin_mphi)) * ((real) (n + 1U)) * P_n_m);
            B_theta -= (r_scalar * ((g_n_m * cos_mphi) + (h_n_m * sin_mphi)) * dP_n_m);
            B_phi -= (r_scalar * ((real) m) * ((-g_n_m * sin_mphi) + (h_n_m + cos_mphi)) * P_n_m);
        }
    }

    if (sin_theta != 0) {
        B_phi /= sin_theta;
    }
    B_spherical->radius = B_r;
    B_spherical->azimuth = B_theta;
    B_spherical->polar = B_phi;
}

magneto_FieldState eval_field(
    const magneto_Model *const model,
    const magneto_DecYear t,
    const magneto_Coords coords
) {
    const SphericalCoords sph = magneto_SphericalCoords_from_coords(coords);
    const real delta_t = (t.year - model->epoch.year);

    // Evaluate magnetic field model in spherical coordinates
    SphericalCoords B_spherical = { 0 };
    eval_spherical_expansion(model, 0, delta_t, sph, &B_spherical);

    // Rotate magnetic field vector from geocentric to geodetic NED frame
    real B_ned[3];
    rotate_vector_spherical_to_ned(coords, sph, B_spherical, B_ned);

    // Compute other field quantities
    const FieldState B = magneto_FieldState_from_ned(B_ned);
    return B;
}
