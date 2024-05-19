#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

extern "C" {
#  include <magneto/magneto.h>
#  include <magneto/model.h>
#  include <magneto/wmm.h>
}

using doctest::Approx;

using real = double;

TEST_CASE("test_constants") {
    CHECK(magneto_PI == Approx(3.14).epsilon(0.01));
    CHECK(magneto_rad_to_deg(magneto_PI / 4) == Approx(45));
    CHECK(magneto_deg_to_rad(60) == Approx(magneto_PI / 3));
}

TEST_CASE("test_dec_year") {
    CHECK(((magneto_DecYear) { 2014.513 }).year == 2014.513);
    const magneto_DecYear t = magneto_DecYear_from_datetime(2020, 3, 1, 7, 6, 5);
    CHECK((int) t.year == 2020);
}

TEST_CASE("test_wmm2020") {
    const magneto_DecYear t = { .year = 2020 };
    const magneto_Coords pos = { .latitude = 80, .longitude = 0, .height = 0 };

    magneto_SphericalCoords sph = magneto_SphericalCoords_from_coords(pos);
    sph.azimuth = 30.0;

    MESSAGE("r = ", sph.radius, ", azi = ", sph.azimuth, ", pol = ", sph.polar);

    const magneto_FieldState B = eval_field(&MODEL_WMM2020, t, pos);

    MESSAGE("B_ned = { ", B.B_ned[0], ", ", B.B_ned[1], ", ", B.B_ned[2], " }");
}
