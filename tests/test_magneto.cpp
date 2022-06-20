#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

extern "C" {
#  include <magneto/magneto.h>
}

using doctest::Approx;

TEST_CASE("test_constants") {
    CHECK(magneto_PI == Approx(3.14).epsilon(0.01));
    CHECK(magneto_rad_to_deg(magneto_PI / 4) == Approx(45));
    CHECK(magneto_deg_to_rad(60) == Approx(magneto_PI / 3));
}

