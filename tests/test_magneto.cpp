#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

extern "C" {
#  include <magneto/magneto.h>
#  include <magneto/model.h>
#  include <magneto/wmm.h>

// Total hack for unit-testing static stuff
#  include <magneto/../../src/magneto.c>
}

using doctest::Approx;

using real = magneto_real;

static_assert(sizeof(real) == sizeof(double), "Tests expecting double-precision");
static_assert(sizeof(real) != sizeof(float), "Tests not expecting single-precision just yet");

TEST_CASE(
    "test_constants"
    * doctest::description("Sanity check on library constants")
) {
    CHECK(magneto_PI == Approx(3.14).epsilon(0.01));
    CHECK(magneto_rad_to_deg(magneto_PI / 4) == Approx(45));
    CHECK(magneto_deg_to_rad(60) == Approx(magneto_PI / 3));
    CHECK(RAD_PER_DEG * DEG_PER_RAD == Approx(1.0).epsilon(1e-14));

    CHECK(magneto_WGS84_A == Approx(magneto_WGS84_B).epsilon(0.1));
    CHECK(magneto_WGS84_F == Approx((magneto_WGS84_A - magneto_WGS84_B) / magneto_WGS84_A).epsilon(1e-10));
    CHECK(magneto_WGS84_F == Approx(1.0 / magneto_WGS84_F_INV).epsilon(1e-14));
    CHECK(magneto_WGS84_E_SQ == Approx(6.69437999014e-3).epsilon(1e-14));

    CHECK((MONTH_MAX + 1U) == ARRAY_SIZE(DAY_MAX));
    CHECK((MONTH_MAX + 1U) == ARRAY_SIZE(DAY_PER_MONTH));
    CHECK(DAYS_IN_YEAR == (DAY_PER_MONTH[MONTH_MAX] + DAY_MAX[MONTH_MAX]));
}

TEST_CASE("test_is_leap_year") {
    CHECK(is_leap_year(2000U));
    CHECK(is_leap_year(2024U));
    CHECK_FALSE(is_leap_year(2019U));
    CHECK_FALSE(is_leap_year(1900U));
    CHECK_FALSE(is_leap_year(2100U));
}

TEST_CASE("test_days_in_year") {
    CHECK(days_in_year(2001U) == 365U);
    CHECK(days_in_year(2000U) == 366U);
}

TEST_CASE("test_days_in_month") {
    CHECK(days_in_month(1999U, 2U) == 28U);
    CHECK(days_in_month(1998U, 8U) == 31U);
    CHECK(days_in_month(2008U, 2U) == 29U);
    CHECK(days_in_month(2012U, 4U) == 30U);
    CHECK(days_in_month(2018U, 0U) == 0U);
    CHECK(days_in_month(2020U, 13U) == 0U);

    uint16_t days = 0U;
    for (size_t i = 1; i < 13U; ++i) {
        CHECK(days == DAY_PER_MONTH[i]);
        days += DAY_MAX[i];
    }
}

TEST_CASE("test_day_of_year") {
    magneto_DateTime t = (magneto_DateTime) {
        .year = 2001U,
        .month = 1U,
        .day = 1U,
        .hour = 0U,
        .minute = 0U,
        .sec = 0U,
    };
    CHECK(day_of_year(t) == 1.0);

    t.month = 9U;
    t.day = 14U;
    CHECK(day_of_year(t) == 257.0);
    t.day += 1;
    CHECK(day_of_year(t) == 258.0);
    t.month += 1U;
    CHECK(day_of_year(t) == 288.0);
    t.year = 2024U;
    CHECK(day_of_year(t) == 289.0);

    t.hour = 12U;
    CHECK(day_of_year(t) == 289.5);
    t.minute = 16U;
    CHECK(day_of_year(t) == Approx(289.5111111111).epsilon(1e-10));
    t.minute = 0U;
    t.sec = 32U;
    CHECK(day_of_year(t) == Approx(289.5003703704).epsilon(1e-10));

    t.month = 0U;
    CHECK_NOTHROW(day_of_year(t));
    t.month = 13U;
    CHECK_NOTHROW(day_of_year(t));
}

TEST_CASE("test_date_time_is_valid") {
    magneto_DateTime t = (magneto_DateTime) {
        .year = 2024U,
        .month = 5U,
        .day = 20U,
        .hour = 23U,
        .minute = 48U,
        .sec = 34U,
    };
    CHECK(magneto_DateTime_is_valid(t));

    t.year = 1582U;
    CHECK_FALSE(magneto_DateTime_is_valid(t));
    t.year = 10000U;
    CHECK_FALSE(magneto_DateTime_is_valid(t));
    t.year = 2071U;
    CHECK(magneto_DateTime_is_valid(t));

    t.month = 0U;
    CHECK_FALSE(magneto_DateTime_is_valid(t));
    t.month = 13U;
    CHECK_FALSE(magneto_DateTime_is_valid(t));
    t.month = 11U;
    CHECK(magneto_DateTime_is_valid(t));

    t.day = 0U;
    CHECK_FALSE(magneto_DateTime_is_valid(t));
    t.day = 31U;
    CHECK_FALSE(magneto_DateTime_is_valid(t));
    t.month = 8U;
    CHECK(magneto_DateTime_is_valid(t));

    t.year = 2032U;
    t.month = 2U;
    t.day = 29U;
    CHECK(magneto_DateTime_is_valid(t));
    t.year = 2079U;
    CHECK_FALSE(magneto_DateTime_is_valid(t));
    t.day = 28U;
    CHECK(magneto_DateTime_is_valid(t));

    t.hour = 25U;
    CHECK_FALSE(magneto_DateTime_is_valid(t));
    t.hour = 13U;
    CHECK(magneto_DateTime_is_valid(t));

    t.minute = 60U;
    CHECK_FALSE(magneto_DateTime_is_valid(t));
    t.minute = 59U;
    CHECK(magneto_DateTime_is_valid(t));

    t.sec = 71U;
    CHECK_FALSE(magneto_DateTime_is_valid(t));
    t.sec = 60U;
    CHECK(magneto_DateTime_is_valid(t));
}

TEST_CASE("test_dec_year_is_valid") {
    magneto_DecYear t;
    t.year = REAL(1987.4);
    CHECK(magneto_DecYear_is_valid(t));
    t.year = REAL(1582.9);
    CHECK_FALSE(magneto_DecYear_is_valid(t));
    t.year = REAL(10000.1);
    CHECK_FALSE(magneto_DecYear_is_valid(t));
    t.year = REAL(2101.0);
    CHECK(magneto_DecYear_is_valid(t));
}

TEST_CASE("test_zero_dec_year_date_time") {
    const magneto_DateTime t0_dt = {};
    const magneto_DecYear  t0_dy = {};

    CHECK(t0_dy.year == 0.0);
    CHECK_FALSE(magneto_DateTime_is_valid(t0_dt));
    CHECK_FALSE(magneto_DecYear_is_valid(t0_dy));

    const magneto_DecYear t1_dy = magneto_DecYear_from_date_time(t0_dt);
    CHECK(t1_dy.year == t0_dy.year);
}

TEST_CASE("test_dec_year_from_date_time") {
    const magneto_DateTime t_dt = (magneto_DateTime) {
        .year = 2020U,
        .month = 3U,
        .day = 1U,
        .hour = 7U,
        .minute = 6U,
        .sec = 5U,
    };
    const magneto_DecYear t_dy = magneto_DecYear_from_date_time(t_dt);
    CHECK((uint16_t) t_dy.year == 2020U);
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
