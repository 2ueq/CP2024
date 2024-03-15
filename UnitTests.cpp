#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("TestCatch2", "[SimplestExample]") {
	REQUIRE(2 * 2 == 4);
};

TEST_CASE("TestCatch2", "[Using SECTION's]") {
	int i = 0;
	SECTION("Comparison 1") {
		i += 1;
		REQUIRE(i == 1);
	}
	SECTION("Comparison 2") {
		REQUIRE(i == 0);
	}
};

TEST_CASE("TestCatch2", "[Floating Point Comparisons]") {
	const double ref = 1.0;
	
	SECTION("WithinAbs") {
		// WithinAbs(double target, double margin) matches if difference relative to target is less-or-equal to the margin
		REQUIRE_THAT(1.2, Catch::Matchers::WithinAbs(ref, 0.2));
	}
	SECTION("WithinRel") {
		// WithinRel(FloatingPoint target, FloatingPoint eps) matches if |arg - target| <= eps * max(|arg|, |target|)
		REQUIRE_THAT(1.2, Catch::Matchers::WithinRel(ref, 0.2));
	}
	SECTION("WithinULP (OK)") {
		// WithinULP(FloatingPoint target, uint64_t maxUlpDiff) matches if floating point numbers that are no more than maxUlpDiff ULPs away from the target value.
		REQUIRE_THAT( std::sin(0.5 * M_PI), Catch::Matchers::WithinULP( ref, 4 ) );
	}
	SECTION("WithinULP (OK/NG)") {
		constexpr double x_OK = 1.0000000000000001;
		constexpr double x_NG = 1.000000000000001;
		constexpr bool use_OK = true;
		double x = use_OK ? x_OK : x_NG;
		REQUIRE_THAT(x, Catch::Matchers::WithinULP( ref, 4 ) );
	}
};

