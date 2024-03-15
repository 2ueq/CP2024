#include<catch2/catch_test_macros.hpp>
#include<catch2/matchers/catch_matchers_floating_point.hpp>
#include<iostream>
#include"Ising system/spins.hpp"

TEST_CASE("IsingSpin","[single spin]"){
    IsingSpin spin;
    SECTION("spin state(initial)"){
        REQUIRE(spin._sz()== 1);
    }

    SECTION("set spin state as up(1)"){
        spin.set_up();
        REQUIRE(spin._sz()==1);
    }

    SECTION("set spin state as down(1)"){
        spin.set_down();
        REQUIRE(spin._sz()==-1);
    }

    SECTION("set spin state as up(2)"){
        spin.set_sz(1);
        REQUIRE(spin._sz()==1);
    }

    SECTION("set spin state as down(2)"){
        spin.set_sz(-1);
        REQUIRE(spin._sz()==-1);
    }
    
    SECTION("spin flip once"){
        spin.flip();
        REQUIRE(spin._sz()==-1);
    }

    SECTION("spin flip twice"){
        spin.flip();
        spin.flip();
        REQUIRE(spin._sz()==1);
    }
};