#include<catch2/catch_test_macros.hpp>
#include<catch2/matchers/catch_matchers_floating_point.hpp>
#include<iostream>
#include"Ising system/spins.hpp"
#include"Ising system/Ising_system.cpp"

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

TEST_CASE("IsingSystem","[10 spins in 1D]")
{
    int n_spins=10;
    IsingSystem spin(n_spins);
    SECTION("initial"){
        REQUIRE(spin.J1()==-1.0);
        REQUIRE(spin._n_spins()==n_spins);
        REQUIRE(spin._maxrep_state()==1023);
    };
    
    SECTION("M 0f #7"){
        spin.set_state_by_code(7);
        REQUIRE(spin.eval_mz()==-4);
    };

    SECTION("E 0f #7"){
        spin.set_state_by_code(7);
        REQUIRE(spin.eval_energy_1D()==-6);
    };

    SECTION("M 0f #77"){
        spin.set_state_by_code(77);
        REQUIRE(spin.eval_mz()==-2);
    };

    SECTION("E 0f #77"){
        spin.set_state_by_code(77);
        REQUIRE(spin.eval_energy_1D()==-2);
    };

    SECTION("M 0f #777"){
        spin.set_state_by_code(777);
        REQUIRE(spin.eval_mz()==-2);
    };

    SECTION("E 0f #777"){
        spin.set_state_by_code(777);
        REQUIRE(spin.eval_energy_1D()==-2);
    };
};