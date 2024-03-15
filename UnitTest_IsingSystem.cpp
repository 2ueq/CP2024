#include<catch2/catch_test_macros.hpp>
#include<catch2/matchers/catch_matchers_floating_point.hpp>
#include<iostream>
#include"Ising system/Ising_system.cpp"

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