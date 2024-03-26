#include<catch2/catch_test_macros.hpp>
#include<catch2/matchers/catch_matchers_floating_point.hpp>
#include<iostream>
#include<vector>
// #include"Ising system/IsingSpins.hpp"
// #include"Ising system/Ising_system.hpp"
// #include"Ising system/IsingSpinOnLattice.hpp"
#include"Ising system/IsingSystem_Square.hpp"

TEST_CASE("IsingSpin","[single spin]"){
    IsingSpin spin;
    SECTION("spin state(initial)"){
        REQUIRE(spin._sz()== 1);
    };

    SECTION("set spin state as up(1)"){
        spin.set_up();
        REQUIRE(spin._sz()==1);
    };

    SECTION("set spin state as down(1)"){
        spin.set_down();
        REQUIRE(spin._sz()==-1);
    };

    SECTION("set spin state as up(2)"){
        spin.set_sz(1);
        REQUIRE(spin._sz()==1);
    };

    SECTION("set spin state as down(2)"){
        spin.set_sz(-1);
        REQUIRE(spin._sz()==-1);
    };
    
    SECTION("spin flip once"){
        spin.flip();
        REQUIRE(spin._sz()==-1);
    };

    SECTION("spin flip twice"){
        spin.flip();
        spin.flip();
        REQUIRE(spin._sz()==1);
    };
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
        REQUIRE(spin.eval_energy_1D()==2);
    };

    SECTION("M 0f #777"){
        spin.set_state_by_code(777);
        REQUIRE(spin.eval_mz()==-2);
    };

    SECTION("E 0f #777"){
        spin.set_state_by_code(777);
        REQUIRE(spin.eval_energy_1D()==-2);
    };

    SECTION("spin position initialization"){
        constexpr int dim=2;
        spin.set_dim(dim);
        for(int i=0;i<n_spins;i++){
            REQUIRE(spin._spin_position(i)==std::vector<int>({0,0}));
            REQUIRE(spin._spin_NN(i)==std::vector<int>({-1}));
            REQUIRE(spin._spin_NN(i,0)==-1);
        };
    };
};

TEST_CASE("IsingSpinOnLattice","single spin"){
    IsingSpinOnLattice spin;
    SECTION("spin position"){
        constexpr int dim=2;
        spin.set_dim(dim);
        REQUIRE(spin._position()==std::vector<int>({0,0}));
        REQUIRE(spin._NN()==std::vector<int>({-1}));    
        REQUIRE(spin._NN(0)==-1);
    }
};

TEST_CASE("IsingSystem_Square","[examples of 6*6 spins]"){
    const std::vector<int>system_size={6,6};
    IsingSystem_Square model(system_size);
    SECTION("basic"){
        REQUIRE(model._n_spins()==36);
        REQUIRE_THAT(model.J1(),Catch::Matchers::WithinULP(-1.0,4));
    };

    SECTION("site_index"){
        const std::vector<int>lattice_coordinate1={3,4};
        REQUIRE(model.site_index(lattice_coordinate1)==27);
        REQUIRE(model._lattice_coordinate(27)==lattice_coordinate1);
    };

    SECTION("neighboring site coordinates"){
        const std::vector<int>lattice_coordinate={3,3};
        const std::vector<int>lattice_coordinate_pos_x={4,3};
        const std::vector<int>lattice_coordinate_pos_y={3,4};
        const std::vector<int> lattice_coordinate_neg_x={2,3};
        const std::vector<int> lattice_coordinate_neg_y={3,2};
        REQUIRE(model.shift_pos_x(lattice_coordinate)==lattice_coordinate_pos_x);
        REQUIRE(model.shift_pos_y(lattice_coordinate)==lattice_coordinate_pos_y);
        REQUIRE(model.shift_neg_x(lattice_coordinate)==lattice_coordinate_neg_x);
        REQUIRE(model.shift_neg_y(lattice_coordinate)==lattice_coordinate_neg_y);
    };

    SECTION("connectivity"){
        constexpr int i=21;
        REQUIRE(model.NN(i,0)==22);
        REQUIRE(model.NN(i,1)==27);
        REQUIRE(model.NN(i,2)==20);
        REQUIRE(model.NN(i,3)==15);
    };

    SECTION("'pi'state:Magnectization and Energy"){
        std::vector<bool>state({1,1,0,1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,0,0,0,0,0,0,1,1,0,1,0,1,1,1,0,0,0,0});
        model.set_state(state);
        REQUIRE(model.eval_Mz()==2);
        REQUIRE_THAT(model.eval_energy(),Catch::Matchers::WithinULP(-4.0,4));
    };

        SECTION("'pi'state:Magnectization and Energy"){
        std::vector<bool>state({0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0});
        model.set_state(state);
        REQUIRE_THAT(model.eval_energy(),Catch::Matchers::WithinULP(-72.0,4));
    };

};

TEST_CASE("IsingSystem_Square","[tests of exact counting]"){
    const std::vector<int>system_size={6,6};
    std::vector<double>beta={0.1,1.0,2.0};
    IsingSystem_Square model(system_size,beta);
    std::vector<bool>pi_state({1,1,0,1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,0,0,0,0,0,0,1,1,0,1,0,1,1,1,0,0,0,0});
    constexpr double energy=-4.0;
    constexpr int magz=2;
    constexpr double w0=1.491824697641270;//exp(-1*0.1*(-4.0))
    constexpr double w1=54.59815003314424;//exp(-1*1.0*(-4.0))
    constexpr double w2=2980.957987041728;//exp(-1*2*(-4.0))
    constexpr double w[3]={w0,w1,w2};

    SECTION("'pi'state:M,E and Boltzman weight"){
        model.set_state(pi_state);
        REQUIRE(model.eval_Mz()==magz);
        REQUIRE_THAT(model.eval_energy(),Catch::Matchers::WithinULP(energy,4));
        for(std::size_t beta_idx=0;beta_idx<beta.size();beta_idx++){
            REQUIRE_THAT(model.weight_unnormalized(beta[beta_idx]),Catch::Matchers::WithinULP(w[beta_idx],4));
        };
    };

    SECTION("'pi'state:single term in the whole sum"){
        model.exactly_evaluate(pi_state);
        for(std::size_t beta_idx=0;beta_idx<beta.size();beta_idx++){
            REQUIRE_THAT(model._exact_energy_Z(beta_idx),Catch::Matchers::WithinULP(w[beta_idx],4));
            REQUIRE_THAT(model._exact_energy_q(beta_idx),Catch::Matchers::WithinULP(energy*w[beta_idx],4));
            REQUIRE_THAT(model._exact_energy_q_sq(beta_idx),Catch::Matchers::WithinULP(energy*energy*w[beta_idx],4));        
            REQUIRE_THAT(model._exact_magz_Z(beta_idx),Catch::Matchers::WithinULP(w[beta_idx],4));
            REQUIRE_THAT(model._exact_magz_q_sq(beta_idx),Catch::Matchers::WithinULP(magz*magz*w[beta_idx],4));
        };  
    };

};

