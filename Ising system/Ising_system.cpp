#include"Ising_system.hpp"
int main(){
    IsingSystem a(10);
    a.set_state_by_code(7);
    std::cout<<"#7的M="<<a.eval_mz()<<std::endl;
    std::cout<<"#7的E="<<a.eval_energy_1D()<<std::endl;

    IsingSystem b(10);
    b.set_state_by_code(77);
    std::cout<<"#77的M="<<b.eval_mz()<<std::endl;
    std::cout<<"#77的E="<<b.eval_energy_1D()<<std::endl;

    IsingSystem c(10);
    c.set_state_by_code(777);
    std::cout<<"#777的M="<<c.eval_mz()<<std::endl;
    std::cout<<"#777的E="<<c.eval_energy_1D()<<std::endl;
    return 0;
};