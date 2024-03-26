#include"IsingSystem_Square.hpp"
#include<iostream>
#include<vector>
using namespace std;

int main(){
    double start=0.05;
    double end=4.0;
    double step=0.05;

    //初始化T（beta-1）
    std::vector<double>Temperature;
    std::vector<double>betaJ;
    for(double i=start;i<end;i+=step){
        Temperature.push_back(i);
        betaJ.push_back(1.0/i);
    };
    
    const std::vector<int>system_size={5,5};
    IsingSystem_Square spin(system_size,betaJ);
    
    spin.exact();
    spin.print_exact();

    return 0;
};

// betaJ={0.1,1,4};
// for(int i = 0; i < 16; i++){
    //     cout << "for configuration idx:" << i << ": energy*weight=" << spin._exact_energy_q(80,i) << endl;
    // }