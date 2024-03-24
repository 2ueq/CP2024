#include"IsingSystem_Square.hpp"
#include<iostream>
#include<vector>
#include

int main(){
    IsingSystem_Square spin;
    int L[4]={2,3,4,5};

    //初始化T（beta-1）
    std::vector<double>betaJ;
    for(double i=0.05;i<=4.0;i+=0.05){
        betaJ.push_back(1/i);
    };
    spin.beta=betaJ;

    for(){
        for(std::size_t beta_idx=0;beta_idx<beta.size();beta_idx++){

        };
    };

    return 0;
};