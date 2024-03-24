#include<iostream>
#include<vector>
#include"IsingSpins.hpp"

class IsingSpinOnLattice:public IsingSpin{
    private:
        std::vector<int>position;
        std::vector<int>NN;
    public:
        IsingSpinOnLattice(){
            set_dim(1);
            NN={-1};
        };
        ~IsingSpinOnLattice(){};
        void set_dim(int dim){
            position.assign(dim,0);
        };

        std::vector<int>_position()const{return position;};
        std::vector<int>_NN()const{return NN;};
        int _NN(const int bond_idx)const{return NN[bond_idx];};

        void set_NN(const int bond_idx,const int site_idx){NN[bond_idx]=site_idx;};
};