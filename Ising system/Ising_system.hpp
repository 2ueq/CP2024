#include<iostream>
#include<cmath>
// #include"IsingSpins.hpp"
#include"IsingSpinOnLattice.hpp"
#include<vector>

class IsingSystem{
protected:
    const double J;
    const int n_spins;
    const long long maxrep_state;//总状态数
    std::vector<IsingSpinOnLattice> spin;
public:
    IsingSystem(const int n_spins_spec):J(-1.0),n_spins(n_spins_spec),maxrep_state(static_cast<long long>(std::pow(2,n_spins))-1){
    spin.resize(n_spins);
};
    virtual ~IsingSystem(){};
    
    double J1()const {return J;};
    int _n_spins()const {return n_spins;};
    long long _maxrep_state()const{return maxrep_state; };
    
    int _sz(const int site_idx)const{return spin[site_idx]._sz();};
    void set_up_spin(const int site_idx){spin[site_idx].set_up();};
    void set_dw_spin(const int site_idx){spin[site_idx].set_down();};
    void set_spin(const int site_idx, int s_spec){spin[site_idx].set_sz(s_spec);};
    void flip_spin(const int site_idx){spin[site_idx].flip();};
    
    void set_state_by_code(long long rep_state){
        int rem=0;
        for(int i=0;i<n_spins;i++){
        rem=rep_state%2;
        if(rem==1){set_up_spin(i);}
        else {set_dw_spin(i);}
        if(rep_state==1){rep_state=0;}
        rep_state=floor(rep_state/2);
        };
    };  
    double eval_mz()const{
        double a=0;
        for(int t=0;t<n_spins;t++){
            a+=_sz(t);
        }
        return a;
    };
    double eval_energy_1D()const{
        double H=0;
        if(n_spins==1&&_sz(0)==1){H=1;}
        if(n_spins==1&&_sz(0)==-1){H=-1;}
        for(int i=0;i<n_spins;i++){
            if(i==n_spins-1){H+=J*_sz(i)*_sz(0);}
            else{H+=J*_sz(i)*_sz(i+1);}
            };
        return H;
    };
    
    void set_dim(int dim){
        for(IsingSpinOnLattice& each:spin)each.set_dim(dim);};
    std::vector<int> _spin_position(const int site_idx)const{
        return spin[site_idx]._position();};
    std::vector<int>_spin_NN(const int site_idx)const{return spin[site_idx]._NN();};
    int _spin_NN(const int site_idx,const int bond_idx){return spin[site_idx]._NN(bond_idx);};

};