#ifndef IsingSystem_Square_hpp
#define IsingSystem_Square_hpp

#include"Ising_system.hpp"
#include<vector>
#include<iostream>
#include<algorithm> 

class IsingSystem_Square:public IsingSystem{
    private:
        std::vector<int>system_size;
        std::vector<double>beta;
        std::vector<double>E;
        std::vector<double>E_sq;
        std::vector<double>Z;
        std::vector<double>M_sq;
        std::vector<double>C;
        double ground_state_energy;
        
//计算每个粒子四周的位置信息
        void setup_NN(){
            for(int site_idx=0;site_idx<n_spins;site_idx++){
                std::vector<int>r=_lattice_coordinate(site_idx);
                spin[site_idx].set_NN(0,site_index(shift_pos_x(r)));
                spin[site_idx].set_NN(1,site_index(shift_pos_y(r)));
                spin[site_idx].set_NN(2,site_index(shift_neg_x(r)));
                spin[site_idx].set_NN(3,site_index(shift_neg_y(r)));
            };
        };
    public:
        IsingSystem_Square(const std::vector<int>system_size_spec,std::vector<double>beta_spec={0}):
            IsingSystem(system_size_spec[0]*system_size_spec[1]),
            system_size(system_size_spec),beta(beta_spec){
                setup_NN();
                E.assign(beta.size(),0);
                E_sq.assign(beta.size(),0);
                M_sq.assign(beta.size(),0);
                Z.assign(beta.size(),0);
                C.assign(beta.size(),0);
                ground_state_energy = eval_energy();
            };
        ~IsingSystem_Square(){};
//位置和坐标转换函数
    int site_index(const std::vector<int>lattice_coordinate)const{
        return lattice_coordinate[1]*system_size[0]+lattice_coordinate[0];
    };
    std::vector<int>_lattice_coordinate(int site_index)const{
        std::vector<int>a(system_size);
        a[1]=site_index/system_size[0];
        a[0]=site_index-a[1]*system_size[0];
        return a;
    };

//移动坐标    
    std::vector<int>shift_pos_x(const std::vector<int>r_spec){
        std::vector<int>r(r_spec);
        r[0]=(r[0]+1)%system_size[0];
        return r;
    };
    std::vector<int>shift_pos_y(const std::vector<int>r_spec){
        std::vector<int>r(r_spec);
        r[1]=(r[1]+1)%system_size[1];
        return r;
    };
    std::vector<int>shift_neg_x(const std::vector<int>r_spec){
        std::vector<int>r(r_spec);
        r[0]=(r[0]-1+system_size[0])%system_size[0];
        return r;
    };
    std::vector<int>shift_neg_y(const std::vector<int>r_spec){
        std::vector<int>r(r_spec);
        r[1]=(r[1]-1+system_size[1])%system_size[1];
        return r;
    };

//给出四周粒子的位置信息
    int NN(const int site_idx,const int bond_idx){
        return spin[site_idx]._NN(bond_idx);
    };

//初始化+计算M和E
    void set_state(const std::vector<bool>state){
        int b=state.size();
        if (b==n_spins) {
            for(int i = 0;i<n_spins;++i) {
                if(state[i]==1) set_up_spin(i);
                else set_dw_spin(i);
            };
        };
    };    
    std::vector<bool>state_by_code(const long long& rep_states){
        int rem=0;
        std::vector<bool>state;
        state.assign(n_spins,0);
        long long rep_state=rep_states;
        for(int i=0;i<n_spins;i++){
            rem=rep_state%2;
            if(rem==1){state[i]=1;}
            else {state[i]=0;}
            if(rep_state==1){rep_state=0;}
            rep_state=floor(rep_state/2);
        };
        return state;
    };
    double eval_Mz()const{
        double a=0;
        for(int t=0;t<n_spins;t++){
            a+=_sz(t);
        };
        return a;
    };
    double eval_energy(){
        double E=0;
        for(int t=0;t<n_spins;t++){
            for(int i=0;i<=3;i++){
                E+=spin[t]._sz()*spin[spin[t]._NN(i)]._sz();
            };
        };
        E*=J/2;
        return E;
    };

//计算玻尔兹曼权函数
    double weight_unnormalized(const std::size_t beta_idx){
        return exp(-1*beta[beta_idx]*eval_energy());
    };

//计算相关参数
    void exactly_evaluate_given(){
        for(std::size_t i=0;i<beta.size();i++){
            E[i]+=_exact_energy_q(i);
            E_sq[i]+=_exact_energy_q_sq(i);
            Z[i]+=weight_unnormalized(i);
            M_sq[i]+=_exact_magz_q_sq(i);
        };
        for(std::size_t i=0;i<beta.size();i++)
        {
            E[i]*=1/Z[i];
            E_sq[i]*= 1/Z[i];
            M_sq[i]*=1/Z[i];
        };
        for(std::size_t i = 0;i<beta.size();i++)
        {
            C[i]=beta[i]*beta[i]*(E_sq[i]-E[i]*E[i]);
        };
    };
    void exactly_evaluate(const std::vector<bool>& state){
        set_state(state);
        exactly_evaluate_given();
    };
    void exactly_evaluate(const long long& rep_state){
        std::vector<bool>state=state_by_code(rep_state);
        exactly_evaluate(state);
    };

    void exact(){
        long long rep_state=0;
        while(rep_state<=maxrep_state){exactly_evaluate(rep_state++);};
        normalize_direct();
    };
    void normalize_direct(){
        double min_C=*std::min_element(C.begin(), C.end());
        double max_C=*std::max_element(C.begin(), C.end());
        double min_M_sq=*std::min_element(M_sq.begin(), M_sq.end());
        double max_M_sq=*std::max_element(M_sq.begin(), M_sq.end());
        for(std::size_t i=0;i<C.size();i++){
            C[i] = (C[i] - min_C) / (max_C - min_C);
            M_sq[i] = (M_sq[i] - min_M_sq) / (max_M_sq - min_M_sq);
        };
    };
    void print_exact()const{
        std::cout<<"Specific Heat:";
        for(double value:C){ 
            std::cout<<value<< " ";
        };
        std::cout<<"."<<std::endl;
        std::cout<<"Magnetization (Squared): ";
        for(double value:M_sq){
            std::cout<<value<<" ";
        };
        std::cout<<"."<<std::endl;
    };

    double _exact_energy_Z(std::size_t beta_idx){
        return weight_unnormalized(beta_idx);
    };
    double _exact_energy_q(std::size_t beta_idx){
        return eval_energy()*_exact_energy_Z(beta_idx);
    };
    double _exact_energy_q_sq(std::size_t beta_idx){
        return eval_energy()*eval_energy()*_exact_energy_Z(beta_idx);
    };
    double _exact_magz_Z(std::size_t beta_idx){
        return weight_unnormalized(beta_idx);
    };
    double _exact_magz_q_sq(std::size_t beta_idx){
        return eval_Mz()*eval_Mz()*_exact_magz_Z(beta_idx);
    };
    double ground_state()const{return ground_state_energy;};
};
#endif