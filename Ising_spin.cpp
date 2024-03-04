#include<iostream>
#include<vector>
#include"spin.hpp"
using namespace std;
int main(int argc, const char * argv[]) 
{
    IsingSpin spin;
    cout << "Spin state (initiate): " << spin._sz() << endl;

    spin.flip();
    cout << "Spin state (initiate): " << spin._sz() << endl;

    spin.flip();
    cout << "Spin state (initiate): " << spin._sz() << endl;

    const int n_spin = 10;
    vector<IsingSpin> spin_array(n_spin);

    const double J = -1;
    double energy = 0;
    for (int i = 1;i<n_spin;i++) {
        energy += spin_array[i-1]._sz() * spin_array[i]._sz();
    }
    energy += spin_array[n_spin - 1]._sz() * spin_array[0]._sz();
    energy *= J;

    cout << "Energy = " << energy << endl;

    return 0;
}