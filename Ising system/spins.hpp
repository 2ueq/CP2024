#ifndef IsingSystem_hpp
#define IsingSystem_hpp

#include<iostream>
#include<cassert>

class IsingSpin{
    private:
        int sz;

    public:
    IsingSpin() {sz = 1;};
    ~IsingSpin() {};

    int _sz() const { return sz; };
    void set_up() {sz = 1;};
    void set_down() {sz = -1;};
    void set_sz(int sz_spec) {
        assert(sz_spec ==1 || sz_spec == -1);
        sz = sz_spec;
    };
    void flip() {sz *= -1;};
};

#endif