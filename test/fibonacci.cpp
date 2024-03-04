#include <iostream>
using namespace std;

int fibonacci(int n){
    int s=0;
    if (n<=0) {return 0;}
    if (n==1) {return 1;}
    else {
        s=fibonacci(n-1)+fibonacci(n-2);
        return s;
    }
}

int main(){
    cout<<fibonacci(7)<<endl;
    cout<<fibonacci(15)<<endl;
    cout<<fibonacci(37)<<endl ;
    return 0;
}


// int main(){
//     int n;
//     cin>>n;
//     cout<<fibonacci(n)<<endl;
//     return 0;
// }