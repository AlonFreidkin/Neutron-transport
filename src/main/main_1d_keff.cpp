#include <iostream>
#include "../../lib/sparsex/src/ellpack.h"
#include "../../lib/sparsex/src/algebra.h"
int main(){
    size_t npoints = 100;
    double d = 1.0;
    double length = 10.0;
    double h = length/(npoints-1);
    double xs_a = 1.0;
    double xs_f = 1.0;
    double nu = 2.0;
    Ellpack A(npoints,npoints,3);
    for (int i = 0;i<npoints;i++){
        double a = -d/(h*h);
        double b = 2*d/(h*h)+xs_a;
        A.insert(i,i-1,a);
        A.insert(i,i,b);
        A.insert(i,i+1,a);
    }
    A.deleteRow(0);
    A.insert(0, 0, 1.0);
    A.deleteRow(npoints-1);
    A.insert(npoints-1, npoints-1, 1.0);
    
    Ellpack F(npoints,npoints,3);
    for (int i = 0;i<npoints;i++){
        double a = nu*xs_f;
        F.insert(i,i,a);
        
    }
    F.deleteRow(0);
    F.insert(0, 0, 1.0);
    F.deleteRow(npoints-1);
    F.insert(npoints-1, npoints-1, 1.0);

    std::vector<double> source(npoints);
    std::vector<double> source_new(npoints);
    std::vector<double> phi(npoints,1.0);
    double keff = 1.0;
    while (true)
    {
        /* code */
        phi[0] = 0;
        phi[npoints-1] = 0;
    
    F.mvp(source,phi);
    for (auto &val : source){
        val /= keff;
    }
    
    A.solve_cg(phi,source);
    F.mvp(source_new, phi);
    
    double power_new = 0, power = 0;
    for (int i = 0;i<npoints;i++){
        power_new += source_new[i];
        power += source[i];
    }
    double keff_new = power_new/power;
    if (std::abs(keff-keff_new)/keff<0.0001) break;
    keff = keff_new;
    }
    std::cout<<"keff: " <<keff<<std::endl;
    for (int i = 0;i<npoints;i++){
        std::cout<<i*h<<" "<<phi[i]<<std::endl;
    }
    
    std::cout<<std::endl;
    return 0;
}
