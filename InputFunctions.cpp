#include "InputFunctions.h"
using namespace std;

double FunctionF(PointD p) {
    double x = p.first;
    double y = p.second; 
    return 2*(x*x+y*y)*(1-2*x*x*y*y)*exp(1-x*x*y*y);
}
double FunctionPhi(PointD p){
    double x = p.first;
    double y = p.second;
    return exp(1-x*x*y*y);
}
