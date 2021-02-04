#include <iostream>  //declaring variables
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "Library.cpp"
using namespace std;

float f(float x){
    //a = sin(pi/8);
    float a = 0.3826;
    float b = pow(a,2);
    float s = pow(sin(x),2);
    float y = pow(1 - b*s,0.5);
    return 1/y;
}

int main(){
    //cout << Simpsons(f,0,1.57079632679,10) << endl;
    float factor = pow((1/9.8),0.5);
    cout << "T: " << 4*factor*Simpsons(f,0,1.57079632679,10) << endl;
    return 0;
    //RESULT  "T: 2.08728 sec"
}
