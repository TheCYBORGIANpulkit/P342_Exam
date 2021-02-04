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

//equation
float f(float x){
    float y = exp(x);
    return (x - 5)*y + 5;
}

int main(){
    float guess = 6; //
    cout << "The value of equation at the root is " << newton_raphson(f,&guess) << endl;
    cout << "The root of the equation is: " << guess << endl;
    //finding wein's constant b
    float b = (19.878*pow(10,-26))/(1.381*pow(10,-23)*guess);
    cout << " Value of Wein's constant b is : " << b << endl;
    return 0;
    /*RESULT :
    The value of equation at the root is 0.000115966
    The root of the equation is: 4.96512
    Value of Wein's constant b is : 0.00289901
 */
}
