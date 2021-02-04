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

float f(float u,float v, float x,int n){
    if(n == 1) return v;
    else if(n == 2) return -9.8;
}

void cases(float h, float *V, float *U){
    cout << "Case when h: " << h << endl;
    float x, u, t, y, v, e, g1, g2, g3, b, b1, b2, b3;
    v = *V;
    u = *U;
    *U = u;
    *V = v;
    int a;
    x = 0;
    u = 1;
    t = 5;
    y = 45;
    cout << "Enter the 1st guess: " << endl;
    cin >> g1;
    v = g1;
    b = y;
    b1 = shoot(&U, &V, x, t, h, 1);
    cout << "b1: " << b1 << endl;
    if (fabs(b1 - b) < 0.00005){
        cout << "  The value of x and respective z are: ";
        e = shoot(&U, &V, x, t, h, 1);
        return;
    }
    else{
        cout << "Enter the 2nd guess :" << endl;
        cin >> g2;
        v = g2;
        b2 = shoot(&U, &V, x, t, h, 1);
        cout << "b2: " << b2 << endl;
    }
    if (fabs(b2 - b) < 0.00005){
        cout << "  The value of x and respective z are " << endl;
        e = shoot(&U, &V, x, t, h, 1);
        return;
    }
    else{
        cout << "g1 = " << g1 << endl;
        cout << "g2 = " << g2 << endl;
        g3 = g2 + (((g2 - g1) * (b - b2)) / (b2 - b1));
        if (b1 - b2 == 0)return;

        cout << "Exact value of guess3: " << g3 << endl;
        v = g3;
        b3 = shoot(&U, &V, x, t, h, 0);
    }
    if (fabs(b3 - b) < 0.00005){
        cout << "For which solution is :" << endl;
        e = shoot(&U, &V, x, t, h, 1);
        return;
    }
    do{
        g1 = g2;
        g2 = g3;
        b1 = b2;
        b2 = b3;
        g3 = g2 + (((g2 - g1) * (b - b2)) / (b2 - b1));
        v = g3;
        b3 = shoot(&U, &V, x, t, h, 0);

    } while (fabs(b3 - b) < 0.00005);
    v = g3;
    e = shoot(&U, &V, x, t, h, 1);
    *V = v;
    *U = u;
}

int main(){
    float V = 0;
    float U = 2;
    cases(0.001,&V, &U);
    return 0;
}
