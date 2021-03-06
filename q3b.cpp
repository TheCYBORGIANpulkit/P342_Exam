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


int main(){
    //cout<<"The matrix M is:"<< endl;
    float arrM[12][2];
    ifstream myfileM;
    myfileM.open("M.txt");
    for(int i = 0;i<12;i++){
        for(int j = 0;j<2;j++){
            myfileM >> arrM[i][j];
            //cout<< arrM[i][j] <<" ";
        }
        //cout<<endl;
    }
    //lnw = ln(w0) - wc*t
    float w0 = 0;
    float wc = 0;
    float r_square = 0;
    least_square_exponential(arrM,&w0,&wc,&r_square);
    cout << "w0: " << w0 << endl;
    cout << "wc: " << -wc << endl;
    cout << "The fitted exponential function will be ln_w = " << w0 << " + (" << wc <<")*t." << endl;
    cout << "The Pearson's r_square is: " << r_square << endl;
    return 0;
    /*
    RESULT:
    w0: 0.790277
    wc: 0.395596
    The fitted exponential function will be ln_w = 0.790277 + (-0.395596)*t.
    The Pearson's r_square is: 0.998236
    */
}
