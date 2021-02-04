#include <iostream>  //declaring variables
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdlib>
#include <cmath>
using namespace std;

//least square fitting for an exponential curve
void least_square_exponential(float arr[12][2], float*A, float*B, float *R){
    float sqx = 0;
    for(int i = 0;i<12;i++){
        sqx = sqx + pow(arr[i][0],2);
    }

    float sqy = 0;
    for(int i = 0;i<12;i++){
        sqy = sqy + pow(log(arr[i][1]),2);
    }

    float sx = 0;
    for(int i = 0;i<12;i++){
        sx = sx + arr[i][0];
    }

    float denom = 12*sqx - pow(sx,2); //denominator in the formula of a and b


    float sy = 0;
    for(int i = 0;i<12;i++){
        sy = sy + log(arr[i][1]);
    }

    float sxy = 0;
    for(int i = 0;i<12;i++){
        sxy = sxy + (arr[i][0])*(log(arr[i][1]));
    }

    ////numerator of a
    float num_a = sy*sqx - sx*sxy;
    float a = *A;
    a = num_a/denom;
    *A = a;

    //numerator of b
    float num_b = 12*sxy - sx*sy;
    float b = *B;
    b = num_b/denom;
    *B = b;
    //calculating S_xx and S_yy
    float S_xx = sqx - pow(sx,2)/12;
    float S_yy = sqy - pow(sy,2)/12;
    //calculating S_xy
    float S_xy = sxy - sx*sy/12;
    //calculating r
    float r_square = *R;
    r_square = pow(S_xy,2)/(S_xx*S_yy);
    *R = r_square;
}


//straight line least square fitting
void least_square_straight_line(float arr[12][2], float*A, float*B, float *R){
    float sqx = 0;
    for(int i = 0;i<12;i++){
        sqx = sqx + pow(arr[i][0],2);
    }

    float sqy = 0;
    for(int i = 0;i<12;i++){
        sqy = sqy + pow(arr[i][1],2);
    }

    float sx = 0;
    for(int i = 0;i<12;i++){
        sx = sx + arr[i][0];
    }

    float denom = 12*sqx - pow(sx,2); //denominator in the formula of a and b


    float sy = 0;
    for(int i = 0;i<12;i++){
        sy = sy + arr[i][1];
    }

    float sxy = 0;
    for(int i = 0;i<12;i++){
        sxy = sxy + (arr[i][0])*(arr[i][1]);
    }

    ////numerator of a
    float num_a = sy*sqx - sx*sxy;
    float a = *A;
    a = num_a/denom;
    *A = a;

    //numerator of b
    float num_b = 12*sxy - sx*sy;
    float b = *B;
    b = num_b/denom;
    *B = b;
    //calculating S_xx and S_yy
    float S_xx = sqx - pow(sx,2)/12;
    float S_yy = sqy - pow(sy,2)/12;
    //calculating S_xy
    float S_xy = sxy - sx*sy/12;
    //calculating r
    float r_square = *R;
    r_square = pow(S_xy,2)/(S_xx*S_yy);
    *R = r_square;
}

//monte carlo in 3D
double monte_carlo_volume(double (*fn)(double,double,double),double a,double b,double c,int N){
    float counter = 0;
    double arrX[N][3];
    ofstream outfile;
        outfile.open("monte_carlo_ellipsoid.csv");
    for(int i = 0;i<N;i++){
        double p = -a + (float(rand())/float(RAND_MAX))*2*a;
        double q = -b + (float(rand())/float(RAND_MAX))*2*b;
        double r = -c + (float(rand())/float(RAND_MAX))*2*c;
        arrX[i][0] = p;
        arrX[i][1] = q;
        arrX[i][2] = r;
        double f = (*fn)(p,q,r);
        if(f < 0){
            outfile<< arrX[i][0] << "," << arrX[i][1] << "," << arrX[i][2]<<  endl;
            counter++;
        }
    }
    outfile.close();
    double x = counter/N;
    double V = x*a*b*c*8;

    return V;
}

void random_walk(double*x, double*y, double* r){
    //starting point
    double X = *x;
    double Y = *y;
    double R = *r;
    //generate a random angle
        R = sqrt(pow(X,2) + pow(Y,2));
        double a = (float(rand())/float(RAND_MAX))*2*3.141592;
        X = X + cos(a);
        Y =  Y + sin(a);
    *x = X;
    *y = Y;
    *r = R;
}

//monte carlo in 1Dimension
double Monte_Carlo(double (*fn)(float), float a, float b, int N){
    double h = (b - a)/N;
    double rand_n;
    double arr[N];
    for (int i = 0; i < N; i++){
      rand_n = (float(rand())/float((RAND_MAX)) * b);
      arr[i]= rand_n;
    }
    double s1 = 0;
    for(int i= 0;i<N;i++){
        double x = pow((*fn)(arr[i]),2);
        s1 = s1 + x;
    }
    double A = s1/N;

    double s2 = 0;
    for(int i= 0;i<N;i++){
        s2 = s2 + (*fn)(arr[i]);
    }
    double B = (pow(s2/N,2));
    double F = h*(s2);
    double sigma = pow((A - B),0.5);
    //cout<<  sigma  << " "<< F << endl;
    return F;
}
//RK4 for only 1st order
void RK4_1storder(float*u, float x, float t, float h,float(*fn)(float, float) )
{
    float U = *u;
	float k11, k12, k13, k14;
	k11 = h * (*fn)(U, x);

	k12 = h * (*fn)(U + 0.5 * k11, x + 0.5 * h);

	k13 = h * (*fn)(U + 0.5 * k12, x + 0.5 * h);

	k14 = h * (*fn)(U + k13, x + h);

	U += (k11 + 2 * k12 + 2 * k13 + k14) / 6;
	*u = U;
}

//RK4 for a system of 1st prder DE converted from 2nd order DE

//converting 2nd order DE in system of 1st order DE
float f(float u,float v, float x,int n){
    if(n == 1) return v;
    else if(n == 2) return -9.8;
}
void RK4(float *u, float *v, float x, float t, float h, float p)
{
    float U = *u;
    float V = *v;
	float k11, k12, k13, k14;
	float k21, k22, k23, k24;
	k11 = h * f(U,V, x, 1);
	k21 = h * f(U, V, x, 2);
	k12 = h * f(U + 0.5 * k11, V + 0.5 * k21, x + 0.5 * h, 1);
	k22 = h * f(U + 0.5 * k11, V + 0.5 * k21, x + 0.5 * h, 2);
	k13 = h * f(U + 0.5 * k12, V + 0.5 * k22, x + 0.5 * h, 1);
	k23 = h * f(U + 0.5 * k12, V + 0.5 * k22, x + 0.5 * h, 2);
	k14 = h * f(U + k13, V + k23, x + h, 1);
	k24 = h * f(U + k13, V + k23, x + h, 2);
	U += (k11 + 2 * k12 + 2 * k13 + k14) / 6;
	V += (k21 + 2 * k22 + 2 * k23 + k24) / 6;
	*u = U;
	*v = V;
}
float shoot(float *u, float *v, float x, float t, float h, int a){
	ofstream outfile;
		outfile.open("shoot4.csv");

	float y, guess,U,V;
	y = *u;
	U = y;
	guess = *v;
	V = guess;
    cout << "x     u" << endl;
	while (abs(x) < abs(t)){
		RK4(&U, &V, x, t, h,2);
		U = U + y;
		V = V + guess;
		x = x + h;
		if (a == 1){
			cout <<  x << "    " << U << endl;
			outfile << t << "," << x <<  "," << U << endl;
		}
	}
	*v = U;
	*u = V;
	return U;
}

void Explicit_Euler(float(*f)(float,float),float xi,float xf,float y, int N){
    float h = (xf - xi)/N;
    float x = xi;
    ofstream outfile;
    outfile.open("explicit1.csv");
    for(int i = 0;i<N;i++){
        outfile << x << "," << y << endl;
        float K1 = h*((*f)(x + h,y));
        y = y + K1;
        x = x + h;
    }
    outfile.close();
}

float Simpsons(float (*fn)(float), float a, float b, int N){
    float h = (b - a)/N;
    float s = 0;
    for(int i = 0;i <= N;i++){
        float  x = (h/3)*((*fn)(a + i*h));
        if(i == 0 || i == N)s = s + x;
        else if(i%2 == 0)s = s + 2*x;
        else s = s + 4*x;
    }
    return s;
}

double Trapezoidal(double (*fn)(float), float a, float b, int N){
    double h = (b - a)/N;
    double s = 0;
    for(int i = 0;i<=N;i++){
        double x = a + i*h;
        double y = (h/2)*((*fn)(x));
        if(i == 0 || i == N)s = s + y;
        else s = s+2*y;
    }
    return s;
}

double Midpoint(double (*fn)(float),float a, float b, int N ){
    double h = (b - a)/N;
    double s = 0;
    for(int i = 0;i<= N;i++){
        double x = (a+i*h + a + (i+1)*h)/2;
        s = s + h*(*fn)(x);
        //cout << s << endl;
    }
    return s;
}

void Bracketing(float (*fn)(float),double a, double b){
    int Count = 0;
    double beta = 0.75;
    for(int j = 0;j< 200;j++){
        if((*fn)(a)*(*fn)(b) > 0){
            if(abs((*fn)(a)) > abs((*fn)(b))){
                b = b + beta*(b - a);
            }
            else a = a - beta*(b - a);
        }
    }
}
double Bisection(float (*fn)(float), double a, double b, double*C){
    Bracketing((*fn),a,b);
    int Count = 0;
    double arrC[200];
    double c = *C;
    //cout << "c:" << c << endl;
    for(int j = 0;j<200;j++){
        double c = (a+b)/2;
        if(abs(a - b) > 0.000001){
            arrC[j] = abs(a - b);
            if((*fn)(a)*(*fn)(c) < 0) b = c;
            else  a = c;
            //Count ++;
        }
        else{
            arrC[j] = abs(a - b);
            ofstream outfile;
            outfile.open("bisection.csv");
            for(int i = 0;i<= j; i++){
                outfile<< arrC[i] << endl;
            }
            outfile.close();
            *C = c;
            //
            return (*fn)(c);
        }
    }
}

double False_position(float (*fn)(float),float a, float b, double *C){
    Bracketing((*fn),a,b);
    double arrF[200];
    double c = *C;
    //cout << "c:" << c << endl;
    for(int j = 0;j<200;j++){
        float c = b - (((b-a)*((*fn)(b)))/((*fn)(b) - (*fn)(a)));
        if(abs((*fn)(c)) > 0.000001){
            arrF[j] = abs((*fn)(c));
            if((*fn)(a)*(*fn)(c) < 0) b = c;
            else a = c;
            //Count ++;
        }
        else{
            arrF[j] = abs((*fn)(c));
            ofstream outfile;
            outfile.open("falsi.csv");
            for(int i =0;i<j;i++){
                outfile<< arrF[j] << endl;
            }
            outfile.close();
            *C = c;
            return (*fn)(c);
        }
    }
}

//Newton_raphson should be called in main function like newton_raphson(f,&guess);
float newton_raphson(float (*fn)(float),float *guess){
    float h = 0.001;
    float x = *guess;
    float arrX[200];

    for(int j = 0;j<200;j++){
        float delf = ((*fn)(x + h) - (*fn)(x - h))/(2*h);
        float a = (*fn)(x)/delf;
        float temp = x - a;
        if(abs(x - temp) > 0.0001){
            x = temp;
            arrX[j] = x - temp;
        }
        else{
            *guess = x; // root
            //for storing the error file
            arrX[j] = x - temp;
            ofstream outfile;
            outfile.open("newton_raphson.csv");
            for(int i = 0;i<= j; i++){
                outfile<< arrX[i] << endl;
            }
            outfile.close();
            //cout << "Please find the attached csv file generated and the plot" << endl;
            return (*fn)(x); //getting the value of function at the root
        }
    }
}

void Part_Pivoting(float arrA[4][4]){
    for(int r = 0;r<3;r++){
        float pivot = arrA[r][r];
        //cout << pivot << endl;
        if(pivot == 0){
           for(int i = r+1;i<4;i++){
                if(abs(int(arrA[i][r])) > abs(int(pivot))){
                    float temp[4];
                    int j = 0;
                    for(int j = 0;j<4;j++){
                        temp[j] = arrA[r][j];
                        arrA[r][j] = arrA[i][j];
                        arrA[i][j] = temp[j];
                    }
                }
                else continue;
           }
        }
        else continue;
    }
}

void Gauss_Jordan(float arrA[3][4]){
    Part_Pivoting(arrA);
    float pivot;
    for(int r = 0;r<3;r++){
        //deciding a pivot
        pivot = arrA[r][r];
        //moving to columns
        for(int c = r;c<4;c++){
            arrA[r][c] = (arrA[r][c])/pivot;
        }
        for(int k = 0;k<3;k++){
            if(k == r || arrA[k][r] == 0){
                continue;
            }
            else{
                float factor = arrA[k][r];
                for(int c = r;c<4;c++){
                   arrA[k][c] = arrA[k][c] - factor*arrA[r][c];
                }
            }
        }
    }
    //make into a new matrix
    float arrX[3][1];
    cout<< "The solution are: " << endl;
    for(int r = 0;r<3;r++){
        arrX[r][0] = arrA[r][3];
        //cout<< arrX[r][0];
        //cout<< endl;
    }

    float x = arrX[0][0];
    float y = arrX[1][0];
    float z = arrX[2][0];
    /*

    cout << "The inverse matrix is: " << endl;
    int** arrN = new int*[3];
    for(int r = 0;r<3;r++){
        arrN[r] = new int[3];
        for(int j = 0;j<3;j++){
            arrN[r][j] = arrA[r][j+3];
            cout<< arrN[r][j]<< " ";
        }
        cout<< endl;
    }
    return arrN;
    */
}

void LU_decomposition(float arrA[4][4]){
    //Partial pivoting at first
    Part_Pivoting(arrA);
    for(int j = 0;j<4;j++){
        if(j == 0){
            //u[0][0] = arrA[0][0];
            //cout<< arrA[0][0] << endl;
            for(int i = 1;i<4;i++){
                arrA[i][j] = arrA[i][j]/arrA[0][0];
                //cout<< arrA[i][j] << endl;
            }
        }
        else{
            for(int i =0;i<4;i++){
                float s = 0;
                float a = 0;
                for(int k = 0;k<i;k++){
                    a = arrA[i][k]*arrA[k][j];
                    s = s + a;
                }
                //cout<< s << endl;
                if(i < j)arrA[i][j] = arrA[i][j] - s;
                else if(i == j){
                    arrA[i][i] = arrA[i][j] - s;
                    //if(i == j)cout<< arrA[i][j] << endl;
                }
                else{
                    float s = 0;
                    float a = 0;
                    for(int k = 0;k<j;k++){
                        a = arrA[i][k]*arrA[k][j];
                        s = s + a;
                    }
                    arrA[i][j] = (arrA[i][j] - s)/arrA[j][j];
                }
            }
        }
    }
}
float Determinant(float arrA[4][4]){
    //Pivoting(arrA);
    LU_decomposition(arrA);
    float p = 1;
    for(int i =0;i<4;i++){
        p = p*arrA[i][i];
    }
    return p;
}
// NOTE: This code assumes that the matrix A is already LU decomposed.
void FB_Substitution(float arrA[4][4],float arrb[4][4]){
    //RHS matrix is arrb[4][n]
    //Making augmented matrix
    int n = 4;
    int m = 4 + n;
    float arrL[4][m];
    for(int i = 0;i<4;i++){
        for(int j = 0;j<4;j++){
            if(i > j) arrL[i][j] = arrA[i][j];
            if(i == j) arrL[i][j] = 1;
            if(i < j) arrL[i][j] = 0;
        }
    }
    float arry[4][n];
    //solving Y
    if(n == 1){
        for(int i = 0;i<4;i++){
            float s = 0;
            for(int j = 0;j<i;j++){
                s = s + arrL[i][j]*arry[j][0];
            }
            arry[i][0] = (arrb[i][0] - s)/arrL[i][i];
        }
        //for backward substitution
        float arrU[4][4];
        for(int i = 0;i<4;i++){
            for(int j = 0;j<4;j++){
                if(i <= j) arrU[i][j] = arrA[i][j];
                if(i > j) arrU[i][j] = 0;
            }
        }
        //solving X
        float arrX[4][1];
        for(int i = 3;i >= 0;i--){
            float s = 0;
            for(int j = i+1;j<4;j++){
                s = s + arrU[i][j]*arrX[j][0];
            }
            //cout<< s << endl;
            arrX[i][0] = (arry[i][0] - s)/arrU[i][i];
        }
    }
    else{
    //solving Y
        for(int i = 0;i<4;i++){
            for(int j = 0;j<n;j++){
                float s = 0;
                for(int k = 0;k<i;k++){
                    s = s + arrL[i][k]*arry[k][j];
                }
                arry[i][j] = (arrb[i][j] - s)/arrL[i][i];
                //cout<< arry[i][j] << " ";
            }
            //cout<< endl;
        }
        //for backward substitution
        float arrU[4][4];
        for(int i = 0;i<4;i++){
            for(int j = 0;j<4;j++){
                if(i < j)arrU[i][j] = arrA[i][j];
                else if(i == j)arrU[i][j] = arrA[i][j];
                else arrU[i][j] = 0;
            }
        }
        //solving X[4][4]
        float arrX[4][4];
        for(int i = 3;i >= 0;i--){
           for(int j = 0;j<n;j++){
                float s = 0;
                for(int k = i+1;k<n;k++){
                    s = s + arrU[i][k]*arrX[k][j];
                }
                arrX[i][j] = (arry[i][j] - s)/arrU[i][i];
            }
        }
        //The INVERSE is arrX[i][j]

    }
}

void Matrix_Product(float arrM[4][4], float arrN[4][4],int n){
    float arrP[n][n], x;
    //int arrQ[3][3];
    for(int i = 0;i<n;i++){
        for(int j = 0;j<n;j++){
            arrP[i][j] = 0;
            //arrQ[i][j] = 0;
        }
    }
    for(int a = 0;a<n;a++){
        for(int b = 0;b<n;b++){
            for(int i = 0;i<n;i++){
                x = float(arrM[a][i])*float(arrN[i][b]);
                //cout<<x<<endl;
                arrP[a][b] = arrP[a][b] + x;
                //cout<< arrP[a][b];
            }
            //cout<< arrP[a][b]<< " ";
            //cout<< endl;
        }
    }

}





