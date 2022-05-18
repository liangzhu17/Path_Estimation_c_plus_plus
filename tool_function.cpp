#include <iostream>
#include <windows.h>
#include <stdio.h>
#include <fstream>
#include <string>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>

using namespace std;
void matrix_transpose(float **, float **, int, int);
float **bspeval(int, int **, float *, float *, int, int, int, int);

void *findspan(int *s, int n, float *u, float *U, int len_u, int len_U) // n,u,U
{
    /*
    findspan: find the knot span index

    Parameters:
    -------------------
    n - number of control points - 1
    p - spline degree
    u - parametric point
    U - knot sequence
    Returns:
    -------------------
    s - knot span index
    */    
    int temp1;   
    for(int j=0;j<len_u;j++)
    {
        if(u[j]==U[n+1])
        {
            s[j]=n;
            continue;
        }              
        for(int i=0;i<len_U;i++)
        {           
            if(U[i]<=u[j])
            {
                temp1=i;
            }            
        }
        s[j] = temp1;
    }  
}

float  **basisfunction(int *iv, float *uv, int p, float *U, int len_uv)
{
    /* 
    basicfunction:  Basis function for B-Spline

    Parameters:
    ---------------------
    iv - knot span
    uv - parametric points
    p  - spline degree
    U  - knot sequence 
    Returns
    ---------------------
    b - Basis functions vector
    */
    int i;  
    float u,temp,saved;
    float **ret; 
    float left[p+1];
    float right[p+1];
    float n[p+1];

    ret=(float **)malloc(len_uv*sizeof(float *));
	for (int l=0;l<len_uv;l++)
		ret[l]=(float *)malloc((p+1)*sizeof(float));
    
    for(int m=0;m<len_uv;m++){
        i = iv[m]+1;
        u = uv[m];        
        n[0] = 1.0;
        for(int j=1;j<p+1;j++){
            left[j] = u - U[i-j];
            right[j] = U[i+j-1] - u;           
            saved = 0.0;
            for(int r=0;r<j;r++){
                temp = n[r] / (right[r + 1] + left[j - r]);
                n[r] = saved + right[r + 1] * temp;                
                saved = left[j - r] * temp;                  
            }
            n[j] = saved;                            
        }        
        for(int k=0;k<p+1;k++){
            ret[m][k] = n[k]; 
        }    
    }
    return ret;
}

float **bspeval(int d, int **c, float *k, float *u, int len_u, int len_U, int c_row, int c_col)
{   /* 
    bspeval:  Evaluate B-Spline at parametric points
    for b-spline: The knot vector always has m+d+1 monotonically increasing knots
    Parameters:
    ---------------------
    d - Degree of the B-Spline.
    c - Control Points, matrix of size
    k - Knot sequence, row vector of size nk.
    u - Parametric evaluation points, row vector of size nu.

    Returns
    ---------------------
    pt - Evaluated points, matrix of size
    */
    float temp3;
    float **temp2,**temp5;    
    float **b,**pt,**p;
    int *s=(int *)malloc(len_u*sizeof(int));  
    findspan(s, c_col-1, u, k,len_u, len_U);

    float **bf = basisfunction(s,u,d,k,len_u);
    int *temp1=(int *)malloc(len_u*sizeof(int));  
    for(int idx=0; idx<len_u; idx++){  
        temp1[idx] = s[idx] - d + 1; 
    }    
    temp2 = (float **)malloc(c_row*sizeof(float*));
    temp5 = (float **)malloc(c_row*sizeof(float*));  

    p = (float **)malloc(c_row*sizeof(float*));
    for(int m=0;m<c_row;m++){
        p[m] = (float*)malloc(len_u*sizeof(float));  
    }
       
    pt = (float **)malloc(len_u*sizeof(float*));
    for(int m=0;m<len_u;m++){
        pt[m] = (float*)malloc(c_row*sizeof(float));   // b_row= len_u, b_col=d+1
    }  
    b = (float **)malloc(len_u*sizeof(float*));
    for(int i=0;i<len_u;i++){
        b[i] = (float*)malloc((d+1)*sizeof(float));
        for(int j=0;j<d+1;j++){
            b[i][j] = *(*(bf + i) + j);
        }
    }
    for(int i=0;i<c_row;i++){  
        for(int j=0;j<len_u;j++){
            p[i][j] = 0.0;
        }        
    }
    for(int i=0;i<d+1;i++){   // column of bf, d+1 
        for(int cr=0; cr<c_row; cr++){  // row of temp2
            temp2[cr] = (float*)malloc(len_u*sizeof(float));  
            temp5[cr] = (float*)malloc(len_u*sizeof(float)); 
            for(int j=0; j<len_u; j++){            
                temp2[cr][j] = b[j][i];
                temp5[cr][j] = c[cr][temp1[j]+i-1];  
            }   
        }                
        for(int pi=0;pi<c_row;++pi){
            for(int pj=0;pj<len_u;++pj){
                temp3 = p[pi][pj] + temp5[pi][pj]*temp2[pi][pj];
                p[pi][pj]= temp3;
            }
        }         
    }
    matrix_transpose(p, pt, c_row, len_u);        
    free(s);
    free(temp1);       
    for(int i=0;i<c_row;i++) 
    {   free(temp2[i]);
        free(temp5[i]);
    }
    free(temp2);
    free(temp5);
    for(int i=0;i<len_u;i++) 
    {   free(b[i]);
        free(bf[i]);
    }
    free(b);
    free(bf);
    free(p);
    return pt;
}

void LineFitLeastSquares(vector<float> data_x, vector<float> data_y, vector<float> &vResult)  
{   /*
    Least Square methods to determin the line of best fit for a set of data
    Parameters:
    -----------------------
    data_x - Values of x-axis from the input points
    data_y - Values of y-axis from the input points
    Returns:
    -----------------------
    vResult - Derived slope and intercept of the line of best fit
    */
    float A = 0.0;
    float B = 0.0;
    float C = 0.0;
    float D = 0.0;
    float k,b = 0.0;
    float data_n = data_x.size();

    for (int i=0;i<data_n;i++)
    {
        A += data_x[i] * data_x[i];
        B += data_x[i];
        C += data_x[i] * data_y[i];
        D += data_y[i];
    }
    
    if (data_n*A-B*B!=0)
    {
        k = (data_n*C-B*D)/(data_n*A-B*B);
        b = (A*D-B*C)/(data_n*A-B*B);
    }
    else
    {
        k = 10.0e8;
        b = 0;
    }
    vResult.push_back(k);
}

float CalculatePathLength(float **p, int pt_num){ 
    /*
    Calculate the approximated spline length using given pt_num points
    Parameters:
    -----------------------
    p - Evaluated points from function bspeval
    pt_num - The number of evaluated points
    Returns:
    -----------------------
    length - The derived curve length
    */
    float length=0.0;
    float tmp=0.0;     
   
    for (int j = 0; j < pt_num-1; j++)
    {
        tmp = sqrt(pow(p[j][0]-p[j+1][0],2) + pow(p[j][1]-p[j+1][1],2));
        length = tmp + length;
    } 
    return length;
}
void matrix_transpose(float **matrix, float **matrix_t, int n, int m)
{   // matrix transpose function, transposed matrix matrix_t
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<m; j++)
        {   
            matrix_t[j][i] = matrix[i][j];
        }
    }
}


/* 
Test the functions **bseval(), *findspan(), **basisfunction
int main(){

    float tmp1 = 1.0/3.0;
    float tmp2 = 1.0/9.0;
    float k[] = {0.0, 0.0, 0.0, 0.0, 0.3, 0.6, 1.0, 1.0, 1.0, 1.0};
    k[4] = tmp1;
    k[5] = tmp1*2;
    int d = 3;
    int n = 2;  
    int m = 6;   
    float u[10];
    int **ct,**c;
    int arr[] = {0,10,20,30,40,50,0,5,-5,5,-5,0};
    for(int q=0;q<10;q++) 
    {
        u[q]=q*tmp2; 
    }
    u[9] = 1.0;
    int len_u = sizeof(u)/sizeof(u[0]);
    int len_U = sizeof(k)/sizeof(k[0]);
    
    ct = (int **)malloc(m*sizeof(int*));
    for(int j=0;j<m;j++){
        ct[j] = (int *) malloc(n*sizeof(int));
    }     
    c = (int **)malloc(n*sizeof(int*));
    for(int j=0;j<n;j++){
        c[j] = (int *) malloc(m*sizeof(int));
    }
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++)
            c[i][j] = arr[i*m+j];
    }   

    float **beva = bspeval(d, c, k, u, len_u, len_U, n,m);    
    for(int i= 0;i<len_u;i++){
        for(int j=0; j<n; j++){
            //cout<< "Values in pt are "<< *((float *)beva+n*i+j)<<endl;  
            cout<<"i=" <<i <<", j=" <<j<< endl;
            cout << beva[i][j]<<endl;
        }           
    }
    system("mode con:cols=100 lines=10000");  // Set number of printed rows/ columns 
    free(ct);
    free(c);
    free(beva);
    
    system("pause");
    return 0;
}*/

