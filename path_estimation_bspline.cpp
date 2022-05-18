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
#include "tool_function.cpp"

using namespace std;

extern void LineFitLeastSquares(vector<float>, vector<float>, vector<float>&);
extern void matrix_transpose(float **, float **, int, int);
extern void *findspan(int *, int, float *, float, int, int);
extern float  **basisfunction(int *, float *, int, float, int);
extern float **bspeval(int, int **, float *, float *, int, int, int, int);
extern float CalculatePathLength(float **, int);

int main(){
    vector<float> vec_pos_A_x;
    vector<float> vec_pos_A_y;
    vector<float> vec_pos_B_x;
    vector<float> vec_pos_B_y;
    vector<float> vec_pos_ABC; 
    vector<float> vec_line_A;
    vector<float> vec_line_B;
    vector<float> mean_vel_AB;
    vector<float> vector_CB;
    vector<float> vector_CA;
    float tmp1 = 1.0/99.0;
    float k[] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};     
    int d = 2;
    int n = 2;    // number of rows of matrix c
    int m = 3;    // number of columns of matrix c
    float u[100]; //sample 100 evaluated points
    int len_u = sizeof(u)/sizeof(u[0]);
    int len_U = sizeof(k)/sizeof(k[0]);    

    int c[n][m];
    const float occlusion_length=35.0; 
    const float pi=3.14159265;
    float tan_angle_betw_AB;
    float angle_betw_AB;    
    float tan_to_angle;
    float bsp_length;
    float cal_path_length;
    float delta_length;

    for(int q=0;q<100;q++) 
    {
        u[q]=q*tmp1; 
    }
    u[99] = 1.0;

    int **cpt = (int **)malloc(n*sizeof(int*));  // pointer of matrix c
    
    // Read CSV file and save them in two vectors, p_x_f1,p_y_f1,p_x_f2,p_y_f2,...,
    ifstream fp("./data_in/GT107_picked.csv"); //Declare an ifstream object, the csv data path
    string line;
    // getline(fp,line);  // Escape the first line
    while (getline(fp,line)){  //Loop to read each row of data
        vector<float> data_line;
        string number;
        istringstream readstr(line);  
        for(int j = 0; j<line.size(); j++){ 
            getline(readstr,number,','); 
            data_line.push_back(atof(number.c_str())); //convert char array to float array 
        }
        if (data_line[2] == 761.0)   // Get point A
        {               
            vec_pos_ABC.push_back(data_line[3]);
            vec_pos_ABC.push_back(data_line[4]);
            c[0][0] = data_line[3];
            c[1][0] = data_line[4];
        }
        if (data_line[2] < 761.0) 
        {
            vec_pos_A_x.push_back(data_line[3]);  
            vec_pos_A_y.push_back(data_line[4]);
        }
        if (data_line[2] == 796.0)  // Get point B
        {
            vec_pos_ABC.push_back(data_line[3]);
            vec_pos_ABC.push_back(data_line[4]);
            c[0][1] = data_line[3];
            c[1][1] = data_line[4];
        }
        if (data_line[2] > 796.0)
        {
            vec_pos_B_x.push_back(data_line[3]);  
            vec_pos_B_y.push_back(data_line[4]);
        }
    }
    
    LineFitLeastSquares(vec_pos_A_x,vec_pos_A_y, vec_line_A);
    vec_line_A.push_back(vec_pos_ABC[1] - vec_line_A[0]*vec_pos_ABC[0]);
    LineFitLeastSquares(vec_pos_B_x,vec_pos_B_y, vec_line_B);  // 7.15655  1951.67, 3.12754  819.817  
    vec_line_B.push_back(vec_pos_ABC[3] - vec_line_B[0]*vec_pos_ABC[2]);
    /* Derive point C, enter_A_pos, enter_B_pos, mean_vel_A, mean_vel_B */
    vec_pos_ABC.push_back((vec_line_A[1]-vec_line_B[1])/(vec_line_B[0]-vec_line_A[0]));
    vec_pos_ABC.push_back((vec_line_A[0]*vec_line_B[1]-vec_line_B[0]*vec_line_A[1])/(vec_line_A[0]-vec_line_B[0]));
    c[0][2] = vec_pos_ABC[4];
    c[1][2] = vec_pos_ABC[5];
    for(int i=0; i<n; i++){
        cpt[i] = (int *) malloc(m*sizeof(int));    
        for(int j = 0; j<m; j++){  
            cpt[i][j] = c[i][j];   
        }
    }
    // Derive the mean velocity magnitude of trajectory A and B 
    mean_vel_AB.push_back(sqrt(pow((vec_pos_A_x[0]-vec_pos_ABC[0]),2)+pow((vec_pos_A_y[0]-vec_pos_ABC[1]),2))/10.0);
    mean_vel_AB.push_back(sqrt(pow((vec_pos_B_x.at(vec_pos_B_x.size()-1)-vec_pos_ABC[2]),2)+pow((vec_pos_B_y.at(vec_pos_B_y.size()-1)-vec_pos_ABC[3]),2))/10.0);

    // calculate the angle between two lines, check if it belongs to a linear/spline matching
    tan_to_angle = abs((vec_line_A[0]-vec_line_B[1])/(1+vec_line_A[0]*vec_line_B[1]));
    if(atan(tan_to_angle)/pi*180.0 < 5.0)
    {
        cout << "The angle between two lines is smaller than 10 degree."<< atan(tan_to_angle)/pi*180.0 << endl;
    }
    else{    // occluded path is approximated as a bspline to calculated the length 
        float **beva = bspeval(d, cpt, k, u, len_u, len_U, n, m);   
        float **pnts;
        pnts = (float **)malloc(len_u*sizeof(float*));
        for(int i=0;i<len_u;i++){
            pnts[i] = (float*)malloc(n*sizeof(float));
            for(int j=0;j<n;j++){
                pnts[i][j] = *(*(beva + i) + j);
                cout<<"i=" <<i <<", j=" <<j<<", "<< "pnts[i][j]="<<pnts[i][j] <<endl;   
            }
        }        
        bsp_length = CalculatePathLength(beva,len_u);
        cal_path_length = 0.5 * (mean_vel_AB[0]+mean_vel_AB[1])*occlusion_length;
        delta_length = abs(cal_path_length - bsp_length);
        cout << "spline length is " << bsp_length << endl;
        cout << "delta length is " << delta_length << endl;	 
        for(int i=0;i<n;i++) 
        {   
            free(cpt[i]);
        }   
        free(cpt);
        for(int i=0;i<len_u;i++) 
        {   
            free(beva[i]);
            free(pnts[i]);
        }   
        free(pnts);
        free(beva);
    }
     
    /* 
    Test result:
    slope, intercept of line A is 6.98239, 1599.61
    slope, intercept of line B 1.93525, 361.231.
    Mean velocity of A is 3.45352, mean velocity of B 5.31534.
    The angle between two lines is greater than 10 degree.
    Arc length is 14.7084m with an estimated error of 0.63707m.
    for bspline curve delta is 8.1479187  */
    //system("mode con:cols=100 lines=10000"); 
    system("pause");
    return 0;
}