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


extern void LineFitLeastSquares(vector<float>, vector<float>, vector<float> &);

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
    vector<float> vector_O1A;
    vector<float> vector_O1B;
    vector<float> vector_O2A;
    vector<float> vector_O2B;
    vector<float> vResult; 
    const float occlusion_length=35.0; 
    const float pi=3.14159265;
    float tan_angle_betw_AB;
    float angle_betw_AB;
    float line_bisector_slope;
    float line_bisector_intercept;
    float length_CB;
    float length_CA;    
    float tan_to_angle;
    float angle_ACB;
    float arc_length;
    float radius;
    float center_o1_x;
    float center_o1_y;
    float center_o2_x;
    float center_o2_y;
    float length_O1A;
    float length_O1B;
    float length_O2A;
    float length_O2B;
    float cos_AO1B;
    float cos_AO2B;
    float cal_path_length;
    float delta_length;

    // Read CSV file and save them in two vectors, p_x_f1,p_y_f1,p_x_f2,p_y_f2,...,
    ifstream fp("./data_in/GT107_picked.csv");  //Declare an ifstream object, the csv data path
    string line;

    while (getline(fp,line)){ // Loop to read each line
        vector<float> data_line;
        string number;
        istringstream readstr(line); 
        // Split each line by ','
        for(int j = 0; j<line.size(); j++){ 
            getline(readstr,number,','); 
            data_line.push_back(atof(number.c_str())); 
        }
        if (data_line[2] == 761.0) // f761 is point A
        {
            vec_pos_ABC.push_back(data_line[3]);
            vec_pos_ABC.push_back(data_line[4]);
        }
        if (data_line[2] < 761.0) 
        {
            vec_pos_A_x.push_back(data_line[3]);  
            vec_pos_A_y.push_back(data_line[4]);
        }
        if (data_line[2] == 796.0)  // is point B
        {
            vec_pos_ABC.push_back(data_line[3]);
            vec_pos_ABC.push_back(data_line[4]);
        }
        if (data_line[2] > 796.0)
        {
            vec_pos_B_x.push_back(data_line[3]);  
            vec_pos_B_y.push_back(data_line[4]);
        }
    }
    LineFitLeastSquares(vec_pos_A_x,vec_pos_A_y, vec_line_A);
    vec_line_A.push_back(vec_pos_ABC[1] - vec_line_A[0]*vec_pos_ABC[0]);
    LineFitLeastSquares(vec_pos_B_x,vec_pos_B_y, vec_line_B);  
    vec_line_B.push_back(vec_pos_ABC[3] - vec_line_B[0]*vec_pos_ABC[2]);
    /* Derive point C, enter_A_pos, enter_B_pos, mean_vel_A, mean_vel_B */
    // add point C
    vec_pos_ABC.push_back((vec_line_A[1]-vec_line_B[1])/(vec_line_B[0]-vec_line_A[0]));
    vec_pos_ABC.push_back((vec_line_A[0]*vec_line_B[1]-vec_line_B[0]*vec_line_A[1])/(vec_line_A[0]-vec_line_B[0]));
    // Calculate the mean velocity in trajectory A and B
    mean_vel_AB.push_back(sqrt(pow((vec_pos_A_x[0]-vec_pos_ABC[0]),2)+pow((vec_pos_A_y[0]-vec_pos_ABC[1]),2))/10.0);
    mean_vel_AB.push_back(sqrt(pow((vec_pos_B_x.at(vec_pos_B_x.size()-1)-vec_pos_ABC[2]),2)+pow((vec_pos_B_y.at(vec_pos_B_y.size()-1)-vec_pos_ABC[3]),2))/10.0);
    for (int i=0; i<vec_line_A.size(); i++)
        cout <<"slope, intercept of line A " <<vec_line_A.at(i) << ' '<<endl;
    for (int i=0; i<vec_line_B.size(); i++)
        cout <<"slope, intercept of line B " <<vec_line_B.at(i) << ' '<<endl;
    for (int i=0; i<mean_vel_AB.size(); i++)
        cout << mean_vel_AB.at(i) << ' '<<endl;
    // calculate the to angle of these two lines
    tan_to_angle = abs((vec_line_A[0]-vec_line_B[1])/(1+vec_line_A[0]*vec_line_B[1]));
    if(atan(tan_to_angle)/pi*180.0 < 10.0)
    {
        cout << "The angle between two lines is smaller than 10 degree."<< atan(tan_to_angle)/pi*180.0 << endl;
    }
    else
    {   cout << "The angle between two lines is greater than 10 degree."<< endl;
        length_CA = sqrt(pow((vec_pos_ABC[0]-vec_pos_ABC[4]),2)+pow((vec_pos_ABC[1]-vec_pos_ABC[5]),2));
        length_CB = sqrt(pow((vec_pos_ABC[2]-vec_pos_ABC[4]),2)+pow((vec_pos_ABC[3]-vec_pos_ABC[5]),2));
        vector_CA.push_back(vec_pos_ABC[0]-vec_pos_ABC[4]);
        vector_CA.push_back(vec_pos_ABC[1]-vec_pos_ABC[5]);
        vector_CB.push_back(vec_pos_ABC[2]-vec_pos_ABC[4]);
        vector_CB.push_back(vec_pos_ABC[3]-vec_pos_ABC[5]);
        tan_angle_betw_AB = (vec_line_A[0]-vec_line_B[1])/(1+vec_line_A[0]*vec_line_B[1]);
        if(atan(tan_angle_betw_AB)<0){
            angle_betw_AB = atan(tan_angle_betw_AB) + pi; // this angle is mapped to [0,pi]
        }
        else{
            angle_betw_AB = atan(tan_angle_betw_AB); // this angle is mapped to [0,pi]
        }        
        line_bisector_slope=(0.5*angle_betw_AB+vec_line_A[0])/(1-0.5*angle_betw_AB+vec_line_A[0]*vec_line_A[0]);
        line_bisector_intercept = vec_pos_ABC[5] - line_bisector_slope * vec_pos_ABC[4];
        angle_ACB = acos((vector_CA[0]*vector_CB[0]+vector_CA[1]*vector_CB[1])/length_CA/length_CB);
        radius = tan(0.5*angle_ACB)*(length_CA+length_CB)/2.0;
        /* Derive circle center O1 and O2, choose smaller angle AOB */ 
        center_o1_x = ((vec_line_A[1]-line_bisector_intercept)+radius*sqrt(1+pow(vec_line_A[0],2)))/(line_bisector_slope-vec_line_A[0]);
        center_o1_y = center_o1_x*line_bisector_slope+line_bisector_intercept;
        center_o2_x = ((vec_line_A[1]-line_bisector_intercept)-radius*sqrt(1+pow(vec_line_A[0],2)))/(line_bisector_slope-vec_line_A[0]);
        center_o2_y = center_o2_x*line_bisector_slope+line_bisector_intercept;
        vector_O1A.push_back(vec_pos_ABC[0]-center_o1_x);
        vector_O1A.push_back(vec_pos_ABC[1]-center_o1_y);
        vector_O1B.push_back(vec_pos_ABC[2]-center_o1_x);
        vector_O1B.push_back(vec_pos_ABC[3]-center_o1_y);
        vector_O2A.push_back(vec_pos_ABC[0]-center_o2_x);
        vector_O2A.push_back(vec_pos_ABC[1]-center_o2_y); 
        vector_O2B.push_back(vec_pos_ABC[2]-center_o2_x);
        vector_O2B.push_back(vec_pos_ABC[3]-center_o2_y);
        length_O1A = sqrt(pow(vec_pos_ABC[0]-center_o1_x,2) + pow(vec_pos_ABC[1]-center_o1_y,2));
        length_O1B = sqrt(pow(vec_pos_ABC[2]-center_o1_x,2)+pow(vec_pos_ABC[3]-center_o1_y,2));
        length_O2A = sqrt(pow(vec_pos_ABC[0]-center_o2_x,2)+pow(vec_pos_ABC[1]-center_o2_y,2));
        length_O2B = sqrt(pow(vec_pos_ABC[2]-center_o2_x,2)+pow(vec_pos_ABC[3]-center_o2_y,2));
        cos_AO1B = (vector_O1A[0]*vector_O1B[0]+vector_O1A[1]*vector_O1B[1])/(length_O1A*length_O1B);
        cos_AO2B = (vector_O2A[0]*vector_O2B[0]+vector_O2A[1]*vector_O2B[1])/(length_O2A*length_O2B);
        if(cos_AO1B>=cos_AO2B)
        {
            arc_length = radius * acos(cos_AO2B);
        }
        else
        {
            arc_length = radius * acos(cos_AO1B); // return value is also positive for obtuse angle, decreasing function
        }        
        // Calculated arc length 
        cal_path_length = 0.5 * (mean_vel_AB[0]+mean_vel_AB[1])*occlusion_length;
        delta_length = abs(cal_path_length - arc_length);
    }
    cout << "arc length is " << arc_length << endl;
    cout << "delta length is " << delta_length << endl;	   
    /* Evaluation result -->
    slope, intercept of line A 6.98239, 1599.61
    slope, intercept of line B 1.93525,  361.231
    Mean velocity in trajectory A 32.6km/h, in B 40.7km/h
    The angle between two lines is greater than 10 degree and it belongs to an arc matching. 
    Arc length is 14.7084m.
    Estimation error is 0.63707m. */
    //system("mode con:cols=100 lines=10000"); 
    system("pause");
    return 0;
}