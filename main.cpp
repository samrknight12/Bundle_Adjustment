#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <Eigen/Dense>
#include <cmath>

using namespace std;
using namespace Eigen;

void getArraySize(string FileName, int &rows, int &cols){
    string line,item;
    vector<float>data;
    rows = 0, cols =0;

    ifstream file(FileName);
    while (getline(file, line))
    {
        rows++;
        if (rows == 1)
        {
            stringstream ss(line);
            while (ss >> item) cols++;
        }
    }
    file.close();
}


void readData(MatrixXd& data, string filename)
{
    ifstream file;
    file.open(filename);
    for (int i = 0; i<data.rows(); i++){
        for (int j = 0; j<data.cols(); j++){
            file >> data(i,j);
        }
    }
}

int findEOPIndex(MatrixXd EOP, float image_id){

    int index;
    for (int i=0; i<EOP.rows();i++){
        if (image_id == EOP(i,0)){
                index = i;
                return index;
            }
    }
}

int findIOPIndex(int camera_ID, MatrixXd IOP){
    int index;
    for (int i=0;i<IOP.rows();i++){
        if (camera_ID == IOP(i,0)){
            index = i;
            return index;
        }
    }
}

bool controlOrTiePoint(MatrixXd control_points, int point_id){
    bool return_val;
    for (int i=0;i<control_points.rows();i++){
        if (point_id == control_points(i,0)){
            return(true);
        }
        else{
            return_val = false;
        }
    }
    return(return_val);
}

int findControlIndex(int point_id, MatrixXd control_points){
    for (int i=0;i<control_points.rows();i++){
        if (point_id == control_points(i,0)){
            return i;
        }
    }
}

int findtieIndex(int point_id, MatrixXd tie_points){
    for (int i=0;i<tie_points.rows();i++){
        if (point_id == tie_points(i,0)){
            return i;
        }
    }
}

void rotationMatrix(MatrixXd& M, MatrixXd EOP, unsigned int cur_eop_idx){
    float degree_2_rad = M_PI / 180;
    float w = EOP(cur_eop_idx,5)*degree_2_rad;
    float phi = EOP(cur_eop_idx,6)*degree_2_rad;
    float kappa = EOP(cur_eop_idx,7)*degree_2_rad;

    M(0,0) = cos(phi)*cos(kappa);
    M(0,1) = (cos(w)*sin(w)) + (sin(w)*sin(phi)*cos(kappa));
    M(0,2) = (sin(w)*sin(kappa)) - (cos(w)*sin(phi)*cos(kappa));

    M(1,0) = -1*cos(phi)*sin(kappa);
    M(1,1) = (cos(w)*cos(kappa)) - (sin(w)*sin(phi)*sin(kappa));
    M(1,2) = (sin(w)*cos(kappa)) + (cos(w)*sin(phi)*sin(kappa));

    M(2,0) = sin(phi);
    M(2,1) = -1*sin(w)*cos(phi);
    M(2,2) = cos(w)*cos(phi);
}

float uNumber(MatrixXd M, MatrixXd ground_point, MatrixXd EOP, unsigned int cur_point_idx, unsigned int cur_eop_idx){
    float U = (M(0,0)*(ground_point(cur_point_idx,1)-EOP(cur_eop_idx,2))) + (M(0,1)* (ground_point(cur_point_idx,2)-EOP(cur_eop_idx,3))) + (M(0,2)* (ground_point(cur_point_idx,3)-EOP(cur_eop_idx,4)));
    return U;
}

float vNumber(MatrixXd M, MatrixXd ground_point, MatrixXd EOP,unsigned int cur_point_idx, unsigned int cur_eop_idx){
    float V = (M(1,0)*(ground_point(cur_point_idx,1)-EOP(cur_eop_idx,2))) + (M(1,1)* (ground_point(cur_point_idx,2)-EOP(cur_eop_idx,3))) + (M(1,2)* (ground_point(cur_point_idx,3)-EOP(cur_eop_idx,4)));
    return V;
}

float wNumber(MatrixXd M, MatrixXd ground_point, MatrixXd EOP,unsigned int cur_point_idx, unsigned int cur_eop_idx){
    float W = (M(2,0)*(ground_point(cur_point_idx,1)-EOP(cur_eop_idx,2))) + (M(2,1)* (ground_point(cur_point_idx,2)-EOP(cur_eop_idx,3))) + (M(2,2)* (ground_point(cur_point_idx,3)-EOP(cur_eop_idx,4)));
    return W;
}

void colinearity_equations(MatrixXd IOP, float U, float V, float W, MatrixXd& colinear_points, unsigned int cur_iop_idx){
    colinear_points(0,0) = IOP(cur_iop_idx,1) - (IOP(cur_iop_idx,3)*(U/W));
    colinear_points(0,1) = IOP(cur_iop_idx,2) - (IOP(cur_iop_idx,3)*(V/W));
}


int main(){
    //Read all input data to eigen fixed matrices
    int check_rows = 0, check_cols = 0;
    getArraySize("check_points.txt", check_rows, check_cols);
    MatrixXd check_points(check_rows,check_cols);
    readData(check_points, "check_points.txt");

    int control_rows = 0, control_cols = 0;
    getArraySize("control_points.txt", control_rows, control_cols);
    MatrixXd control_points(control_rows,control_cols);
    readData(control_points, "control_points.txt");

    int EOP_rows = 0, EOP_cols = 0;
    getArraySize("EOP.txt", EOP_rows, EOP_cols);
    MatrixXd EOP(EOP_rows,EOP_cols);
    readData(EOP, "EOP.txt");

    int IOP_rows = 0, IOP_cols = 0;
    getArraySize("IOP.txt", IOP_rows, IOP_cols);
    MatrixXd IOP(IOP_rows,IOP_cols);
    readData(IOP, "IOP.txt");

    int img_rows = 0, img_cols = 0;
    getArraySize("image_point_observations.txt", img_rows, img_cols);
    MatrixXd img_points(img_rows,img_cols);
    readData(img_points, "image_point_observations.txt");

    int tie_rows = 0, tie_cols = 0;
    getArraySize("tie_points.txt", tie_rows, tie_cols);
    MatrixXd tie_points(tie_rows,tie_cols);
    readData(tie_points, "tie_points.txt");

    //Create weight matrix P size: number of x,y img point observations 20 img point observations each has 2 x,y and there are 2 images (20x2x2) = 80x80
    MatrixXd P(img_points.rows()*2,img_points.rows()*2);

    float sig_img = 0.000005; //stdev of img point observations
    for(int i=0;i<P.rows();i++){
        for(int j=0;j<P.cols();j++){
            if (i==j){
                P(i,j)=(1/(sig_img*sig_img));
            }
            else{
                P(i,j)=0;
            }
        }

    }

    for (int i=0;i<img_points.rows();i++){
        bool is_tie_point;
        bool is_control_point;
        unsigned int point_id = img_points(i,0);
        unsigned int image_id = img_points(i,1);

        unsigned int cur_eop_idx = findEOPIndex(EOP,image_id); //find the EOP index where the image ID of the EOP matches the image ID of the image point
        unsigned int cur_iop_idx = findIOPIndex(EOP(cur_eop_idx,1),IOP); //find the IOP index where the IOP camera ID matches the current indexed EOP's camera ID

        MatrixXd M(3,3);
        MatrixXd colinear_points(1,2);

        is_control_point = controlOrTiePoint(control_points,point_id); //find if the image point corresponds to a control point or tie point
        if (is_control_point == false){
            // Tie point loop. All observations are considered for both the Ao and Ae design matrices.
            is_tie_point = true;
            unsigned int cur_tie_idx = findtieIndex(point_id, tie_points);

            rotationMatrix(M,EOP,cur_eop_idx);
            float U = uNumber(M, tie_points, EOP, cur_tie_idx, cur_eop_idx);
            float V = vNumber(M, tie_points, EOP, cur_tie_idx, cur_eop_idx);
            float W = wNumber(M, tie_points, EOP, cur_tie_idx, cur_eop_idx);

            colinearity_equations(IOP,U,V,W,colinear_points,cur_iop_idx);
        }
        else{
            // Control point loop. GCP are considered constant so the observation is only included in the Ae matrix not the Ao.
            is_tie_point = false;
            unsigned int cur_control_idx = findControlIndex(point_id, control_points);

            rotationMatrix(M,EOP,cur_eop_idx);
            float U = uNumber(M, control_points, EOP, cur_control_idx, cur_eop_idx);
            float V = vNumber(M, control_points, EOP, cur_control_idx, cur_eop_idx);
            float W = wNumber(M, control_points, EOP, cur_control_idx, cur_eop_idx);

            colinearity_equations(IOP,U,V,W,colinear_points,cur_iop_idx);

        }

    }
}
