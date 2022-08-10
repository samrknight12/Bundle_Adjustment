#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <Eigen/Dense>

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


void readData(MatrixXf& data, string filename)
{
    ifstream file;
    file.open(filename);
    for (int i = 0; i<data.rows(); i++){
        for (int j = 0; j<data.cols(); j++){
            file >> data(i,j);
        }
    }
}


int main(){

    int check_rows = 0, check_cols = 0;
    getArraySize("check_points.txt", check_rows, check_cols);
    MatrixXf check_points(check_rows,check_cols);
    readData(check_points, "check_points.txt");

    int control_rows = 0, control_cols = 0;
    getArraySize("control_points.txt", control_rows, control_cols);
    MatrixXf control_points(control_rows,control_cols);
    readData(control_points, "control_points.txt");

    int EOP_rows = 0, EOP_cols = 0;
    getArraySize("EOP.txt", EOP_rows, EOP_cols);
    MatrixXf EOP(EOP_rows,EOP_cols);
    readData(EOP, "EOP.txt");

    int IOP_rows = 0, IOP_cols = 0;
    getArraySize("IOP.txt", IOP_rows, IOP_cols);
    MatrixXf IOP(IOP_rows,IOP_cols);
    readData(IOP, "IOP.txt");

    int img_rows = 0, img_cols = 0;
    getArraySize("image_point_observations.txt", img_rows, img_cols);
    MatrixXf img_points(img_rows,img_cols);
    readData(img_points, "image_point_observations.txt");



}
