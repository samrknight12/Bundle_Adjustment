#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

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

int main(){

    int check_rows = 0, check_cols = 0;
    getArraySize("check_points.txt", check_rows, check_cols);
    int control_rows = 0, control_cols = 0;
    getArraySize("control_points.txt", control_rows, control_cols);
    int EOP_rows = 0, EOP_cols = 0;
    getArraySize("EOP.txt", EOP_rows, EOP_cols);
    int IOP_rows = 0, IOP_cols = 0;
    getArraySize("IOP.txt", IOP_rows, IOP_cols);
    int img_rows = 0, img_cols = 0;
    getArraySize("image_point_observations.txt", img_rows, img_cols);
    cout<<img_rows<<endl;

}
