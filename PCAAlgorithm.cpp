#include <iostream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

// Xhanti Singatha
// SNGXHA002

using namespace std;
using namespace Eigen;

vector<double> calculate_covariance(vector<double> x1, vector<double> y1);

vector<double> calculate_covariance(vector<double> x1, vector<double> y1){
                       vector<double> output;
                       
                       double yrow1, xcol1;
                       double yrow2, xcol2;
                                             
                       for (int m = 0; m < x1.size(); m++){
                           yrow1 += (x1[m]*x1[m]);      
                       }
                       
                       output.push_back(yrow1);
                       
                       for (int n = 0; n < x1.size(); n++){
                           xcol1 += (x1[n]*y1[n]);      
                       }
                       
                       output.push_back(xcol1);
                       
                       for (int o = 0; o < x1.size(); o++){
                           yrow2 += (y1[o]*x1[o]);      
                       }
                       
                       output.push_back(yrow2);
                       
                       for (int p = 0; p < x1.size(); p++){
                           xcol2 += (y1[p]*y1[p]);      
                       }
                       
                       output.push_back(xcol2);
                       
                       return output;
                       

}


int main(){

    //vector<double> x = {2.5, 0.5, 2.2, 1.9, 3.1, 2.3, 2.0, 1.0, 1.5, 1.1};
    //vector<double> y = {2.4, 0.7, 2.9, 2.2, 3.0, 2.7, 1.6, 1.1, 1.6, 0.9};
    
    ifstream myfile;
    myfile.open("2018-AvgRainfall_mm_.txt");
    string line = "";
    string line1 = "";
        
    vector<string> names;
    vector<double> x;
    vector<double> y;
    
    string name;
    string xstr;
    string ystr;

    int counter = 0;
    
    vector<string> temp;
    stringstream ss;
    stringstream ss1;
    if (myfile.is_open()){
        //first line
        getline(myfile, line);
        
        // second line
        getline(myfile, line);
        ss << line;
        string var, num;
        ss >> name >> xstr >> ystr >> var >> num;
        x.push_back(atof(xstr.c_str()));
        y.push_back(atof(ystr.c_str()));
        
        while (!myfile.eof()){
               getline(myfile, line1);
               ss1 << line1;
               ss1 >> name >> xstr >> ystr;
               x.push_back(atof(xstr.c_str()));
               y.push_back(atof(ystr.c_str()));
        }
        myfile.close();    
    }
    
    else {
        cout << "Unable to open file";
    }
    
    cout << "\n";
    
    int sample_size = x.size();
    
    double xmean = 0.0;
    double xtotal = 0.0;
    for (int i = 0; i < x.size(); i++ ){
         xtotal += x[i];
    }
    
    xmean = xtotal/x.size();
    
    double ymean = 0.0; 
    double ytotal = 0.0;
    for (int j = 0; j < y.size(); j++){
        ytotal += y[j];
    }
    
    ymean = ytotal/y.size();
    
    vector<double> xReCenter;
    for (int k = 0; k < x.size(); k++){
        xReCenter.push_back(x[k] - xmean);
    }
    
    vector<double> yReCenter; 
    for (int l = 0; l < y.size(); l++){
        yReCenter.push_back(y[l] - ymean);
    }
    
    vector<double> copy = calculate_covariance(xReCenter, yReCenter);
    
    cout << "Covariance Matrix:" << "\n";
    for (int q = 0; q < copy.size(); q++){
        if (q == 2){
           cout << "\n";
        }
        cout << copy[q]*(1.0/(sample_size-1)) << " ";
    }
    
    cout << "\n";
    cout << "\n";
    
    MatrixXd a(2,2);
    a(0,0) = copy[0];
    a(0,1) = copy[1];
    a(1,0) = copy[2];
    a(1,1) = copy[3];
    
    EigenSolver<MatrixXd> eigensolver(a);
    
    cout << "Eigenvalues for the covariance matrix:" << "\n"; 
    cout << eigensolver.eigenvalues();
    
    cout << "\n";
    cout << "\n";
    
    cout << "Eigenvectors for the covariance matrix as columns:" << "\n";
    cout  << eigensolver.eigenvectors();
    
}