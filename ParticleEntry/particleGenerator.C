//Author: Sagar G. Nayak
//Affiliation:UQIDAR, IITDelhi, Univerisity of Queensland


#include<iostream>
#include<math.h>
#include<fstream>

using namespace std;

int main(){
  const double DIAMETER=0.02; 
  double RADIUS=0.5*DIAMETER, radius;
  double circumference, dtheta, theta, X, Y;

  double Lx = 3.0, Ly = 2.0;
  int num_surf_pts;
  int num_in_x=40;
  int num_in_y=60;
  double y_com_far[num_in_y], y_com_ner[num_in_y], x_com[num_in_x], y_com[num_in_y];
  x_com[0] =  0.02;//    - (num_in_x*0.5*3.0*RADIUS) + 1.5*RADIUS;     // 0.615 ;//+ 3.0*0.5*DIAMETER;//+ 6.0*DIAMETER/16.0 + RADIUS;
  y_com[0] =  0.5*Ly - (num_in_y*0.5*3.0*RADIUS) + 1.5*RADIUS;//0.715 ;//+ 3.0*0.5*DIAMETER; 
  for(int i=1;i<num_in_y;i++){
      y_com[i] = y_com[i-1] + 3.0*0.5*DIAMETER; 
  }

  for(int i=1;i<num_in_x;i++)
    x_com[i] = x_com[i-1];// + 3.0*0.5*DIAMETER ;//+ 2.0*DIAMETER/16.0;


  int circle_num=0;
  std::ofstream comFile("comFile.dat");
  for(int j=0;j<num_in_x;j++){
    for(int i=0;i<num_in_y;i++){
      radius = RADIUS;
      comFile << x_com[j] << "\t";
      comFile << y_com[i] << "\n";
    }
  }
    comFile.close();
}
