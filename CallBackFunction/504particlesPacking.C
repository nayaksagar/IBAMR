//Author: Sagar G. Nayak
//Affiliation:UQIDAR, IITDelhi, Univerisity of Queensland


#include<iostream>
#include<math.h>
#include<fstream>

using namespace std;

int main(){
        const double DIAMETER=0.0625; 
        double RADIUS=0.5*DIAMETER, radius;
        double circumference, dtheta, theta, X, Y;

        double Lx = 2.0, Ly = 2.0;
        const int Nx = 64*4, Ny = 64*4;
        double dx = Lx/Nx, dy = Ly/Ny;
        int num_surf_pts;
        int num_in_x=28;
        int num_in_y=18;
        double x_com_far[num_in_x], x_com_ner[num_in_x], y_com[num_in_y];
        x_com_far[0] = 0.0 + 6.0*DIAMETER/16.0 + RADIUS;
        x_com_ner[0] = 0.0 + 4.0*DIAMETER/16.0 + RADIUS;
        y_com[0] =  Ly - 6.0*DIAMETER/16.0 - RADIUS;
        
        for(int i=1;i<num_in_x;i++){
          x_com_far[i] = x_com_far[i-1] + DIAMETER + 2.0*DIAMETER/16.0;
          x_com_ner[i] = x_com_ner[i-1] + DIAMETER + 2.0*DIAMETER/16.0;
        }

        for(int i=1;i<num_in_y;i++)
          y_com[i] = y_com[i-1] - DIAMETER - 2.0*DIAMETER/16.0;


        int circle_num=0;
        std::ofstream comFile("comFile.dat");
        for(int j=0;j<num_in_y;j++){
          for(int i=0;i<num_in_x;i++){
            radius = RADIUS;
            if(j%2==0)
             comFile << x_com_far[i] << "\t";
            else
             comFile << x_com_ner[i] << "\t";
             
             comFile << y_com[j] << "\n";
            
            // char filename[50];
            // sprintf(filename,"p_%d.vertex",circle_num);
            // ofstream particle_file;
            // particle_file.open(filename);
            
            // int total_num_pts=0;
            // while(radius>0.0/*0.75*RADIUS*/){
            //     circumference = 2.0*M_PI*radius;
            //     num_surf_pts = round(circumference/dx);
            //     total_num_pts = total_num_pts + num_surf_pts;
            //     radius = radius - dx;
            // }

            // particle_file<<total_num_pts<<endl;
            // radius = RADIUS;
            // while(radius>0.0/*0.75*RADIUS*/){
            // circumference = 2.0*M_PI*radius;
            // num_surf_pts = round(circumference/dx);
            // dtheta = (2.0*M_PI)/num_surf_pts;
            // for(int k=1; k<=num_surf_pts; k++){
            //     theta = (k-1)*dtheta;
            //     if(j%2==0)
            //     X = x_com_far[i] + radius*cos(theta);
            //     else
            //     X = x_com_ner[i] + radius*cos(theta);
                
            //     Y = y_com[j] + radius*sin(theta);
            //     particle_file << X << "\t" << Y << "\n";
            // }
            // radius = radius - dx;
            // }
            // circle_num ++;
            // particle_file.close();
          }
        }
    comFile.close();
}
