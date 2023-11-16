#include<iostream>
#include<math.h>
#include<fstream>

using namespace std;

int main(){
        double RADIUS=0.05, radius;
        double circumference, dtheta, theta, X, Y;

        double Lx = 1.0, Ly = 1.0;
        const int Nx = 32*4, Ny = 32*4;
        double dx = Lx/Nx, dy = Ly/Ny;
        int num_surf_pts;
        int num_in_x=10;
        int num_in_y=5;
        double x_com[num_in_x], y_com[num_in_y];
        x_com[0] = 0.55;
        y_com[0] = 1.95;
        for(int i=1;i<num_in_x;i++)
          x_com[i] = x_com[i-1]-2.0*RADIUS;

        for(int i=1;i<num_in_y;i++)
          y_com[i] = y_com[i-1]-2.0*RADIUS;

        int circle_num=0;
        
        
        for(int i=0;i<num_in_x;i++){
            for(int j=0;j<num_in_y;j++){
                radius = RADIUS;
                char filename[50];
                sprintf(filename,"particle_%d.vertex",circle_num);
                ofstream particle_file;
                particle_file.open(filename);
                
                int total_num_pts=0;
		while(radius>0.0/*0.75*circle.radius*/){
			circumference = 2.0*M_PI*radius;
			num_surf_pts = round(circumference/dx);
			total_num_pts = total_num_pts + num_surf_pts;
		        radius = radius - dx;
		}
		
                particle_file<<total_num_pts<<endl;
                radius = RADIUS;
                while(radius>0.0/*0.75*circle.radius*/){
                circumference = 2.0*M_PI*radius;
                num_surf_pts = round(circumference/dx);
                dtheta = (2.0*M_PI)/num_surf_pts;
                for(int k=1; k<=num_surf_pts; k++){
                    theta = (k-1)*dtheta;
                    X = x_com[i] + radius*cos(theta);
                    Y = y_com[j] + radius*sin(theta);
                    particle_file << X << "\t" << Y << "\n";
                }
                radius = radius - dx;
                }
                circle_num ++;
                particle_file.close();
            }
        }
}
