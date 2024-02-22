//Author: Sagar G. Nayak
//Affiliation:  IIT Delhi, University of Queensland

/* This code takes the grain data from "grain_info.txt" and reconstructs the porous media*/

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <sstream>

struct Circle {
    double x;
    double y;
    double radius;
};

int main() {
    double domainSize = 2.0; // Size of the square domain
    double targetPorosity = 0.6; // Target porosity (between 0 and 1)
    double minRadius = 0.05; // Minimum circle radius
    double maxRadius = 0.1; // Maximum circle radius

    std::vector<Circle> porousMedia;// = generatePorousMedia(domainSize, targetPorosity, minRadius, maxRadius);
    std::ifstream readGrainFile("grain_info.txt");
    std::string line;
    Circle grain;
    while(std::getline(readGrainFile,line))
    {
            std::istringstream iss(line);
            if (!(iss >> grain.x >> grain.y >> grain.radius)) {
                // Handle the case where the line does not contain three valid numbers
                std::cerr << "Error reading line: " << line << std::endl;
                continue; // Skip to the next iteration
            }
            porousMedia.push_back(grain);
    }

    //generate lagrangian points
    double Lx = 3.0;
    double Ly = 2.0;
    int Nx = 96*4;
    int Ny = 64*4; 
    double dx = Lx/Nx;
    double dy = Ly/Ny;
    

    double radius, circumference, dtheta, theta, X, Y;
    int num_surf_pts;
    int circle_num=0;

    std::ofstream porous_file;
    porous_file.open("porous.vertex");
    for(const Circle& circle: porousMedia){
        radius = circle.radius;
        circle_num++;
        while(radius>0.0/*0.75*circle.radius*/){
           circumference = 2.0*M_PI*radius;
           num_surf_pts = round(circumference/dx);
           dtheta = (2.0*M_PI)/num_surf_pts;
           for(int i=1; i<=num_surf_pts; i++){
             theta = (i-1)*dtheta;
             X = circle.x + radius*cos(theta);
             Y = circle.y + radius*sin(theta);
             porous_file << X << "\t" << Y << "\n";
           }
           radius = radius - dx;
        }
    }
    porous_file.close();
    
     return 0;
}
