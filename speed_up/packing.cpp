#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>

struct Circle {
    double x;
    double y;
    double radius;
};

// Function to check if two circles overlap
bool circlesOverlap(const Circle& c1, const Circle& c2) {
    double dx = c1.x - c2.x;
    double dy = c1.y - c2.y;
    double distance = std::sqrt(dx * dx + dy * dy);
    return distance < (c1.radius + c2.radius);
}

// Function to check if a circle is within the boundaries of the domain
bool circleWithinBounds(const Circle& circle, double domainSize) {
    return (circle.x - circle.radius >= 0 && circle.x + circle.radius <= domainSize &&
            circle.y - circle.radius >= 0 && circle.y + circle.radius <= domainSize);
}

// Function to generate a random circle within the domain
Circle generateRandomCircle(double domainSize, double minRadius, double maxRadius) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<double> disRadius(minRadius, maxRadius);
    std::uniform_real_distribution<double> disPos(0, domainSize);

    Circle circle;
    circle.radius = disRadius(gen);
    circle.x = disPos(gen);
    circle.y = disPos(gen);

    return circle;
}

// Function to generate a porous media using circle packing
std::vector<Circle> generatePorousMedia(double domainSize, double targetPorosity, double minRadius, double maxRadius) {
    std::vector<Circle> circles;
    double totalArea = domainSize * domainSize;
    double targetArea = totalArea * (1.0 - targetPorosity);

    while (true) {
        Circle circle = generateRandomCircle(domainSize, minRadius, maxRadius);

        bool isValid = true;
        for (const Circle& existingCircle : circles) {
            if (circlesOverlap(existingCircle, circle) || !circleWithinBounds(circle, domainSize)) {
                isValid = false;
                break;
            }
        }

        if (isValid) {
            circles.push_back(circle); 

            // Check if target area is reached
            double currentArea = 0;
            for (const Circle& existingCircle : circles) {
                currentArea += M_PI * existingCircle.radius * existingCircle.radius;
            }
            if (currentArea >= targetArea) {
                break;
            }
        }
    }

    return circles;
}

int main() {
    double domainSize = 2.0; // Size of the square domain
    double targetPorosity = 0.6; // Target porosity (between 0 and 1)
    double minRadius = 0.05; // Minimum circle radius
    double maxRadius = 0.1; // Maximum circle radius

    std::vector<Circle> porousMedia = generatePorousMedia(domainSize, targetPorosity, minRadius, maxRadius);

    // Print the generated circles
    // for (const Circle& circle : porousMedia) {
        // std::cout << "Circle: (x=" << circle.x << ", y=" << circle.y << ", radius=" << circle.radius << ")\n";
        // std::cout << "[" << circle.x << ", " << circle.y << ", " << circle.radius << "],\n";
    // }

    //generate lagrangian points
    double Lx = 3.0;
    double Ly = 2.0;
    int Nx = 48*4;
    int Ny = 32*4;
    double dx = Lx/Nx;
    double dy = Ly/Ny;
    

    double radius, circumference, dtheta, theta, X, Y;
    int num_surf_pts;
    double pm_pos_x = 1.0, pm_pos_y=0.0; //coords of the left bottom corner of the porous media
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
             porous_file << pm_pos_x+X << "\t" << pm_pos_y+Y << "\n";
           }
           radius = radius - dx;
        }
    }
    porous_file.close();
    
    std::ofstream structuredata;
    structuredata.open("grain_info.txt");
    structuredata << circle_num << "\n";
    for(const Circle& circle:porousMedia){
        structuredata << pm_pos_x+circle.x << "\t" << pm_pos_y+circle.y << "\t" << circle.radius << "\n";
    }
    structuredata.close();
    
     return 0;
}
