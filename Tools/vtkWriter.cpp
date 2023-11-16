#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
int main() {
    double x, y, r;
   
    for(int i = 1; i <= 200; i++){ 
    std::ifstream readtxt("vtk/" + std::to_string(i) + ".txt");
    std::string line;

    std::vector<std::vector<double>> coordinates;
    
    while(std::getline(readtxt,line)){
	std::istringstream iss(line);
        
 	if(!(iss >> x >> y >> r))
	    continue;
	
        std::vector<double> coord = {x,y,1.0};
	coordinates.push_back(coord);
     }	


    // Output VTK file
    std::ofstream vtkFile("coordinates" + std::to_string(i) + ".vtk");
    if (!vtkFile.is_open()) {
        std::cerr << "Error opening file.\n";
        return 1;
    }

    // Write VTK header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "vtk output\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET POLYDATA\n";

    // Write points
    vtkFile << "POINTS " << coordinates.size() << " float\n";
    for (const auto& coord : coordinates) {
        vtkFile << coord[0] << " " << coord[1] << " " << 0.0 << "\n";
    }
    vtkFile << "POINT_DATA " << coordinates.size() << std::endl;
    vtkFile << "\n" << "SCALARS rho float" << std::endl;
    vtkFile << "LOOKUP_TABLE default" << std::endl;
    for (const auto& coord : coordinates) {
        vtkFile << coord[2] << std::endl;
    }

    //// Write vertices
    //vtkFile << "VERTICES 1 " << coordinates.size() << "\n";
    //for (size_t i = 0; i < coordinates.size(); ++i) {
    //    vtkFile << i << " ";
    //}

    }
    return 0;
}

