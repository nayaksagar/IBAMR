#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<filesystem>

int main(int argc, char *argv[]){
   int no_of_particles = atoi(argv[1]);
   double time, xcom, ycom;
   double radius = 0.01;
   int particle_id;
   double empty;

   std::ifstream file("SedimentingCylinder_COM_coordinates");
   int lineCount = 0;
   std::string line;
   while(std::getline(file, line)){
	      lineCount ++;
   }
   file.close();
  
   int repeat= (int)(lineCount/(no_of_particles+3));  
   std::ifstream parent("SedimentingCylinder_COM_coordinates");
    for(int i=1; i<=repeat; i++){
      parent >> time;
      std::cout << time << std::endl; 
      std::ofstream child(std::to_string(i) + ".txt");
      for(int j=1;j<=no_of_particles;j++){
        parent >> particle_id;
        parent >> xcom;
        parent >> ycom;
        child << time << "\n";
        child << xcom << "\t" << ycom << "\t" << radius << "\n";
      }
      child.close();
    }
    parent.close();
    

    //write vtk files
    double x, y, r;

    if(!(std::filesystem::exists("vtk") && std::filesystem::is_directory("vtk")))
      system("mkdir vtk");
      
    for(int i = 1; i <= repeat; i++){ 
      std::ifstream readtxt(std::to_string(i) + ".txt");
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
      std::ofstream vtkFile("vtk/coordinates" + std::to_string(i) + ".vtk");
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
    }

    // writing shell script 
    std::ofstream bash("result.sh");
    bash << "rm *.mp4" << std::endl;
    bash << "for i in {1.." << repeat << "};" << std::endl;
    bash << "do" << std::endl;
    bash << " gnuplot -e \" title = system('head -n 1 $i.txt'); set title title; set xrange[0:2.98];\n set yrange[0:2];\n set size ratio 0.66667; \n set terminal jpeg size 1280,1024 font \',20\'; \n plot '$i.txt' u 1:2:3 w circles lc \"black\" notitle, \'../grain_info.txt\' u 1:2:3 w circles lc \"black\" notitle \" > pic$i.jpeg;" << std::endl;
    bash << "done\n" << std::endl;
    bash << "ffmpeg -i pic%d.jpeg -pix_fmt yuv444p movie.mp4\n" << std::endl;
    bash << "rm *.jpeg " << std::endl;
    bash << "mv *.txt vtk " << std::endl;
    
    bash.close(); 
    system("bash result.sh");
   // system("rm result.sh");
}
