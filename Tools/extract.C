#include<iostream>
#include<fstream>

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
        child << xcom << "\t" << ycom << "\t" << radius << "\n";
      }
      child.close();
    }
    parent.close();
   
    std::ofstream bash("result.sh");
//    bash << "rm *.mp4" << std::endl;
    bash << "for i in {1.." << repeat << "};" << std::endl;
    bash << "do" << std::endl;
    bash << " gnuplot -e \" set xrange[0:2.98];\n set yrange[0:2];\n set size ratio 0.66667; \n set terminal jpeg size 1280,1024; \n plot '$i.txt' u 1:2:3 w circles lc \"black\" notitle, \'../grain_info.txt\' u 1:2:3 w circles lc \"black\" notitle \" > pic$i.jpeg;" << std::endl;
    bash << "done\n" << std::endl;
    bash << "ffmpeg -i pic%d.jpeg -pix_fmt yuv444p movie.mp4\n" << std::endl;
    //bash << "rm *.jpeg *.txt" << std::endl;
    
    bash.close(); 
    system("bash result.sh");
    system("rm result.sh");
}
