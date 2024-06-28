#include<iostream>
#include<fstream>
#include<vector>
#include<sstream>
#include<filesystem>
#include<numeric>
#include<algorithm>
#include<cmath>

// int findClosestIndex(const std::vector<double>& data, double value) {
//     auto it = std::lower_bound(data.begin(), data.end(), value);
//     if (it == data.end()) {
//         return data.size() - 1;
//     }
//     if (it == data.begin()) {
//         return 0;
//     }
//     if (std::fabs(*it - value) < std::fabs(*(it - 1) - value)) {
//         return std::distance(data.begin(), it);
//     } else {
//         return std::distance(data.begin(), it - 1);
//     }
// }

// struct row{
//     double firstCol;
//     int secondCol;
// }

int main(int argc, char *argv[]){
   int no_of_particles = atoi(argv[1]);
   double time, xcom, ycom;
   double radius = 0.1;
   int particle_id;
   double empty;
   double constrictionStartsAt = 4.0;
   double allowanceToDetectClogging = 4.0 * radius;
   double Vf = 1.0; //in cm/s
   
   std::vector<std::vector<int>> allData;
   std::vector<double> timeData;
   for (int fileCount = 1; fileCount <= 6; fileCount++)
   {
        std::string filename = "SedimentingCylinder_COM_coordinates_" + std::to_string(fileCount);

        std::ifstream file(filename);
        int lineCount = 0;
        std::string line;
        while(std::getline(file, line)){
                lineCount ++;
        }
        file.close();
        
        int producedParticleCount;
        int repeat= (int)(lineCount/(no_of_particles+3));
        int recordEvery = 10;  
        std::ifstream parent(filename);
        std::vector<int> producedParticles;
        for(int i=1; i<=repeat; i++){
            parent >> time;
            if (fileCount == 1)
                timeData.push_back((time * Vf)/(2.0 * radius));
            producedParticleCount = 0;
            for(int j = 1; j <= no_of_particles; j++)
            {
                parent >> particle_id;
                parent >> xcom;
                parent >> ycom;
                if (xcom > constrictionStartsAt + allowanceToDetectClogging)
                    producedParticleCount ++;
            }
                producedParticles.push_back(producedParticleCount);
        }
        allData.push_back(producedParticles);
        parent.close();
   }
   std::cout << "hello" << std::endl;

   std::cout << allData[0].size();
   int sum[allData[0].size()] = {0}, maxm[allData[0].size()] = {0}, minm[allData[0].size()] = {0};
   double mean[allData[0].size()] = {0.0};
   for(int fileCount = 0; fileCount < 6; fileCount++)
   {
        int i = 0;
        for(const auto& particleInTrial : allData[fileCount])
        {
            mean[i] += particleInTrial/6.0;
            if (fileCount == 0)
            {
                maxm[i] = allData[fileCount][i];
                minm[i] = allData[fileCount][i];
            }
            else
            {
                maxm[i] = fmax(particleInTrial, allData[fileCount-1][i]);
                minm[i] = fmin(particleInTrial, allData[fileCount-1][i]);
            }
            i++;           
        }
   }
    std::ofstream plotFile("particlesVsTime.txt");
    for(int i = 0; i < allData[0].size(); ++i)
        plotFile << timeData[i] << "\t" << mean[i] << "\t" << maxm[i] << "\t" << minm[i] << std::endl;
    plotFile.close();
}