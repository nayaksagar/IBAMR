#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>

int main(int argc, char *argv[]) {
    
    int num_of_structures = atoi(argv[1])  ;

    bool with_porous = true;
    

    std::cout<<"num_of_structs = "<< num_of_structures << "\t porous media = "<< with_porous << std::endl;
    /*****************************beginning***************************************/
    std::ifstream beginningFile("inputGenerator/beginning.txt");  // Open the source file for reading
    std::ofstream destinationFile("input2D");  // Open the destination file for writing

    if (beginningFile.is_open() && destinationFile.is_open()) {
        std::string line;
        while (std::getline(beginningFile, line)) {
            destinationFile << line << std::endl;  // Write the line to the destination file
        }

    } else {
        std::cout << "Failed to open the file(s)." << std::endl;
    }

    beginningFile.close();  // Close the source file
    /*****************************beginning*****************************************/




    /*********************************rho_solid***************************************/
    destinationFile << "rho_solid = ";
    for(int i=0; i<num_of_structures; i++){
            if(i==num_of_structures-1)
              destinationFile << "RHO \n";
            else
              destinationFile << "RHO,";
    }
    /*********************************rho_solid***************************************/
    
    
    /*********************************middle***************************************/
    std::ifstream middleFile("inputGenerator/middle.txt");
    if (middleFile.is_open() && destinationFile.is_open()) {
        std::string line;
        while (std::getline(middleFile, line)) {
            destinationFile << line << std::endl;  // Write the line to the destination file
        }

    } else {
        std::cout << "Failed to open the file(s)." << std::endl;
    }

    middleFile.close();  // Close the source file
    /*********************************middle***************************************/

    
    /*********************************ConstraintIB & IBStandardInitializer***************************************/
    destinationFile << "num_structures = " << num_of_structures << "\n";

   destinationFile << "ConstraintIBKinematics { \n\n";
   destinationFile << " structure_names = ";
   for(int i = 0;i<num_of_structures;i++){
     if(with_porous){
        if(i==num_of_structures-1)
        destinationFile << " \"porous\" \n";
        else
        destinationFile << " \"p_"<<i<<"\",";
     }
     else{
        if(i==num_of_structures-1)
        destinationFile << " \"p_"<<i<<"\"\n ";
        else
        destinationFile << " \"p_"<<i<<"\",";
     }
   }

   for(int i = 0;i<num_of_structures;i++){
     if(with_porous){
        if(i==num_of_structures-1){
            destinationFile << " porous {  \n	structure_names         \t = \"porous\" \n";
            destinationFile << "\tstructure_levels                 =  MAX_LEVELS - 1 \n";
            destinationFile << "\tcalculate_translational_momentum = 0,0,0 \n";
            destinationFile << "\tcalculate_rotational_momentum    = 0,0,0 \n";
            destinationFile << "\tlag_position_update_method       = \"CONSTRAINT_VELOCITY\" \n";
            destinationFile << "\ttagged_pt_identifier             = MAX_LEVELS - 1, 0 \n }\n\n}\n\n";
        }
        else{
            destinationFile << " p_"<<i<<" {  \n	structure_names         \t = \"p_"<<i<<"\" \n";
            destinationFile << "\tstructure_levels                 =  MAX_LEVELS - 1 \n";
            destinationFile << "\tcalculate_translational_momentum = 1,1,0 \n";
            destinationFile << "\tcalculate_rotational_momentum    = 0,0,0 \n";
            destinationFile << "\tlag_position_update_method       = \"CONSTRAINT_VELOCITY\"\n";
            destinationFile << "\ttagged_pt_identifier             = MAX_LEVELS - 1, 0 \n } \n\n";
        }
     }
     else{
        if(i==num_of_structures-1){
            destinationFile << " p_"<<i<<" {  \n	structure_names         \t = \"p_"<<i<<"\" \n";
            destinationFile << "\tstructure_levels                 =  MAX_LEVELS - 1 \n";
            destinationFile << "\tcalculate_translational_momentum = 1,1,0 \n";
            destinationFile << "\tcalculate_rotational_momentum    = 0,0,1 \n";
            destinationFile << "\tlag_position_update_method       = \"CONSTRAINT_VELOCITY\" \n";
            destinationFile << "\ttagged_pt_identifier             = MAX_LEVELS - 1, 0 \n }\n\n}\n\n";
        }
        else{
            destinationFile << " p_"<<i<<" {  \n	structure_names         \t = \"p_"<<i<<"\" \n";
            destinationFile << "\tstructure_levels                 =  MAX_LEVELS - 1 \n";
            destinationFile << "\tcalculate_translational_momentum = 1,1,0 \n";
            destinationFile << "\tcalculate_rotational_momentum    = 0,0,1 \n";
            destinationFile << "\tlag_position_update_method       = \"CONSTRAINT_VELOCITY\"\n";
            destinationFile << "\ttagged_pt_identifier             = MAX_LEVELS - 1, 0 \n } \n\n";
        }
     }

   }

/*
   destinationFile << "IBStandardInitializer { \n\n";
   destinationFile << " max_levels      = MAX_LEVELS \n";
   destinationFile << " structure_names = ";
   for(int i = 0;i<num_of_structures;i++){
     if(with_porous){
        if(i==num_of_structures-1)
        destinationFile << " \"porous\" \n";
        else
        destinationFile << " \"p_"<<i<<"\",";
     }
     else{
        if(i==num_of_structures-1)
        destinationFile << " \"p_"<<i<<"\"\n";
        else
        destinationFile << " \"p_"<<i<<"\",";
     }
   }
   
   for(int i = 0;i<num_of_structures;i++){
     if(with_porous){
        if(i==num_of_structures-1){
            destinationFile << "porous {  \n";
            destinationFile << "\t level_number            = MAX_LEVELS - 1 \n}\n\n}\n";
        }
        else{
            destinationFile << "p_"<<i<<" {  \n	level_number         \t =  MAX_LEVELS - 1 \n} \n";
        }
     }
     else{
        if(i==num_of_structures-1){
            destinationFile << "p_"<<i<<" {  \n";
            destinationFile << "\t level_number            = MAX_LEVELS - 1 \n}\n\n}\n";
        }
        else{
            destinationFile << "p_"<<i<<" {  \n";
            destinationFile << "\t level_number            =  MAX_LEVELS - 1 \n} \n";
        }
     }

   }*/
   
   destinationFile << "IBRedundantInitializer { \n\n";
   destinationFile << " max_levels      = MAX_LEVELS \n } \n";
    /*********************************ConstraintIB & IBStandardInitializer***************************************/



    /*********************************end***************************************/
    std::ifstream endFile("inputGenerator/end.txt");
    if (endFile.is_open() && destinationFile.is_open()) {
        std::string line;
        while (std::getline(endFile, line)) {
            destinationFile << line << std::endl;  // Write the line to the destination file
        }

        std::cout << "File copied successfully." << std::endl;
    } else {
        std::cout << "Failed to open the file(s)." << std::endl;
    }

    endFile.close();  // Close the source file 
    /*********************************end***************************************/

    destinationFile.close();  // Close the destination file
    
    
    /************************* Generate visit macro file ***********************/
   /* std::ofstream visitMacro;
    visitMacro.open("visit.txt");
    
    //Take the name of the pwd
    std::string path;
    std::filesystem::path currentPath = std::filesystem::current_path();
    path = currentPath.string();
    path += "/viz_FC";
    
  //  char path[100] = "/media/c2fd/c2fd/Work/Simulations/Porous/Code/2dPM/504particles/viz_FC";
    
	//visitMacro<< "def user_macro_load_levels_lag():\n";
    visitMacro<< "OpenDatabase(\"localhost:"<<path<<"/dumps.visit\", 0)\n";
    visitMacro<< "AddPlot(\"Pseudocolor\", \"U_magnitude\", 1, 1)\n\n";
    visitMacro<< "DrawPlots()\n\n";
    visitMacro<< "OpenDatabase(\"localhost:"<<path<<"/lag_data.visit\", 0)\n";
    for(int structure=0;structure<num_of_structures;structure++){
        // if(structure==num_of_structures-1){
        //     visitMacro << "AddPlot(\"Mesh\", \"porous_vertices\", 1, 1)\n";
        //     visitMacro << "MeshAtts = MeshAttributes()\n";
        //     visitMacro << "MeshAtts.legendFlag = 1 \n";
        //     visitMacro << "MeshAtts.lineWidth = 0 \n";
        //     visitMacro << "MeshAtts.meshColor = (192,192,192,255) \n";
        //     visitMacro << "MeshAtts.meshColorSource = MeshAtts.MeshCustom  # Foreground, MeshCustom, MeshRandom \n";
        //     visitMacro << "MeshAtts.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom, OpaqueRandom \n";
        //     visitMacro << "MeshAtts.opaqueMode = MeshAtts.Auto  # Auto, On, Off \n";
        //     visitMacro << "MeshAtts.pointSize = 0.05 \n";
        //     visitMacro << "MeshAtts.opaqueColor = (255, 255, 255, 255) \n";
        //     visitMacro << "MeshAtts.smoothingLevel = MeshAtts.NONE  # NONE, Fast, High \n";
        //     visitMacro << "MeshAtts.pointSizeVarEnabled = 0 \n";
        //     visitMacro << "MeshAtts.pointSizeVar = \"default\"\n";
        //     visitMacro << "MeshAtts.pointType = MeshAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere \n";
        //     visitMacro << "MeshAtts.showInternal = 0 \n";
        //     visitMacro << "MeshAtts.pointSizePixels = 2 \n";
        //     visitMacro << "MeshAtts.opacity = 1 \n";
        //     visitMacro << "SetPlotOptions(MeshAtts) \n\n\n";
        // }
        // else{
            visitMacro << "AddPlot(\"Mesh\", \"p_"<<structure<<"_vertices\", 1, 1)\n";
            visitMacro << "MeshAtts = MeshAttributes()\n";
            visitMacro << "MeshAtts.legendFlag = 1 \n";
            visitMacro << "MeshAtts.lineWidth = 0 \n";
            visitMacro << "MeshAtts.meshColor = (0, 0, 0, 0) \n";
            visitMacro << "MeshAtts.meshColorSource = MeshAtts.MeshCustom  # Foreground, MeshCustom, MeshRandom \n"; 
            visitMacro << "MeshAtts.opaqueColorSource = MeshAtts.Background  # Background, OpaqueCustom, OpaqueRandom \n"; 
            visitMacro << "MeshAtts.opaqueMode = MeshAtts.Auto  # Auto, On, Off \n";
            visitMacro << "MeshAtts.pointSize = 0.05 \n";
            visitMacro << "MeshAtts.opaqueColor = (255, 255, 255, 255) \n";
            visitMacro << "MeshAtts.smoothingLevel = MeshAtts.NONE  # NONE, Fast, High \n";
            visitMacro << "MeshAtts.pointSizeVarEnabled = 0 \n";
            visitMacro << "MeshAtts.pointSizeVar = \"default\"\n";
            visitMacro << "MeshAtts.pointType = MeshAtts.Point  # Box, Axis, Icosahedron, Octahedron, Tetrahedron, SphereGeometry, Point, Sphere \n";
            visitMacro << "MeshAtts.showInternal = 0 \n";
            visitMacro << "MeshAtts.pointSizePixels = 2 \n";
            visitMacro << "MeshAtts.opacity = 1 \n";
            visitMacro << "SetPlotOptions(MeshAtts) \n\n\n";
            
            if(structure > 0 && structure % 10 == 0 && structure != num_of_structures-1)
              visitMacro << "DrawPlots()\n";
        // }
    }

    visitMacro<< "DrawPlots()\n";

    visitMacro.close();
*/
    return 0;
}
