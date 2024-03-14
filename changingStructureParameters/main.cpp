// Filename: main.cpp
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE

// Config files
#include <IBAMR_config.h> 
#include <IBTK_config.h>
#include <SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>
#include <LoadBalancer.h>
#include <StandardTagAndInitialize.h>

// Headers for application-specific algorithm/data structure objects
#include <ConstraintIBMethod.h>
#include <ibamr/IBExplicitHierarchyIntegrator.h>
#include <ibamr/IBStandardForceGen.h>
#include <ibamr/IBRedundantInitializer.h>
#include <ibamr/INSCollocatedHierarchyIntegrator.h>
#include <ibamr/INSStaggeredHierarchyIntegrator.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/LData.h>
#include <ibtk/muParserCartGridFunction.h>
#include <ibtk/muParserRobinBcCoefs.h>

// Application
#include "CartGridBodyForce.h"
#include "ForceProjector.h"
#include "RigidBodyKinematics.h"

#include <limits>
#include <chrono>

int max_finest_ln;
std::vector<double> xcom;
std::vector<double> ycom;
int num_structures;

// Function prototypes
void output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                 Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
                 LDataManager* l_data_manager,
                 const int iteration_num,
                 const double loop_time,
                 const string& data_dump_dirname);

void
generate_structure(const unsigned int& strct_num,
                     const int& ln,
                     int& num_vertices,
                     std::vector<IBTK::Point>& vertex_posn,
                     void* /*ctx*/)
{
    if ( strct_num == num_structures-1){
        int vertex_id = 0; 
       //read porous data from file
       std::ifstream porousData;
       porousData.open("porous.vertex");
       porousData >> num_vertices;
       vertex_posn.resize(num_vertices);
       for (int i = 1; i<=num_vertices; i++){
        IBTK::Point& X = vertex_posn[vertex_id];
        porousData >> X(0);
        porousData >> X(1);
            vertex_id++;
       }
       porousData.close();
    }
    else{
        num_vertices = 1;
        vertex_posn.resize(num_vertices);
        IBTK::Point& X = vertex_posn[0];
        X(0) = xcom[strct_num];
        X(1) = ycom[strct_num];
    } 
  
    return;
} // generate_structure

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int
main(int argc, char* argv[])
{
    // Initialize PETSc, MPI, and SAMRAI.
    PetscInitialize(&argc, &argv, NULL, NULL);
    SAMRAI_MPI::setCommunicator(PETSC_COMM_WORLD);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();
    
    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        Pointer<AppInitializer> app_initializer = new AppInitializer(argc, argv, "IB.log");
        Pointer<Database> input_db = app_initializer->getInputDatabase();

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && !app_initializer->getVisItDataWriter().isNull();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_postproc_data = app_initializer->dumpPostProcessingData();
        const int postproc_data_dump_interval = app_initializer->getPostProcessingDataDumpInterval();
        const string postproc_data_dump_dirname = app_initializer->getPostProcessingDataDumpDirectory();
        if (dump_postproc_data && (postproc_data_dump_interval > 0) && !postproc_data_dump_dirname.empty())
        {
            Utilities::recursiveMkdir(postproc_data_dump_dirname);
        }

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        Pointer<INSHierarchyIntegrator> navier_stokes_integrator = new INSStaggeredHierarchyIntegrator(
            "INSStaggeredHierarchyIntegrator",
            app_initializer->getComponentDatabase("INSStaggeredHierarchyIntegrator"));

        num_structures = input_db->getIntegerWithDefault("num_structures", 1);
        Pointer<ConstraintIBMethod> ib_method_ops = new ConstraintIBMethod(
            "ConstraintIBMethod", app_initializer->getComponentDatabase("ConstraintIBMethod"), num_structures);
        Pointer<IBHierarchyIntegrator> time_integrator =
            new IBExplicitHierarchyIntegrator("IBHierarchyIntegrator",
                                              app_initializer->getComponentDatabase("IBHierarchyIntegrator"),
                                              ib_method_ops,
                                              navier_stokes_integrator);

        Pointer<CartesianGridGeometry<NDIM> > grid_geometry = new CartesianGridGeometry<NDIM>(
            "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        Pointer<PatchHierarchy<NDIM> > patch_hierarchy = new PatchHierarchy<NDIM>("PatchHierarchy", grid_geometry);

        Pointer<StandardTagAndInitialize<NDIM> > error_detector =
            new StandardTagAndInitialize<NDIM>("StandardTagAndInitialize",
                                               time_integrator,
                                               app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        Pointer<BergerRigoutsos<NDIM> > box_generator = new BergerRigoutsos<NDIM>();
        Pointer<LoadBalancer<NDIM> > load_balancer =
            new LoadBalancer<NDIM>("LoadBalancer", app_initializer->getComponentDatabase("LoadBalancer"));
        Pointer<GriddingAlgorithm<NDIM> > gridding_algorithm =
            new GriddingAlgorithm<NDIM>("GriddingAlgorithm",
                                        app_initializer->getComponentDatabase("GriddingAlgorithm"),
                                        error_detector,
                                        box_generator,
                                        load_balancer);

        // Configure the IB solver.
        Pointer<IBRedundantInitializer> ib_initializer = new IBRedundantInitializer(
            "IBRedundantInitializer", app_initializer->getComponentDatabase("IBRedundantInitializer"));
        tbox::Array<std::string> structure_names =
            app_initializer->getComponentDatabase("ConstraintIBKinematics")->getStringArray("structure_names"); //has nothing to do with constraintIBKinematics, takes only structure names
        //read COM of all structures at once
        double XCOM, YCOM;
        std::vector<std::string> struct_list_vec;
        std::ifstream strctFile("comFile.dat");
        for(int i = 0; i < num_structures - 1; ++i){ // -1 because last structure is porous media
            struct_list_vec.push_back(structure_names[i]);
            strctFile >> XCOM >> YCOM;
            xcom.push_back(XCOM);
            ycom.push_back(YCOM);
        }
        strctFile.close();
	    struct_list_vec.push_back(structure_names[num_structures-1]); //uncomment this only when there is porous media

        max_finest_ln = input_db->getInteger("MAX_LEVELS") - 1; // change this according to structure if you want
        ib_initializer->setStructureNamesOnLevel(max_finest_ln, struct_list_vec); //takes structure names
	    ib_initializer->registerInitStructureFunction(generate_structure);
        ib_method_ops->registerLInitStrategy(ib_initializer);
        Pointer<IBStandardForceGen> ib_force_fcn = new IBStandardForceGen();
        ib_method_ops->registerIBLagrangianForceFunction(ib_force_fcn);

        // Create Eulerian initial condition specification objects.
        if (input_db->keyExists("VelocityInitialConditions"))
        {
            Pointer<CartGridFunction> u_init = new muParserCartGridFunction(
                "u_init", app_initializer->getComponentDatabase("VelocityInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerVelocityInitialConditions(u_init);
        }

        if (input_db->keyExists("PressureInitialConditions"))
        {
            Pointer<CartGridFunction> p_init = new muParserCartGridFunction(
                "p_init", app_initializer->getComponentDatabase("PressureInitialConditions"), grid_geometry);
            navier_stokes_integrator->registerPressureInitialConditions(p_init);
        }

        // Create Eulerian boundary condition specification objects (when necessary).
        const IntVector<NDIM>& periodic_shift = grid_geometry->getPeriodicShift();
        vector<RobinBcCoefStrategy<NDIM>*> u_bc_coefs(NDIM);
        if (periodic_shift.min() > 0)
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                u_bc_coefs[d] = NULL;
            }
        }
        else
        {
            for (unsigned int d = 0; d < NDIM; ++d)
            {
                ostringstream bc_coefs_name_stream;
                bc_coefs_name_stream << "u_bc_coefs_" << d;
                const string bc_coefs_name = bc_coefs_name_stream.str();

                ostringstream bc_coefs_db_name_stream;
                bc_coefs_db_name_stream << "VelocityBcCoefs_" << d;
                const string bc_coefs_db_name = bc_coefs_db_name_stream.str();

                u_bc_coefs[d] = new muParserRobinBcCoefs(
                    bc_coefs_name, app_initializer->getComponentDatabase(bc_coefs_db_name), grid_geometry);
            }
            navier_stokes_integrator->registerPhysicalBoundaryConditions(u_bc_coefs);
        }

        // Create Eulerian body force function specification objects.
        if (input_db->keyExists("ForcingFunction"))
        {
            Pointer<CartGridFunction> f_fcn = new muParserCartGridFunction(
                "f_fcn", app_initializer->getComponentDatabase("ForcingFunction"), grid_geometry);
            time_integrator->registerBodyForceFunction(f_fcn);
        }

        // initialize gravitational force projector object.
        ForceProjector* ptr_gravityforce =
            new IBTK::ForceProjector("GravityForceProjector",
                                     ib_method_ops->getLDataManager(),
                                     patch_hierarchy,
                                     app_initializer->getComponentDatabase("ForceProjector"),
                                     ib_method_ops,
                                     "STAGGERED");
        ib_method_ops->registerPreProcessSolveFluidEquationsCallBackFunction(&callForceProjectorCallBackFunction,
                                                                             static_cast<void*>(ptr_gravityforce));

        // initialize Eulerian body force object.
        Pointer<CartGridFunction> ptr_cartgravityforce =
            new IBTK::CartGridBodyForce(ptr_gravityforce->getEulerianForcePatchDataIndex());
        time_integrator->registerBodyForceFunction(ptr_cartgravityforce);

        // Set up visualization plot file writers.
        Pointer<VisItDataWriter<NDIM> > visit_data_writer = app_initializer->getVisItDataWriter();
        Pointer<LSiloDataWriter> silo_data_writer = app_initializer->getLSiloDataWriter();
        if (uses_visit)
        {
            ib_initializer->registerLSiloDataWriter(silo_data_writer);
            ib_method_ops->registerLSiloDataWriter(silo_data_writer);
            time_integrator->registerVisItDataWriter(visit_data_writer);
        }

        // Initialize hierarchy configuration and data on all patches.
        time_integrator->initializePatchHierarchy(patch_hierarchy, gridding_algorithm);
        
        // Create ConstraintIBKinematics objects
        vector<Pointer<ConstraintIBKinematics> > ibkinematics_ops_vec;
        for (int i = 0; i < structure_names.size(); ++i)
        {
            Pointer<RigidBodyKinematics> ib_kinematics_op; 
            ib_kinematics_op = new RigidBodyKinematics(
                structure_names[i],
                app_initializer->getComponentDatabase("ConstraintIBKinematics")->getDatabase(structure_names[i]),
                ib_method_ops->getLDataManager(),
                patch_hierarchy);
            
            ibkinematics_ops_vec.push_back(ib_kinematics_op);
        }

        // register ConstraintIBKinematics objects with ConstraintIBMethod.
        ib_method_ops->registerConstraintIBKinematics(ibkinematics_ops_vec);
        ib_method_ops->initializeHierarchyOperatorsandData();
        // associate volume element with force projector.
        ptr_gravityforce->associateVolumeElement(ib_method_ops->getVolumeElement()[0]);
        
        // Deallocate initialization objects.
        ib_method_ops->freeLInitStrategy();
        ib_initializer.setNull();
        app_initializer.setNull();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);
        
        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            time_integrator->setupPlotData();
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            silo_data_writer->writePlotData(iteration_num, loop_time);
        }
        
        // inactivate redundant structures
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        int batch = input_db->getIntegerWithDefault("batch", 60);
        double frequency = input_db->getDoubleWithDefault("frequency",0.6);

        std::vector<int> initially_inactive_structures;
        for(int i = batch; i < num_structures - 1; i++){
            initially_inactive_structures.push_back(i);
        }
        ib_method_ops->getLDataManager()->inactivateLagrangianStructures(initially_inactive_structures,finest_ln);
        int row = 1;

        std::vector<bool> isactive(num_structures - 1, false) ;
        std::vector<int> activeParticles;
        for(int i = 0; i < batch; ++i)
        {
            isactive[i] = true;
            activeParticles.push_back(i);
        }

        std::vector<std::vector<bool>> isConstrainable(num_structures - 1, std::vector<bool>(50,false)); // last structure is PM
        using StructureParameters = ConstraintIBKinematics::StructureParameters;

        std::vector<std::vector<double>> old_particle_COM = ib_method_ops->getCurrentStructureCOM();
        std::vector<std::vector<double>> old_comVel = ib_method_ops->getCurrentCOMVelocity();

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        double dt = 0.0;
        std::ofstream clocktime("time.dat");
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            auto start = std::chrono::high_resolution_clock::now();

            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();
            

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            dt = time_integrator->getMaximumTimeStepSize();
            loop_time += dt;
            
            /**************************************************Activation******************************************************/
            // Activation of particles at the inlet
            if(loop_time > (frequency*row) && loop_time < (frequency*row + dt) && (row+1)*batch < num_structures){
              std::vector<int> activateThese;
              for(int particle_idx = row*batch; particle_idx < (row+1)*batch; particle_idx++){
                activateThese.push_back(particle_idx);
                isactive[particle_idx] = true;
                activeParticles.push_back(particle_idx);
              }
              ib_method_ops->getLDataManager()->activateLagrangianStructures(activateThese,finest_ln);
              activateThese.clear();
              row++;
            }
            /************************************************************************************************************************/

            /*************************************************Deactivation*********************************************************/
            // Deactivation of particles at the outlet
            std::vector<std::vector<double>> particle_COM = ib_method_ops->getCurrentStructureCOM();
            std::vector<int> deactivateThese;
            double particle_rad = input_db -> getDouble("particle_radius");
            // std::vector<double> acceleration(particle_COM.size()-1, 0.0) ; 
            double acceleration;
            std::vector<std::vector<double>> comVel = ib_method_ops->getCurrentCOMVelocity();

            for(auto particle_idx = activeParticles.begin(); particle_idx!=activeParticles.end();)
            {
                if((particle_COM[*particle_idx][0]+particle_rad) > 3.0)
                {
                    deactivateThese.push_back(*particle_idx);
                    isactive[*particle_idx] = false;
                    particle_idx = activeParticles.erase(particle_idx);
                }
                else
                {
                    double vel = std::sqrt(std::pow(comVel[*particle_idx][0],2) + std::pow(comVel[*particle_idx][1],2));
                    double old_vel = std::sqrt(std::pow(old_comVel[*particle_idx][0],2) + std::pow(old_comVel[*particle_idx][1],2));
                    acceleration = (vel - old_vel)/dt;
                    if(particle_COM[*particle_idx][0] > 0.5 && std::fabs(acceleration) < 0.001 && std::fabs(vel) < 0.001)
                    {
                        isConstrainable[*particle_idx].erase(isConstrainable[*particle_idx].begin());
                        isConstrainable[*particle_idx].push_back(true);
                        bool allTrue = true;
                        for(bool temp:isConstrainable[*particle_idx]) //if the particle was constrainable for 10 consecutive time steps, its motion will be constrained forever
                        {
                            if(!temp) 
                            {
                                allTrue = false;
                                break;
                            } 
                        }
                        if(allTrue)
                        {
                            StructureParameters& struct_param = ibkinematics_ops_vec[*particle_idx]->getStructureParameters();
                            struct_param.resetStructureIsSelfTranslating();
                            std::cout << *particle_idx << " constrained" << std::endl;
                        }
                    }
                    ++particle_idx;

                }
            }
            ib_method_ops->getLDataManager()->inactivateLagrangianStructures(deactivateThese, finest_ln);
            old_comVel = comVel;

            // for(unsigned int particle_idx = 0; particle_idx < particle_COM.size() - 1; particle_idx++){ 
                
            //     if(isactive[particle_idx] && particle_COM[particle_idx][0] > 0.5 && std::fabs(acceleration[particle_idx]) < 0.001 && std::fabs(vel) < 0.001)
            //     {
            //         isConstrainable[particle_idx].erase(isConstrainable[particle_idx].begin());
            //         isConstrainable[particle_idx].push_back(true);
            //         bool allTrue = true;
            //         for(bool temp:isConstrainable[particle_idx]) //if the particle was constrainable for 10 consecutive time steps, its motion will be constrained forever
            //         {
            //             if(!temp) 
            //             {
            //                 allTrue = false;
            //                 break;
            //             } 
            //         }
            //         if(allTrue)
            //         {
            //             StructureParameters& struct_param = ibkinematics_ops_vec[particle_idx]->getStructureParameters();
            //             struct_param.resetStructureIsSelfTranslating();
            //             std::cout << particle_idx << " constrained" << std::endl;
            //         }
            //     }
            // }            
            // ib_method_ops->getLDataManager()->inactivateLagrangianStructures(deactivateThese, finest_ln);
            // old_comVel = comVel;                    
            /************************************************************************************************************************/                   


            // Advance the hierarchy
            time_integrator->advanceHierarchy(dt);
            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // print out timer data, and store hierarchy data for post
            // processing.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                time_integrator->setupPlotData();
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
                silo_data_writer->writePlotData(iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\n";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
            if (dump_postproc_data && (iteration_num % postproc_data_dump_interval == 0 || last_step))
            {
                output_data(patch_hierarchy,
                            navier_stokes_integrator,
                            ib_method_ops->getLDataManager(),
                            iteration_num,
                            loop_time,
                            postproc_data_dump_dirname);
            }
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start);
            clocktime << iteration_num << "\t" << duration.count() << std::endl;
        }
        clocktime.close();

        // Cleanup Eulerian boundary condition specification objects (when
        // necessary).
        for (unsigned int d = 0; d < NDIM; ++d) delete u_bc_coefs[d];

    } // cleanup dynamically allocated objects prior to shutdown

    SAMRAIManager::shutdown();
    PetscFinalize();
    return 0;
}

void
output_data(Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
            Pointer<INSHierarchyIntegrator> navier_stokes_integrator,
            LDataManager* l_data_manager,
            const int iteration_num,
            const double loop_time,
            const string& data_dump_dirname)
{
    plog << "writing hierarchy data at iteration " << iteration_num << " to disk" << endl;
    plog << "simulation time is " << loop_time << endl;

    // Write Cartesian data.
    string file_name = data_dump_dirname + "/" + "hier_data.";
    char temp_buf[128];
    sprintf(temp_buf, "%05d.samrai.%05d", iteration_num, SAMRAI_MPI::getRank());
    file_name += temp_buf;
    Pointer<HDFDatabase> hier_db = new HDFDatabase("hier_db");
    hier_db->create(file_name);
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    ComponentSelector hier_data;
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getVelocityVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    hier_data.setFlag(var_db->mapVariableAndContextToIndex(navier_stokes_integrator->getPressureVariable(),
                                                           navier_stokes_integrator->getCurrentContext()));
    patch_hierarchy->putToDatabase(hier_db->putDatabase("PatchHierarchy"), hier_data);
    hier_db->putDouble("loop_time", loop_time);
    hier_db->putInteger("iteration_num", iteration_num);
    hier_db->close();

    // Write Lagrangian data.
    const int finest_hier_level = patch_hierarchy->getFinestLevelNumber();
    Pointer<LData> X_data = l_data_manager->getLData("X", finest_hier_level);
    Vec X_petsc_vec = X_data->getVec();
    Vec X_lag_vec;
    VecDuplicate(X_petsc_vec, &X_lag_vec);
    l_data_manager->scatterPETScToLagrangian(X_petsc_vec, X_lag_vec, finest_hier_level);
    file_name = data_dump_dirname + "/" + "X.";
    sprintf(temp_buf, "%05d", iteration_num);
    file_name += temp_buf;
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    VecView(X_lag_vec, viewer);
    PetscViewerDestroy(&viewer);
    VecDestroy(&X_lag_vec);
    return;
} // output_data
