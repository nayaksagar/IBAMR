// Filename: ForceProjector.cpp
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

/////////////////////////////////////// INCLUDES ///////////////////////////////////////////

#include "ForceProjector.h"

#include "Box.h"
#include "Patch.h"
#include "PatchLevel.h"
#include "VariableDatabase.h"
#include "ibamr/namespaces.h"
#include "tbox/Array.h"
#include "tbox/PIO.h"

namespace IBTK
{
void
callForceProjectorCallBackFunction(const double current_time, const double new_time, const int cycle_num, void* ctx)
{
    static ForceProjector* ptr_forceprojector = static_cast<ForceProjector*>(ctx);
    if (cycle_num == 0)
    {
        ptr_forceprojector->calculateLagrangianBodyForce(new_time, current_time);
        ptr_forceprojector->calculateEulerianBodyForce(new_time, current_time);
    }
    return;

} // callForceProjectorCallBackFunction

ForceProjector::ForceProjector(const std::string& object_name,
                               LDataManager* lag_data_manager,
                               Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                               Pointer<Database> input_db,
                               ConstraintIBMethod* constraint_ib_method,
                               const std::string solver_type)
    : d_object_name(object_name),
      d_lag_data_manager(lag_data_manager),
      d_patch_hierarchy(patch_hierarchy),
      d_constraint_ib_method(constraint_ib_method),
      d_solver_type(solver_type),
      d_lag_pp_repulsion_force(),
      d_grav_const(NDIM)
{
    // put some default values.
    d_rho_fluid = 1.0;

    // Initialize  variables & variable contexts associated with Eulerian forces.
    VariableDatabase<NDIM>* var_db = VariableDatabase<NDIM>::getDatabase();
    d_body_force_context = var_db->getContext(d_object_name + "::BODYFORCE");
    if (d_solver_type == "STAGGERED")
        d_body_force_var = new SideVariable<NDIM, double>(d_object_name + "::BodyForce_sc_var");
    if (d_solver_type == "COLLOCATED")
        d_body_force_var = new CellVariable<NDIM, double>(d_object_name + "::BodyForce_cc_var", NDIM);
    d_body_force_idx = var_db->registerVariableAndContext(
        d_body_force_var, d_body_force_context, d_lag_data_manager->getGhostCellWidth());

    getFromInput(input_db);

    return;

} // ForceProjector

ForceProjector::~ForceProjector()
{
    // intentionally left blank
    return;

} //~ForceProjector

void
ForceProjector::getFromInput(Pointer<Database> input_db)
{
    d_rho_fluid = input_db->getDoubleWithDefault("rho_fluid", d_rho_fluid);
    d_grav_const = input_db->getDoubleArray("gravitational_constant");

    // Potential model
    d_eps_P = input_db->getDouble("eps_P");
    d_eps_W = input_db->getDouble("eps_W");
    d_xi = input_db->getDouble("xi");
    num_of_particles = input_db->getInteger("num_of_particles"); 
    particle_rad = input_db->getDouble("particle_radius"); //bug 

    return;

} // getFromInput

void
ForceProjector::registerLagrangianQuantityName(const std::string& lag_quantity_name)
{
    registerLagrangianQuantitiesName(std::vector<std::string>(1, lag_quantity_name));

    return;

} // registerLagrangianQuantityName

void
ForceProjector::registerLagrangianQuantitiesName(const std::vector<std::string>& lag_quantities_name)
{
    const unsigned size = lag_quantities_name.size();
    for (unsigned i = 0; i < size; ++i)
    {
        d_lag_quantities_name.push_back(lag_quantities_name[i]);
    }

    return;

} // registerLagrangianQuantitiesName

void
ForceProjector::associateVolumeElement(const double vol_lag_pt)
{
    d_vol_lag_pt = vol_lag_pt;

    return;

} // associateVolumeElement

void
ForceProjector::calculateLagrangianBodyForce(const double /*new_time*/, const double /*current_time*/)
{   
    //read grain_info data
    int num_grains;
    std::ifstream graindata;
    graindata.open("grain_info.txt");
    graindata >> num_grains;
    double grain_com_x[num_grains];
    double grain_com_y[num_grains];
    double grain_rad[num_grains];
    for(int i=0; i<num_grains; i++){
       graindata >> grain_com_x[i];
       graindata >> grain_com_y[i];
       graindata >> grain_rad[i];
    }
    graindata.close();

    const int coarsest_ln = 0;
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    d_lag_force.clear();
    d_lag_pp_repulsion_force.clear();

    d_lag_force.resize(finest_ln + 1);
    d_lag_pp_repulsion_force.resize(finest_ln + 1);
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_lag_data_manager->levelContainsLagrangianData(ln)) continue;
        d_lag_force[ln] = d_lag_data_manager->createLData(d_object_name + "::lag_force_data", ln, NDIM, false);
        d_lag_pp_repulsion_force[ln] =
            d_lag_data_manager->createLData(d_object_name + "::lag_pp_repulsion_force_data", ln, NDIM, false);
    }

    // Get the centers of mass of the structures
    std::vector<std::vector<double> > structure_COM = d_constraint_ib_method->getCurrentStructureCOM();

    // Positions of each particle
    std::vector<Pointer<LData> > x_coord(finest_ln + 1, Pointer<LData>(NULL));
    x_coord[finest_ln] = d_lag_data_manager->getLData("X", finest_ln);

    // Set Lagrangian gravitational force.
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        if (!d_lag_data_manager->levelContainsLagrangianData(ln)) continue;

        // Get pointer to LData corresponding to lagrangian force.
        boost::multi_array_ref<double, 2>& X_data = *x_coord[ln]->getLocalFormVecArray();
        boost::multi_array_ref<double, 2>& F_data = *d_lag_force[ln]->getLocalFormVecArray();
        const Pointer<LMesh> mesh = d_lag_data_manager->getLMesh(ln);
        const std::vector<LNode*>& local_nodes = mesh->getLocalNodes();

        for (std::vector<LNode*>::const_iterator cit = local_nodes.begin(); cit != local_nodes.end(); ++cit)
        {
            const LNode* const node_idx = *cit;
            const int local_idx = node_idx->getLocalPETScIndex();
            const double* const X = &X_data[local_idx][0];
            double* const F = &F_data[local_idx][0];
            
            IBTK::Vector3d particle_COM;
              
            for(int particle_index=0; particle_index<num_of_particles; particle_index++){                
                for (int d = 0; d < 3; ++d) particle_COM[d] = structure_COM[particle_index][d];
                
                //check whether lagrangian point lies within the particle
                bool within_particle = false;
                if(std::pow(X[0] - particle_COM[0], 2.0) + std::pow(X[1] - particle_COM[1], 2.0) <= 
                        std::pow(particle_rad + d_xi/4.0, 2.0)){
                    within_particle = true;
                }
                
                IBTK::Vector3d Fpg, Fpp, Fpw;

                if(within_particle){
                    
                    for (int d = 0; d < NDIM; ++d) F[d] = d_rho_fluid * d_grav_const[d] * d_vol_lag_pt;

                    //if it lies within the particle calculate force acting on it due to neighboring grains
                    if(num_grains > 0){
                        for(int grain_index=0; grain_index<num_grains; grain_index++){
                            IBTK::Vector3d grain_COM;
                            grain_COM[0] = grain_com_x[grain_index];
                            grain_COM[1] = grain_com_y[grain_index];
                            Fpg  = computeParticleGrainForce(particle_COM, grain_COM, grain_rad[grain_index]);
                            for (int d = 0; d < NDIM; ++d) F[d] += Fpg[d] * d_vol_lag_pt;
                            // std::cout<<"grain force being computed"<<std::endl;
                        }
                    }

                    //if it lies within the particle calculate force acting on it due to neighboring particles
                    if(num_of_particles > 1){
                        // std::cout<<"no. of particles = "<< num_of_particles<<std::endl;
                        for(int neigh_part_index=0; neigh_part_index<num_of_particles; neigh_part_index++){
                            if (neigh_part_index == particle_index)
                            continue;
                            else{
                                IBTK::Vector3d neigh_part_COM;
                                neigh_part_COM[0] = structure_COM[neigh_part_index][0];
                                neigh_part_COM[1] = structure_COM[neigh_part_index][1];
                                Fpp  = computeParticleParticleForce(particle_COM, neigh_part_COM); //bug: all particles are assumed to be of same radius
                                for (int d = 0; d < NDIM; ++d) F[d] += Fpp[d] * d_vol_lag_pt;
                                // std::cout<<"particle force being computed"<<std::endl;
                            }
                        }
                    }
                
                    //check if particle is near any wall, if yes, which wall (left, right or bottom)
                    double left = 0.0, right = 2.0, bottom = 0.0, wall_pos; //bug:top wall is not modelled
                    int wall;
                    bool near_wall = false;
                    if(fabs(left - particle_COM[0]) < particle_rad + d_xi){
                         wall_pos = left;
                         wall = 0;
                         near_wall = true;
                    }
                    else if(fabs(right - particle_COM[0]) < particle_rad + d_xi){
                        wall_pos = right;
                        wall = 0;
                        near_wall = true;
                    }
                    else if(fabs(bottom - particle_COM[1]) < particle_rad + d_xi){
                         wall_pos = bottom;
                         wall = 1;
                         near_wall = true;
                    }
                    if(near_wall){
                        Fpw = computeParticleWallForce(particle_COM, wall_pos, wall);
                        for (int d = 0; d < NDIM; ++d) F[d] += Fpw[d] * d_vol_lag_pt;
                    }

                }
                
            }
        }
        x_coord[ln]->restoreArrays();
        d_lag_force[ln]->restoreArrays();
        d_lag_pp_repulsion_force[ln]->restoreArrays();
    } // all levels

    return;

} // calculateLagrangianBodyForce

void
ForceProjector::calculateEulerianBodyForce(const double /*new_time*/, const double current_time)
{
    // allocate patch data for Eulerian forcing.
    const int coarsest_ln = 0;
    const int finest_ln = d_patch_hierarchy->getFinestLevelNumber();
    for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
    {
        Pointer<PatchLevel<NDIM> > level = d_patch_hierarchy->getPatchLevel(ln);
        if (level->checkAllocated(d_body_force_idx)) level->deallocatePatchData(d_body_force_idx);
        level->allocatePatchData(d_body_force_idx, current_time);

        for (PatchLevel<NDIM>::Iterator p(level); p; p++)
        {
            Pointer<Patch<NDIM> > patch = level->getPatch(p());
            if (d_solver_type == "STAGGERED")
            {
                Pointer<SideData<NDIM, double> > body_force_data = patch->getPatchData(d_body_force_idx);
                body_force_data->fill(0.0);
            }
            else if (d_solver_type == "COLLOCATED")
            {
                Pointer<CellData<NDIM, double> > body_force_data = patch->getPatchData(d_body_force_idx);
                body_force_data->fill(0.0);
            }
            else
            {
                TBOX_ERROR("ForceProjector::calculateEulerianBodyForce() "
                           << "UNKNOWN SOLVER ENCOUNTERED"
                           << std::endl);
            }

        } // iterate over patches

    } // all levels.

    // spread the lagrangian force from finest level to the finest level.
    std::vector<Pointer<LData> > F_data(finest_ln + 1, Pointer<LData>(NULL));
    std::vector<Pointer<LData> > X_data(finest_ln + 1, Pointer<LData>(NULL));

    // Fill in the above vectors at the finest level.
    F_data[finest_ln] = d_lag_force[finest_ln];
    X_data[finest_ln] = d_lag_data_manager->getLData("X", finest_ln);

    // Spread the deformation velocities.
    d_lag_data_manager->spread(d_body_force_idx, F_data, X_data, (RobinPhysBdryPatchStrategy*)NULL);

    return;

} // calculateEulerianBodyForce

IBTK::Vector3d
ForceProjector::computeParticleParticleForce(IBTK::Vector3d particle_COM, IBTK::Vector3d neigh_part_COM)
{
    IBTK::Vector3d Fpp;

    IBTK::Vector3d COM_diff = particle_COM - neigh_part_COM;
    const double dij = std::sqrt(COM_diff.dot(COM_diff));
    if (dij > particle_rad + particle_rad + d_xi)  //bug: particle_rad and neigh_part_rad are same
    {
        Fpp.setZero();
    }
    else if (dij <= particle_rad + particle_rad + d_xi)
    {
        d_cij = d_rho_fluid*3.14159265*particle_rad*particle_rad*980.0;
        Fpp = d_cij / d_eps_P * pow((dij - particle_rad - particle_rad - d_xi) / d_xi, 2.0) * (COM_diff) / dij;
    }

    return Fpp;
}

IBTK::Vector3d
ForceProjector::computeParticleGrainForce(IBTK::Vector3d particle_COM, IBTK::Vector3d grain_COM, double grain_rad)
{
    IBTK::Vector3d Fpg;
    
    IBTK::Vector3d COM_diff = particle_COM - grain_COM;
    const double dij = std::sqrt(COM_diff.dot(COM_diff));
    if (dij > particle_rad + grain_rad + d_xi)
    {
        Fpg.setZero();
    }
    else if (dij <= particle_rad + grain_rad + d_xi)
    {
        d_cij = d_rho_fluid*3.14159265*particle_rad*grain_rad*980.0;
        Fpg = d_cij / d_eps_P * pow((dij - particle_rad - grain_rad- d_xi) / d_xi, 2.0) * (COM_diff) / dij;
    }

    return Fpg;
}

IBTK::Vector3d
ForceProjector::computeParticleWallForce(IBTK::Vector3d particle_COM, double wall_pos, int wall){
    IBTK::Vector3d Fpw;
    
    double COM_diff = 2.0*(particle_COM[wall] - wall_pos);
    double dij = std::sqrt(COM_diff*COM_diff);
    if (dij > 2.0*particle_rad + d_xi)
    {
        Fpw.setZero();
    }
    else if (dij <= 2.0*particle_rad + d_xi)
    {
        for(int d=0;d<3;d++){
          if(d==wall){
            d_cij = d_rho_fluid*3.14159265*particle_rad*980.0;
            Fpw[d] = d_cij / d_eps_W * pow((dij - 2.0*particle_rad - d_xi) / d_xi, 2.0) * COM_diff/dij;
          }
          else{
            Fpw[d] = 0.0;
          }
        }
    }

    return Fpw;
}
} // namespace IBTK
