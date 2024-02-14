# IBAMR
Some standard codes and tools

#1200 particles
Simulation of 1200 particles moving through a porous media in an inlet-outlet setup. Attempt has been made to define each particle with only one Lagrangian point, so the computational time is reduced. Grains of porous media still contain multiple L points. Code captures particles moving out of the domain.

#CallBackFunction
Usually IBAMR creates multiple output files for each structure (in Dump--...). These file contains information about structure's COM, transloational and rotational momentum etc. When number of structures increases beyond 143, this creates problem there will be lot files opened on the go. Also giving those many vertex files is cumbersome. This code creates the vertices of the structures during execution and also prints CoM data related to all the structures in one single file (SedimentingCylinde...). That data is then read using another script (Tools/extract.C)

