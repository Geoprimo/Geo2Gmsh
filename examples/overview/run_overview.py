#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import gmsh
import sys
import os
sys.path.append(os.path.abspath("../.."))

import Geo2Gmsh

def main():

    ##///////////////////////////////////////start///////////////////////////////

    gmsh.initialize()
    gmsh.clear() 
    gmsh.model.add("Overview_example")

    exag = 1
    Geo2Gmsh.counter = 0
    Geo2Gmsh.counter1 = 0

    Geo2Gmsh.create_surface(
        surf_id=1,
        file_name="layer_1.txt",
        samx=448,
        samy=448,
        v_ex=1,
        show_color=False
    )

    Geo2Gmsh.create_surface(
        surf_id=2,
        file_name="layer_2.txt",
        samx=448,
        samy=448,
        v_ex=1,
        show_color=False
    )

    volumes = Geo2Gmsh.volume_generation(
        num_loaded_surfaces = 2
    )

    fault_1 = Geo2Gmsh.add_fault(
        file_name="fault_1.txt",
        v_ex=1,
        surf_id = 1,
        fault_id = 1,
        dip = 65.,
        dip_dir = 280.,
        fault_len = 5.
    )

    well_1 = Geo2Gmsh.add_well(
        file_name="well_1.txt",
        v_ex=1,
        well_id=1
    )

    Geo2Gmsh.local_refinement(
        element_type = "well",
        element_list = well_1,
        sampling =100,
        Size_Min =0.5,
        Size_Max =2.0,
        Dist_Min =3,
        Dist_Max =5
    )


    Geo2Gmsh.local_refinement(
        element_type = "fault",
        element_list = fault_1,
        sampling =100,
        Size_Min = 0.8,
        Size_Max = 4.,
        Dist_Min = 3,
        Dist_Max = 5
    )

    Geo2Gmsh.local_refinement(
        element_type = "surface",
        element_list = [2],
        sampling =1000,
        Size_Min = 1.2,
        Size_Max = 4,
        Dist_Min = 2.,
        Dist_Max = 3.
    )


    # #//////// Section for field refinement /////////// 
    Geo2Gmsh.counter += 1
    gmsh.model.mesh.field.add("Min", Geo2Gmsh.counter)
    gmsh.model.mesh.field.setNumbers(Geo2Gmsh.counter, "FieldsList", [field for field in range(2, Geo2Gmsh.counter, 2)])
    gmsh.model.mesh.field.setAsBackgroundMesh(Geo2Gmsh.counter)

    # #/////////////////////////////////////////////////////////////////////////////


    # #/////////////////////changes are optional//////////////////////////

    #Mesh size for areas where local refinement is not specified. In this case 2 units.
    gmsh.option.setNumber("Mesh.MeshSizeMin",0.4)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 4. )
    # For certain errors changing the meshing algoritms can help. In this case algor. 2. (see Gmsh documentation)
    gmsh.option.setNumber("Mesh.Algorithm3D", 2)  
    #////////////////////////////////////////////////////////
    gmsh.option.setNumber("Geometry.Tolerance", 0.4) #This allows Gmsh recognize close points (distance>0.1) as different nodes.

    #////////////////////physical group assignment////////////////////////////

    Geo2Gmsh.physical_group(
        element_type = "well",
        element_list = well_1
    )

    Geo2Gmsh.physical_group(
        element_type = "fault",
        element_list = fault_1
    )

    Geo2Gmsh.physical_group(
        element_type = "surface",
        element_list = [2]
    )

    Geo2Gmsh.physical_group(
        element_type = "volume",
        element_list = volumes)


    gmsh.model.geo.synchronize() #Super important. Not delete this line!!
    gmsh.model.mesh.generate(3) #Super important. Not delete this line!!

    #///////// visulization (optional) ////////////////////
    gmsh.option.setNumber("Mesh.SurfaceEdges", 1)
    #gmsh.option.setNumber("Mesh.VolumeEdges", 1)
    #gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
    gmsh.option.setNumber("Mesh.VolumeFaces", 1)  

    #///////// output format ////////////
    gmsh.write("outputs/overview.msh")
    gmsh.write("outputs/overview.vtk")


    #///// if working with sfepy the original vtk must be modified////
    with open("outputs/overview.vtk", "r") as f:
        lines = f.readlines()

    with open("outputs/overview.vtk", "w") as f:
        for line in lines:
            if "SCALARS CellEntityIds" in line:
                line = "SCALARS mat_id long\n"
            f.write(line)

    gmsh.fltk.run()
    gmsh.clear() 
    gmsh.finalize()

if __name__ == "__main__":
    main()


# In[ ]:




