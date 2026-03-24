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
    gmsh.model.add("Colombia_Case")

    # Here, we define a scalar factor and apply it to all functions.
    exag = 1

    # Global counters. DO NOT DELETE!
    Geo2Gmsh.counter = 0
    Geo2Gmsh.counter1 = 0

    ##/////////////-Section for functions and mesh settings-//////////////////

    Geo2Gmsh.create_surface(
        surf_id=1,
        file_name="topo.txt",
        samx=150,
        samy=115,
        v_ex=exag,
        show_color=False
    )

    Geo2Gmsh.create_surface(
        surf_id=2,
        file_name="sed_bottom.txt",
        samx=150,
        samy=115,
        v_ex=exag,
        show_color=False
    )

    Geo2Gmsh.create_surface(
        surf_id=3,
        file_name="moho.txt",
        samx=150,
        samy=115,
        v_ex=exag,
        show_color=False
    )

    volumes = Geo2Gmsh.volume_generation(
        num_loaded_surfaces = 3
        )

    fault_1 = Geo2Gmsh.add_fault(
        file_name="fault_1.txt",
        v_ex=exag,
        surf_id = 2,
        fault_id = 1,
        dip = 90.,
        dip_dir = 45.,
        fault_len = 10000.,
    )

    fault_2 = Geo2Gmsh.add_fault(
        file_name="fault_2.txt",
        v_ex=exag,
        surf_id = 2,
        fault_id = 2,
        dip = 90.,
        dip_dir = 45.,
        fault_len = 10000.,
    )

    wells = []
    for i in range(1, 12):
        well = Geo2Gmsh.add_well(
            file_name=f"well_{i}.txt",
            v_ex=exag,
            well_id=i
        )
        wells.append(well)

    Geo2Gmsh.local_refinement(
        element_type = "well",
        element_list = wells[0]+wells[1]+wells[2]+wells[3]+wells[4]+wells[5]+wells[6]+wells[7]+wells[8]+wells[9]+wells[10],
        sampling =100,
        Size_Min =500,
        Size_Max =8000,
        Dist_Min =5000,
        Dist_Max =30000
    )

    Geo2Gmsh.local_refinement(
        element_type = "fault",
        element_list = fault_1 + fault_2,
        sampling =800,
        Size_Min = 1000,
        Size_Max = 5000,
        Dist_Min = 10000,
        Dist_Max = 30000
    )

    Geo2Gmsh.local_refinement(
        element_type = "surface",
        element_list = [2],
        sampling =800,
        Size_Min = 2000,
        Size_Max = 5000,
        Dist_Min = 10000,
        Dist_Max = 40000
    )


    # #//////// Section for field refinement. Do not modify nor delete! ///////////
    Geo2Gmsh.counter += 1
    gmsh.model.mesh.field.add("Min", Geo2Gmsh.counter)
    gmsh.model.mesh.field.setNumbers(Geo2Gmsh.counter, "FieldsList", [field for field in range(2, Geo2Gmsh.counter, 2)])
    gmsh.model.mesh.field.setAsBackgroundMesh(Geo2Gmsh.counter)
    # #/////////////////////////////////////////////////////////////////////////////


    # #/////////////////////changes are optional//////////////////////////

    #Mesh size for areas where local refinement is not specified. In this case 15000 units.
    gmsh.option.setNumber("Mesh.MeshSizeMax", 15000 )
    # For certain errors changing the meshing algoritms can help. In this case algor. 2. (see Gmsh documentation)
    gmsh.option.setNumber("Mesh.Algorithm3D", 2)  
    #gmsh.option.setNumber("Geometry.Tolerance", 1000) # If this option is properly configured, geometries such as wedge layering can be represented
    #////////////////////////////////////////////////////////



    #////////////////physical group assignment/////////////////////

    Geo2Gmsh.physical_group(
        element_type = "well",
        element_list = wells[0]+wells[1]+wells[2]+wells[3]+wells[4]+wells[5]+wells[6]+wells[7]+wells[8]+wells[9]+wells[10]
    )

    Geo2Gmsh.physical_group(
        element_type = "fault",
        element_list = fault_1 + fault_2
    )

    Geo2Gmsh.physical_group(
        element_type = "surface",
        element_list = [1]
    )

    Geo2Gmsh.physical_group(
        element_type = "surface",
        element_list = [3]
    )

    Geo2Gmsh.physical_group(
        element_type = "volume",
        element_list = volumes
    )


    gmsh.model.geo.synchronize() #Super important. Not delete this line!!
    gmsh.model.mesh.generate(3) #Super important. Not delete this line!!

    #///////// visulization (optional)
    gmsh.option.setNumber("Mesh.SurfaceEdges", 1)
    #gmsh.option.setNumber("Mesh.VolumeEdges", 1)
    gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
    #gmsh.option.setNumber("Mesh.VolumeFaces", 1)  

    #///////// output format ////////////
    gmsh.write("outputs/test_llanos.msh")
    gmsh.write("outputs/test_llanos.vtk")


    #///// if working with sfepy the original vtk must be modified////
    with open("outputs/test_llanos.vtk", "r") as f:
        lines = f.readlines()

    with open("outputs/test_llanos.vtk", "w") as f:
        for line in lines:
            if "SCALARS CellEntityIds" in line:
                line = "SCALARS mat_id long\n"
            f.write(line)

    gmsh.fltk.run()
    gmsh.clear() 
    gmsh.finalize()

if __name__ == "__main__":
    main()

