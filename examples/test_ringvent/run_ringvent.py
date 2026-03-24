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
    gmsh.model.add("Ringvent_Case")

    exag = 1
    Geo2Gmsh.counter = 0
    Geo2Gmsh.counter1 = 0

    Geo2Gmsh.create_surface(
        surf_id=1,
        file_name="topo.txt",
        samx=136,
        samy=97,
        v_ex=exag,
        show_color=False
    )

    Geo2Gmsh.create_surface(
        surf_id=2,
        file_name="interface.txt",
        samx=136,
        samy=97,
        v_ex=exag,
        show_color=False
    )

    Geo2Gmsh.create_surface(
        surf_id=3,
        file_name="bottom.txt",
        samx=136,
        samy=97,
        v_ex=exag,
        show_color=False
    )

    volumes = Geo2Gmsh.volume_generation(
        num_loaded_surfaces = 3
        )

    wells = []
    for i in range(1, 10):
        well = Geo2Gmsh.add_well(
            file_name=f"well_{i}.txt",
            v_ex=exag,
            well_id=i,
        )
        wells.append(well)

    Geo2Gmsh.local_refinement(
        element_type = "well",
        element_list = wells[0]+wells[1]+wells[2]+wells[3]+wells[4]+wells[5]+wells[6]+wells[7]+wells[8],
        sampling =100,
        Size_Min =10,
        Size_Max =50,
        Dist_Min =150,
        Dist_Max =600
    )


    Geo2Gmsh.local_refinement(
        element_type = "surface",
        element_list = [2],
        sampling =800,
        Size_Min = 20,
        Size_Max = 50,
        Dist_Min = 60,
        Dist_Max = 100
    )


    #//////// block for field refinement ///////////
    Geo2Gmsh.counter += 1
    gmsh.model.mesh.field.add("Min", Geo2Gmsh.counter)
    gmsh.model.mesh.field.setNumbers(Geo2Gmsh.counter, "FieldsList", [field for field in range(2, Geo2Gmsh.counter, 2)])
    gmsh.model.mesh.field.setAsBackgroundMesh(Geo2Gmsh.counter)

    #///////////////////////////////////////////////
    gmsh.option.setNumber("Mesh.MeshSizeMax", 80 )
    gmsh.option.setNumber("Mesh.MeshSizeMin", 1)
    gmsh.option.setNumber("Mesh.Algorithm3D", 2)  # Change algoritm could work
    gmsh.option.setNumber("Geometry.Tolerance", 0.01)  # e.g., 0.01 meter

    #Physical entities

    Geo2Gmsh.physical_group(
        element_type = "well",
        element_list = wells[0]+wells[1]+wells[2]+wells[3]+wells[4]
    )

    Geo2Gmsh.physical_group(
        element_type = "well",
        element_list = wells[5]+wells[6]+wells[7]+wells[8]
    )

    Geo2Gmsh.physical_group(
        element_type = "surface",
        element_list = [1]
    )

    Geo2Gmsh.physical_group(
        element_type = "surface",
        element_list = [2]
    )

    Geo2Gmsh.physical_group(
        element_type = "volume",
        element_list = volumes
    )
    gmsh.model.geo.synchronize() #Super important. Not delete this line!!
    gmsh.model.mesh.generate(3)

    #///////// visulization (optional)
    gmsh.option.setNumber("Mesh.SurfaceEdges", 1)
    #gmsh.option.setNumber("Mesh.VolumeEdges", 1)
    gmsh.option.setNumber("Mesh.SurfaceFaces", 1)
    #gmsh.option.setNumber("Mesh.VolumeFaces", 1)  

    #///////// output format ////////////
    gmsh.write("outputs/test_ringvent.msh")
    gmsh.write("outputs/test_ringvent.vtk")


    #///// if working with sfepy the original vtk must be modified////
    with open("outputs/test_ringvent.vtk", "r") as f:
        lines = f.readlines()

    with open("outputs/test_ringvent.vtk", "w") as f:
        for line in lines:
            if "SCALARS CellEntityIds" in line:
                line = "SCALARS mat_id long\n"
            f.write(line)

    gmsh.fltk.run()
    gmsh.clear() 
    gmsh.finalize()

if __name__ == "__main__":
    main()

