#!/usr/bin/env python
# coding: utf-8

# In[4]:


import gmsh
import Geo2Gmsh

def main():

    gmsh.initialize()
    gmsh.model.add("Overview_example")

    # Global counters. DO NOT DELETE!
    Geo2Gmsh.counter = 0
    Geo2Gmsh.counter1 = 0
    
    
    ##/////////////-Section for functions and mesh settings-//////////////////
    
    Geo2Gmsh.create_surface(
        surf_id=1,
        file_name="layer_1.txt",
        samx=50,
        samy=50,
        v_ex=1,
        show_color=False
    )
    
    Geo2Gmsh.create_surface(
        surf_id=2,
        file_name="layer_2.txt",
        samx=50,
        samy=50,
        v_ex=1,
        show_color=False
    )
    
    volumes = Geo2Gmsh.volume_generation(
        num_loaded_surfaces = 2
    )
    
    well_1 = Geo2Gmsh.add_well(
        file_name="well_1.txt",
        v_ex=1,
        well_id=1
    )
    
    fault_4 = Geo2Gmsh.add_fault(
        file_name="fault_4.txt",
        v_ex=1,
        surf_id = 1,
        fault_id = 4,
        dip = 55.,
        dip_dir = 310.,
        fault_len = 5.
    )
    
    Geo2Gmsh.local_refinement(
        element_type = "well",
        element_list = well_1,
        sampling =50,
        Size_Min =0.3,
        Size_Max =0.8,
        Dist_Min =2,
        Dist_Max =3.
    )
    
    Geo2Gmsh.local_refinement(
        element_type = "fault",
        element_list = fault_4,
        sampling =1000,
        Size_Min = 0.8,
        Size_Max = 1.,
        Dist_Min = 1,
        Dist_Max = 2
    )
    
    Geo2Gmsh.local_refinement(
        element_type = "surface",
        element_list = [2],
        sampling =1000,
        Size_Min = 0.4,
        Size_Max = 0.8,
        Dist_Min = 2.,
        Dist_Max = 3.
    )
    
    # #//////// Section for field refinement. Do not modify nor delete! ///////////
    Geo2Gmsh.counter += 1 
    gmsh.model.mesh.field.add("Min", Geo2Gmsh.counter)
    gmsh.model.mesh.field.setNumbers(Geo2Gmsh.counter, "FieldsList", [field for field in range(2, Geo2Gmsh.counter, 2)])
    gmsh.model.mesh.field.setAsBackgroundMesh(Geo2Gmsh.counter)
    # #/////////////////////////////////////////////////////////////////////////////
    
    
    # #/////////////////////changes are optional//////////////////////////
    
    #Mesh size for areas where local refinement is not specified. In this case 2 units.
    gmsh.option.setNumber("Mesh.MeshSizeMax", 2 )
    # For certain errors changing the meshing algoritms can help. In this case algor. 2. (see Gmsh documentation)
    gmsh.option.setNumber("Mesh.Algorithm3D", 2)  
    #////////////////////////////////////////////////////////
    
    #////////////////////physical group assignment////////////////////////////
    
    Geo2Gmsh.physical_group(
        element_type = "well",
        element_list = well_1
    )
    
    Geo2Gmsh.physical_group(
        element_type = "fault",
        element_list = fault_4
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
    
    
    #///////// output format ////////////
    gmsh.write("outputs/Example.msh")
    gmsh.write("outputs/Example.vtk")
    
    
    #///// if working with sfepy the original vtk must be modified////
    with open("outputs/Example.vtk", "r") as f:
        lines = f.readlines()
    
    with open("outputs/Example.vtk", "w") as f:
        for line in lines:
            if "SCALARS CellEntityIds" in line:
                line = "SCALARS mat_id long\n"
            f.write(line)
            
    gmsh.clear() 
    gmsh.finalize()

if __name__ == "__main__":
    main()


# In[ ]:




