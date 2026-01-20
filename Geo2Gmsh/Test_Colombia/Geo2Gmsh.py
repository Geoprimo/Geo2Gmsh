#!/usr/bin/env python
# coding: utf-8

import numpy as np
import gmsh

counter = 0
conunter1 = 0

def create_surface(
    surf_id: int,
    file_name: str,
    samx: int,
    samy: int,
    v_ex: float = 1,
    show_color:bool = False
):
    """
    Create a surface mesh in Gmsh from a text file with x, y, z coordinates.
    txt struc:

    |   id    |    x    |    y    |    z    |
    
    n rows x 4 columns!

    Parameters:
        surf_id (int): Surface identifier, used for tagging.
        file_name (str): Path to text file with x,y,z data.
        samx (int): Number of columns in the sampled grid (along X-axis).
        samy (int): Number of rows in the sampled grid (along Y-axis).
        v_ex (float): Vertical exaggeration factor applied to z-values.
        show_color (bool): Whether to display the surface with color mapping based on z values.
    """
    Nx = samx - 1
    Ny = samy -1
    
    # Load data and extract dimensions:
    data = np.loadtxt(f"layers/{file_name}", skiprows=1)
    id = (surf_id-1) * 10


    # Defining extend from txt file
    min_x = np.min(data[:, 1])
    max_x = np.max(data[:, 1])
    min_y = np.min(data[:, 2])
    max_y = np.max(data[:, 2])

    
    z = data[:, 3] * v_ex
    Z = np.reshape(z, [samy, samx])

    # Creacion de arrays en ejes x e y
    x = np.linspace(min_x, max_x, num=samx, endpoint=True)   
    y = np.linspace(min_y, max_y, num=samy, endpoint=True)

    # GRID based on data:
    X, Y = np.meshgrid(x, y, indexing='xy')
    Y = Y[::-1, :]
    
    # Helper function to return a node tag given indices i and j:
    def tag(i, j):
        return ((Ny + 1) * i + j + 1) 
        
    # Coordinates list
    coords = []
    # tag list
    nodes = []
    # traingular-element connectivity  (3 node tags per triangle):
    tris = []
    # line element connectivity on the four boundaries
    lin = [[], [], [], []]
    
    # point element connectivity at the four corners
    pnt = [tag(0, Ny), tag(Nx, Ny), tag(Nx, 0), tag(0, 0)]

    # Populate the previously created vectors:
    for i in range( Nx + 1):    #column
        for j in range( Ny, -1, -1 ):      #row
            nodes.append(tag(i, j))
            coords.extend([X[j, i], Y[j, i], Z[j, i]])
            if i > 0 and j > 0:
                tris.extend([tag(i - 1, j - 1), tag(i, j - 1), tag(i - 1, j)])
                tris.extend([tag(i, j - 1), tag(i, j), tag(i - 1, j)])
            if (i == 0 or i == Nx) and j > 0:
                lin[3 if i == 0 else 1].extend([tag(i, j - 1), tag(i, j)])
            if (j == 0 or j == Ny) and i > 0:
                lin[0 if j == 0 else 2].extend([tag(i - 1, j), tag(i, j)])
                
    # Create 4 discrete points for the four corners of the terrain:
    for i in range(4):
        gmsh.model.addDiscreteEntity(0, i + 1 + id)
        
    # Modify the coordinates of the points created earlier
    gmsh.model.setCoordinates(1 + id, min_x, min_y, coords[3 * tag(0, 0) - 1])
    gmsh.model.setCoordinates(2 + id, max_x, min_y, coords[3 * tag(Nx, 0) - 1])
    gmsh.model.setCoordinates(3 + id, max_x, max_y, coords[3 * tag(Nx, Ny) - 1])
    gmsh.model.setCoordinates(4 + id, min_x, max_y, coords[3 * tag(0, Ny) - 1])

    # Create 4 boundary curves of the terrain using the previous points:
    for i in range(4):
        gmsh.model.addDiscreteEntity(1, i + 1 + id, [i + 1 + id, i + 2 + id if i < 3 else 1 + id])

    # Create a surface using the previously created curves:
    gmsh.model.addDiscreteEntity(2, 1 + id, [1 + id, 2 + id, 3 + id, 4 + id])

    # Add all nodes to the surface:
    gmsh.model.mesh.addNodes(2, 1 + id, nodes, coords)

    for i in range(4):
        # point elements
        gmsh.model.mesh.addElementsByType(i + 1 + id, 15, [], [pnt[i]])
        # line elements (2 nodes)
        gmsh.model.mesh.addElementsByType(i + 1 + id, 1, [], lin[i])
    # triangle elements (3 nodos)
    gmsh.model.mesh.addElementsByType(1 + id, 2, [], tris)


    
    # Generate a color-mapped view based on Z (elevation)

    if show_color == True:
        view_tag = gmsh.view.add(f"height_{surf_id}")
        z_only = [[z] for z in coords[2::3]]  # xtract every third value from coords = Z
        gmsh.view.addModelData(
            view_tag,
            0,                          # time step. Do not delete
            gmsh.model.getCurrent(),   # current model name
            "NodeData",                # data type
            nodes,                     # nodes tag (list)
            z_only                     # Scalar values per node (Z)
    )

    print(f"Superficie {surf_id } successfully created.")

def volume_generation(
    num_loaded_surfaces: int = 2
):
    """
    Generates 3D volumes from a stack of loaded geological surfaces. Returns a list of the created volumes.

    Parameters:
    -----------
    num_loaded_surfaces : int
        Total number of 2D surfaces loaded into the model.
        Must be >= 2. One volume is created between each pair of consecutive surfaces.
    """
    # Reclassify the nodes on the curves and points (since we previously added them all to the surface using 'addNodes' for simplicity)
    gmsh.model.mesh.reclassifyNodes()

    # Create a geometry from the discrete curves and surface, so that they can be remeshed later:
    gmsh.model.mesh.createGeometry()

    vertical_union = {}
    surfaces = {}
    volume = []

    for j in range (num_loaded_surfaces - 1):  # Number of surfaces minus 1
        for i in range(1,5):
            #Creating vertical lines conection between each surface
            k= i + j * 10
            tag = gmsh.model.geo.addLine(k, k + 10)
            vertical_union[f"corner_{k}"] = tag
        
        #Creating lateral surfaces 
        for i in range(1,5):
            k = (i + j * 10)
            l = np.array([3 , 4, 2, 3, 1, 2, 4, 1]) + (j * 10)
            l = l.reshape(-1,2)
            loop = gmsh.model.geo.addCurveLoop([-k, -vertical_union[f"corner_{l[i-1,0]}"], (k + 10), vertical_union[f"corner_{l[i-1,1]}"] ])
            tag = gmsh.model.geo.addPlaneSurface([loop])
            surfaces [f"surface_{k}"] = tag

    #creating the volume
    for i in range (num_loaded_surfaces - 1):
        loop = gmsh.model.geo.addSurfaceLoop([1 + 10 * (i+1), surfaces[f"surface_{10*(i) + 1}"], surfaces[f"surface_{10*(i) + 2}"], 
                                              surfaces[f"surface_{10*(i) + 3}"], surfaces[f"surface_{10*(i) + 4}"], 10*(i) + 1])
        tag = gmsh.model.geo.addVolume([loop])
        volume.append(tag)
        print(f"Volumes {i+1} succesfully created.")
    return volume

def add_well(
    file_name: str,
    v_ex: float = 1,
    surf_id: int = 1,
    well_id:int = 1
):

    """
    Create a well with a trajectory defined by a txt file.
    txt contains points along the well. Returns a list of lines tags belonging to the well.
    txt struc:

    |   id    |    x    |    y    |    z    |
     
    n rows x 4 columns!
    
    Parameters:
        file_name (str): Path to text file with x,y,z data.
        v_ex (float): Vertical exaggeration factor.
        surf_id: ID of the top surface (e.g., topography). This is the surface that contains the top of the well.  
        well_id: ID of the well user-defined.
    """
    # Text file with x, y, z for each point belonging to the well trajectory. 
    well_txt = np.loadtxt(f"wells/{file_name}", skiprows=1)

    x = well_txt[:, 1]  
    y = well_txt[:, 2]
    z = well_txt[:, 3] * v_ex
    
    # Create a list using well points
    well_points  = []
    for xi, yi, zi in zip(x, y, z):
        tag = gmsh.model.geo.addPoint(xi,yi,zi)
        well_points.append(tag)
        
    # Generate lines by connecting points from the list well_points
    well_line  = []
    for i in range(len(well_points) - 1):
        tag = gmsh.model.geo.addLine(well_points[i], well_points[i+1])
        well_line.append(tag)

    gmsh.model.geo.synchronize()
    
    # Embedding points and lines to the mesh
    gmsh.model.mesh.embed(0, [well_points[0]], 2, 1+10*(surf_id-1))
    gmsh.model.mesh.embed(1, well_line, 3, surf_id)
    gmsh.model.geo.synchronize()
    print(f"Well {well_id} succesfully created.")
    return well_line
    
    
def add_fault(
    file_name: str,
    v_ex: float = 1,
    surf_id: int = 1,
    fault_id: int = 1,
    dip: float = 90.,
    dip_dir: float = 90.,
    fault_len: float = 1000.,
):
    """
    Create a fault that intercepts one specified surface only. Returns a list of surface tags belonging to the fault.
    txt struc:

    |   id    |    x    |    y    |    z    |
     
    n rows x 4 columns!

    Parameters:
        file_name (str): Path to text file with x,y,z data.
        v_ex (float): Vertical exaggeration factor.
        surf_id (int): ID of the top surface (e.g., topography). This is the surface that is being intercepted by the fault.
        fault_id (int): ID of the fault user-defined.
        dip (float): dip, Common parameter in geology for fault description.
        dip_dir (float): dip direction, Common parameter in geology for fault description.
        fault_len (float): fault length.
    """
    # Loading TXT file of the fault curve at surface
    fault_txt = np.loadtxt(f"faults/{file_name}", skiprows=1)

    x = fault_txt[:, 1]  
    y = fault_txt[:, 2]
    z_top = fault_txt[:, 3] * v_ex

    # Add points of the fault curve at surface
    top_points  = []
    for xi, yi, zi in zip(x, y, z_top):
        tag = gmsh.model.geo.addPoint(xi,yi,zi)
        top_points.append(tag)

    line_top = []

    # Generate a line by connecting points from the list top_points
    for i in range(len(top_points) - 1):
        line_tag = gmsh.model.geo.addLine(top_points[i], top_points[i+1])
        line_top.append(line_tag)

    # Here extrusion capabilities are used to generate the fault surface
    to_extrude = [(1, tag) for tag in line_top]


    #Calculating parameters required for extrusion from dip/dip direction and length of faults
    
    #vertical component dz
    dz = - (np.sin(np.radians(dip)) * fault_len)
    
    #Calculating horizontal component of vector: (dx,dy)
    hor_comp = np.cos(np.radians(dip)) * fault_len


    if 0 <= dip_dir <= 90:
        alpha = np.radians(dip_dir)
        dx = np.sin(alpha) * hor_comp
        dy = np.cos(alpha) * hor_comp
    elif 90 < dip_dir <= 180:
        alpha = np.radians(180 - dip_dir)
        dx = np.sin(alpha) * hor_comp
        dy = - (np.cos(alpha) * hor_comp)
    elif 180 < dip_dir <= 270:
        alpha = np.radians(dip_dir - 180)
        dx = - (np.sin(alpha) * hor_comp)
        dy = - (np.cos(alpha) * hor_comp)
    elif 270 < dip_dir < 360:
        alpha = np.radians(360 - dip_dir)
        dx = - (np.sin(alpha) * hor_comp)
        dy =  np.cos(alpha) * hor_comp
    else:
        print("dip_dir should be 0-359°")
        
    
    # Linear extrusion:
    extruded = gmsh.model.geo.extrude(to_extrude, dx, dy, dz)
    sect_fault = [e[1] for e in extruded if e[0] == 2]

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(0, top_points, 2, 1+10*(surf_id-1))
    gmsh.model.mesh.embed(1, line_top, 2, 1+10*(surf_id-1))
    gmsh.model.mesh.embed(2, sect_fault, 3, surf_id)
    gmsh.model.geo.synchronize()
    print(f"Fault {fault_id} succesfully created.")
    return sect_fault

def local_refinement(
    element_type: str,
    element_list: list,
    sampling: int = 100,
    Size_Min: int = 6000,
    Size_Max: int = 25000,
    Dist_Min: int = 10000,
    Dist_Max: int = 30000
):
    """
    Gmsh supports the generation of locally refined meshes, the local_refinement function integrates these capabilities and facilitates their application to the     surfaces, wells and faults previously loaded.

    Parameters:
        element_type (str): Specify the type of element. Choose from “well”, “fault” or “surface”. 
        element_list (list): For wells and faults, provide the list returned by the corresponding function. For surfaces, specify the surface ID (surf_id)                                    enclosed in square brackets (e.g., [1]). 
        sampling (int): Defines how many points are sampled along the feature’s geometry to guide local mesh refinement.  
        Size_Min (int): Specifies the minimum element size allowed during mesh refinement. 
        Size_Max (int): Specifies the maximum element size allowed during mesh generation. 
        Dist_Min (int): Defines the minimum radius around a feature within which the specified Size_Min is strictly enforced.
        Dist_Max (int): Defines the maximum distance (radius) over which mesh refinement can influence element size. 
    """
    gmsh.model.mesh.reclassifyNodes()
    gmsh.model.geo.removeAllDuplicates()
    gmsh.model.geo.synchronize()
    
    global counter
    
    if element_type.lower() == "well":
        
        gmsh.model.mesh.field.add("Distance", counter + 1)
        gmsh.model.mesh.field.setNumbers(counter + 1, "CurvesList", element_list)
        gmsh.model.mesh.field.setNumber(counter + 1, "Sampling", sampling)

        gmsh.model.mesh.field.add("Threshold", counter + 2)
        gmsh.model.mesh.field.setNumber(counter + 2, "InField", counter + 1)
        gmsh.model.mesh.field.setNumber(counter + 2, "SizeMin", Size_Min)
        gmsh.model.mesh.field.setNumber(counter + 2, "SizeMax", Size_Max)
        gmsh.model.mesh.field.setNumber(counter + 2, "DistMin", Dist_Min)
        gmsh.model.mesh.field.setNumber(counter + 2, "DistMax", Dist_Max)
        counter+=2

        print(f"Mesh succesfully refined around wells.")

    elif element_type.lower() == "fault":

        gmsh.model.mesh.field.add("Distance", counter + 1)
        gmsh.model.mesh.field.setNumbers(counter + 1, "SurfacesList", element_list)
        gmsh.model.mesh.field.setNumber(counter + 1, "Sampling", sampling)

        gmsh.model.mesh.field.add("Threshold", counter + 2)
        gmsh.model.mesh.field.setNumber(counter + 2, "InField", counter + 1)
        gmsh.model.mesh.field.setNumber(counter + 2, "SizeMin", Size_Min)
        gmsh.model.mesh.field.setNumber(counter + 2, "SizeMax", Size_Max)
        gmsh.model.mesh.field.setNumber(counter + 2, "DistMin", Dist_Min)
        gmsh.model.mesh.field.setNumber(counter + 2, "DistMax", Dist_Max)
        counter+=2

        print(f"Mesh succesfully refined around faults.")

    elif element_type.lower() == "surface":
        surface_id = []
        
        for val in element_list:
            surface_tag = 1 + 10 * (val - 1)
            surface_id.append(surface_tag)
        

        gmsh.model.mesh.field.add("Distance", counter + 1)
        gmsh.model.mesh.field.setNumbers(counter + 1, "SurfacesList", surface_id)
        gmsh.model.mesh.field.setNumber(counter + 1, "Sampling", sampling)

        gmsh.model.mesh.field.add("Threshold", counter + 2)
        gmsh.model.mesh.field.setNumber(counter + 2, "InField", counter + 1)
        gmsh.model.mesh.field.setNumber(counter + 2, "SizeMin", Size_Min)
        gmsh.model.mesh.field.setNumber(counter + 2, "SizeMax", Size_Max)
        gmsh.model.mesh.field.setNumber(counter + 2, "DistMin", Dist_Min)
        gmsh.model.mesh.field.setNumber(counter + 2, "DistMax", Dist_Max)
        counter+=2

        print(f"Mesh succesfully refined around surface {element_list}.")

    gmsh.model.geo.synchronize()

def physical_group(
    element_type: str,
    element_list: list
):
    """
    Geo2Gmsh offers an adapted version of Gmsh’s physical group assignment to automatically assign physical group IDs to wells, faults, surfaces or vaolumes.

    Parameters:
        element_type (str): Specify the type of element. Choose from “well”, “fault”. “Surface” or “volume”.  
        element_list (list): For wells and faults, provide the list returned by the corresponding function. For surfaces and volumes, specify the surface ID                                (surf_id) or volume ID enclosed in square brackets (e.g., [1]). 
    """
    gmsh.model.geo.synchronize()
    
    global counter1
    if element_type.lower() == "well":

        physical_curve_id = gmsh.model.addPhysicalGroup(1, element_list, counter1+1)
        print("physical_id_(wells):", physical_curve_id)
        counter1+=1

    elif element_type.lower() == "fault":

        physical_fault_id = gmsh.model.addPhysicalGroup(2, element_list, counter1+1)
        print("physical_id_(faults):", physical_fault_id)
        counter1+=1

    elif element_type.lower() == "surface":
        surface_id = []
        for val in element_list:
            surface_tag = 1 + 10 * (val - 1)
            surface_id.append(surface_tag)
        physical_surface_id = gmsh.model.addPhysicalGroup(2, surface_id, counter1+1)
        print("physical_id_(surfaces):", physical_surface_id)
        counter1+=1
        
    elif element_type.lower() == "volume":
        for val in element_list:
            physical_volume_id = gmsh.model.addPhysicalGroup(3, [val], counter1+1)
            print(f"physical_id_(volume {val}):", physical_volume_id)
            counter1+=1
        
    gmsh.model.geo.synchronize()



