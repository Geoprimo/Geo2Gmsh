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
    Create a surface mesh in Gmsh from a text file with coordinates.

    Parameters:
        surf_id (int): Surface identifier, used for tagging.
        file_name (str): Path to text file with x,y,z data.
        samx (int): Number of samples in x-direction.
        samy (int): Number of samples in y-direction.
        v_ex (float): Vertical exaggeration factor.
        base_level (float): Base elevation level to offset z-coordinates.
        show_color (bool): Whether to create a colored view of elevation.
    """
    Nx = samx - 1
    Ny = samy -1
    
    # Cargar datos y extraer dimensiones
    data = np.loadtxt(f"layers/{file_name}", skiprows=1)
    id = (surf_id-1) * 10


    #Defining extend from txt file
    min_x = np.min(data[:, 1])
    max_x = np.max(data[:, 1])
    min_y = np.min(data[:, 2])
    max_y = np.max(data[:, 2])

    
    z = data[:, 3] * v_ex
    Z = np.reshape(z, [samy, samx])

    # Creacion de arrays en ejes x e y
    x = np.linspace(min_x, max_x, num=samx, endpoint=True)   
    y = np.linspace(min_y, max_y, num=samy, endpoint=True)

    #Creacion de puntos en malla
    X, Y = np.meshgrid(x, y, indexing='xy')
    Y = Y[::-1, :]
    
    # Helper function to return a node tag given three indices i, j and k:
    def tag(i, j):
        return ((Ny + 1) * i + j + 1) 
        
    # coordenadas de todos los nodos
    coords = []
    # tag de los nodos
    nodes = []
    # conectividad de elementos triangualares  (3 tags de nodo por triangulo):
    tris = []
    #conectividad de los elementos linea en los 4 limites
    lin = [[], [], [], []]
    
    #conectividad de los elementos punto en las cuatro esquinas
    pnt = [tag(0, Ny), tag(Nx, Ny), tag(Nx, 0), tag(0, 0)]

    #Llenado de los vectores creados anteriormente:
    for i in range( Nx + 1):    #columnas
        for j in range( Ny, -1, -1 ):      #filas
            nodes.append(tag(i, j))
            coords.extend([X[j, i], Y[j, i], Z[j, i]])
            if i > 0 and j > 0:
                tris.extend([tag(i - 1, j - 1), tag(i, j - 1), tag(i - 1, j)])
                tris.extend([tag(i, j - 1), tag(i, j), tag(i - 1, j)])
            if (i == 0 or i == Nx) and j > 0:
                lin[3 if i == 0 else 1].extend([tag(i, j - 1), tag(i, j)])
            if (j == 0 or j == Ny) and i > 0:
                lin[0 if j == 0 else 2].extend([tag(i - 1, j), tag(i, j)])
                
    # Crear 4 puntos discretos para las cuadro esquinas del terreno:
    for i in range(4):
        gmsh.model.addDiscreteEntity(0, i + 1 + id)
        
    # Modifico las coordenadas de los puntos creadores anteriormente
    gmsh.model.setCoordinates(1 + id, min_x, min_y, coords[3 * tag(0, 0) - 1])
    gmsh.model.setCoordinates(2 + id, max_x, min_y, coords[3 * tag(Nx, 0) - 1])
    gmsh.model.setCoordinates(3 + id, max_x, max_y, coords[3 * tag(Nx, Ny) - 1])
    gmsh.model.setCoordinates(4 + id, min_x, max_y, coords[3 * tag(0, Ny) - 1])

    # Crear 4 curvas de borde del terreno con los puntos anteriores:
    for i in range(4):
        gmsh.model.addDiscreteEntity(1, i + 1 + id, [i + 1 + id, i + 2 + id if i < 3 else 1 + id])

    # Crear una superficie con las curvas creadas anteriormente:
    gmsh.model.addDiscreteEntity(2, 1 + id, [1 + id, 2 + id, 3 + id, 4 + id])

    # Agregar todos los nodos en la superficie (for simplicity... see below):
    gmsh.model.mesh.addNodes(2, 1 + id, nodes, coords)

    for i in range(4):
        #elementos punto
        gmsh.model.mesh.addElementsByType(i + 1 + id, 15, [], [pnt[i]])
        #elementos linea 2 nodos
        gmsh.model.mesh.addElementsByType(i + 1 + id, 1, [], lin[i])
    #elementos triangulo 3 nodos
    gmsh.model.mesh.addElementsByType(1 + id, 2, [], tris)


    
    # Crear una vista coloreada por el valor de Z (elevación)

    if show_color == True:
        view_tag = gmsh.view.add(f"Z_surface_{surf_id}")
        z_only = [[z] for z in coords[2::3]]  # Extraer cada tercer valor de coords = Z
        gmsh.view.addModelData(
            view_tag,
            0,                          # paso de tiempo (0 si no hay simulación)
            gmsh.model.getCurrent(),   # nombre del modelo actual
            "NodeData",                # tipo de datos
            nodes,                     # tags de nodo
            z_only                     # valores escalares por nodo (Z)
    )

    print(f"Superficie {surf_id } successfully created.")

def volume_generation(
    num_loaded_surfaces: int = 2
):
    """
    Generates 3D volumes from a stack of loaded geological surfaces.

    Parameters:
    -----------
    num_loaded_surfaces : int
        Total number of 2D surfaces loaded into the model.
        Must be >= 2. One volume is created between each pair of consecutive surfaces.
    """
    # Reclasificar los nodos en las curvas y puntos (since we put them all on the surface before with `addNodes' for simplicity)
    gmsh.model.mesh.reclassifyNodes()

    # Crear una geometria para las curvas discretas y superficie, de forma que se puede remallarlas
    # luego:
    gmsh.model.mesh.createGeometry()

    vertical_union = {}
    surfaces = {}
    volume = []

    for j in range (num_loaded_surfaces - 1):  #cantidad de superficies menos 1
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
    txt contains points along the well. 
    txt struc:

    |   id    |    x    |    y    |    z    |
     
    n rows x 4 columns!
    
    Parameters:
        file_name (str): Path to text file with x,y,z data.
        v_ex (float): Vertical exaggeration factor.
        surf_id: ID of the top surface (e.g., topography). This is the surface that contains the top of the well.  
        well_id: ID of the well user-defined.
    """
    #Archivo de texto con x, y, z_top y z_base. Estos puntos pueden ser extraidos de Qgis por muestreo a lo largo de la falla en cada superfice
    well_txt = np.loadtxt(f"wells/{file_name}", skiprows=1)

    x = well_txt[:, 1]  
    y = well_txt[:, 2]
    z = well_txt[:, 3] * v_ex
    
 # #Agrega punto del tope
    well_points  = []
    for xi, yi, zi in zip(x, y, z):
        tag = gmsh.model.geo.addPoint(xi,yi,zi)
        well_points.append(tag)
    
    well_line  = []
    for i in range(len(well_points) - 1):
        tag = gmsh.model.geo.addLine(well_points[i], well_points[i+1])
        well_line.append(tag)

    gmsh.model.geo.synchronize()
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
    Create a fault that intercepts surface only.

    Parameters:
        file_name (str): Path to text file with x,y,z data.
        v_ex (float): Vertical exaggeration factor.
        lc: (float): longitud caracteristica
    """
    #Archivo de texto con x, y, z_top y z_base. Estos puntos pueden ser extraidos de Qgis por muestreo a lo largo de la falla en cada superfice
    fault_txt = np.loadtxt(f"faults/{file_name}", skiprows=1)

    x = fault_txt[:, 1]  
    y = fault_txt[:, 2]
    z_top = fault_txt[:, 3]*v_ex

    #Agrega los puntos del tope
    top_points  = []
    for xi, yi, zi in zip(x, y, z_top):
        tag = gmsh.model.geo.addPoint(xi,yi,zi)
        top_points.append(tag)

    #Agrega los puntos de la base
    # base_points  = []
    # for xi, yi, zi in zip(x, y, z_base):
    #     tag = gmsh.model.geo.addPoint(xi,yi,zi)
    #     base_points.append(tag)

    line_top = []
    #line_base = []

    #creando una linea con puntos del tope y otra con puntos de base
    for i in range(len(top_points) - 1):
        line_tag = gmsh.model.geo.addLine(top_points[i], top_points[i+1])
        line_top.append(line_tag)
    
        # line_tag = gmsh.model.geo.addLine(base_points[i], base_points[i+1])
        # line_base.append(line_tag)

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
        
    
    # Extrusión lineal:
    extruded = gmsh.model.geo.extrude(to_extrude, dx, dy, dz)
    sect_fault = [e[1] for e in extruded if e[0] == 2]

    # vertical_line = []
    # #creando lineas verticales
    # for i in range(len(top_points)):
    #     tag = gmsh.model.geo.addLine(top_points[i], base_points[i])
    #     vertical_line.append(tag)
    

    # sect_fault = [] 
    # for i in range(len(top_points) - 1):   
    #     loop = gmsh.model.geo.addCurveLoop([line_top[i], vertical_line[i+1], -line_base[i], -vertical_line[i]])
    #     tag = gmsh.model.geo.addPlaneSurface([loop])
    #     sect_fault.append(tag)
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



