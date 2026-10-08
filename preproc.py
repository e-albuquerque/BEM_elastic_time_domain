import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import meshio

def input_data():
    """
    This function defines input parameters for the Boundary Element Method (BEM)
    simulation applied to the plane elasticity equation with radial integration for domain
    integrals.

    Returns:
        dict: A dictionary containing the following key-value pairs:
            * bound_cond (dict): Boundary conditions for each edge:
                *
            * E (float): Material elastic modulus
            * nu (float): Poisson's ratio
            * file_name (str): Name of the .msh file (Gmsh mesh)
            * qpoint (str): String representing the coordinates of the heat source point
    """
    # Directly define values here for scripting (modify as needed)
    # Boundary conditions are specified in local coordinate system (normal-tangent)
    # type_x = type of the boundary condition in x direction
    # type_y = type of the boundary condition in y direction
    #    (type_x and type_y: 0 = displacement is know while 1 = traction is known)

    bound_cond = {'botton': {'type_x': 1, 'value_x': 0., 'tfun': 1, 'type_y': 0, 'value_y': 0., 'tfun': 1},
                  'right': {'type_x': 0, 'value_x': 0., 'tfun': 1, 'type_y': 1, 'value_y': 0., 'tfun': 1},
        'top': {'type_x': 1, 'value_x': 0., 'tfun': 1, 'type_y': 1, 'value_y': 25., 'tfun': 1},
        'left': {'type_x': 0, 'value_x': 0., 'tfun': 1, 'type_y': 1, 'value_y': 0., 'tfun': 1}
    }

    E = 260.0e8 # # Elastic modulus
    nu = 0.3 # Poisson ration
    file_name = 'placa'
    rho = 5000. # density
    tstep = 2.0E-4 # time step
    ntfunctions = 1 # number of time functions

    return {'bound_cond': bound_cond, 'E': E, 'nu': nu, 'file_name': file_name,\
            'rho' : rho, 'tstep' : tstep, 'ntfunctions' : ntfunctions}


def compute_inodes(file_name, bound_cond):
    """
    This function reads a mesh file and extracts information about nodes, elements,
    and boundary conditions.

    Args:
      file_name: Name of the mesh file.
      bound_cond: Dictionary containing boundary conditions for each edge.

    Returns:
      A dictionary containing information about nodes, elements, and boundary
      conditions.
    """
    mesh = meshio.read(file_name + '.msh')
    coordinates = mesh.points
    line_elements = mesh.cells_dict['line3']
    segments = mesh.cell_data_dict['gmsh:physical']['line3']

    bc_info = {}
    for key, value in bound_cond.items():
        segment = mesh.field_data[key][0]
        bc_info[key] = {
            'type_x': value['type_x'],
            'value_x': value['value_x'],
            'segment': segment,
            'type_y': value['type_y'],
            'value_y': value['value_y']
        }

    boundary_nodes = np.unique(line_elements)

    return {
        'nodes': boundary_nodes,
        'coordinates': coordinates,
        'elements': line_elements,
        'segments': segments,
        'bc_info': bc_info
    }

def comp_xy_plots(elem, nodes):
    """
    Calculates the nodes and outward normal vector for each boundary element.

    This function takes element connectivity information (`elem`) and node coordinates (`nodes`)
    to determine the geometric center of each boundary element and the unit normal vector
    pointing outward from that element.

    Args:
        elem (np.ndarray): A 2D array defining the connectivity of the boundary elements.
                          Each row contains the indices of the two nodes that form an element.
                          Shape: (num_elements, 2).
        nodes (np.ndarray): A 2D array containing the coordinates (x, y) of each node in the mesh.
                           Shape: (num_nodes, 2).

    Returns:
        node_des (np.ndarray): A 2D array containing the coordinates (x, y) of the discontinuous nodes
                              of each boundary element. Shape: (2*num_elements, 2).
        normal (np.ndarray): A 2D array containing the unit outward normal vectors for each
                             boundary element. Shape: (num_elements, 2).
    """
    npp = 8                              # Number of plot points
    eta = np.linspace(-1., 1., npp)        # Plot points

    num_elements = elem.shape[0]

    xplot = np.zeros((num_elements,npp))        # Store x-coordinates for plotting
    yplot = np.zeros((num_elements,npp))        # Store y-coordinates for plotting

    # Iterate over each boundary element
    for i in range(num_elements):
        # Get node indices for the current element
        node1, node3, node2 = elem[i]

        # Get coordinates of the two nodes forming the element
        x1, y1 = nodes[node1][0:2]
        x2, y2 = nodes[node2][0:2]
        x3, y3 = nodes[node3][0:2]


        # Iterate over Gauss points for numerical integration
        for k in range(npp):
            # Shape functions and their derivatives
            N1 = 0.5 * eta[k] * (eta[k] - 1)
            N2 = 1 - eta[k]**2
            N3 = 0.5 * eta[k] * (eta[k] + 1)

            # Calculate coordinates and normal vector at the Gauss point
            x, y = N1 * x1 + N2 * x2 + N3 * x3, N1 * y1 + N2 * y2 + N3 * y3
            ivet = npp * i + k  # Index for storing plot coordinates
            xplot[i,k], yplot[i,k] = x, y

    return xplot,yplot

# Read information about input data
inp_data = input_data()

# Format input data
computed_data = compute_inodes(inp_data['file_name'], inp_data['bound_cond'])

xplot,yplot = comp_xy_plots(computed_data['elements'], computed_data['coordinates'])



def show_geometry(nodes,xplot,yplot,elem):
    xnode=nodes[:,0]
    ynode=nodes[:,1]
    ax = plt.axes()
    plt.plot(xnode[elem[:,0]],ynode[elem[:,0]],"gd",markersize=8)	# Plot the node of the elements
    plt.plot(xnode,ynode,"ro",markersize=2)	# Plot the node of the elements
    plt.axis("equal")

    n_el = elem.shape[0]
    for i in range(n_el):
        x=xplot[i,:]
        y=yplot[i,:]
        line_segments = LineCollection([list(zip(x, y))],
                                        linestyles='solid')
        ax.add_collection(line_segments)

    plt.savefig("geometry.png")

# Re-get the global `inp_data` dictionary to avoid shadowing
global_input_data = input_data()

inp_data_filename='str12.dat'

# Time function

time_function = np.array([[0., 1.],[1.E10, 1.]]) # Heaviside (step) function

# plot time function
plt.plot(time_function[:,0],time_function[:,1])
plt.xlabel('t')
plt.ylabel('f(t)')
plt.show()



def get_midnode_edge(computed_data, edge):
    # Get the segment ID for the specified edge
    segment_id = computed_data['bc_info'][edge]['segment']

    # Filter elements that belong to this segment
    mask = computed_data['segments'] == segment_id

    # Get global indices of the elements for this edge
    global_indices = np.where(mask)[0]
    n_elements = len(global_indices)

    if n_elements % 2 != 0:
        # Ímpar: o nó central é o nó do meio do elemento central
        local_mid_idx = n_elements // 2
        global_mid_idx = global_indices[local_mid_idx]
        mid_element = computed_data['elements'][global_mid_idx]

        # Para Gmsh line3 [start, end, mid], o nó do meio é o índice 2
        mid_node = mid_element[2]
    else:
        # Par: o nó central da aresta é o nó compartilhado pelos dois elementos centrais.
        # Pegamos o elemento da esquerda (onde o nó é o nó final/terceiro nó).
        local_mid_idx = (n_elements // 2) - 1
        global_mid_idx = global_indices[local_mid_idx]
        mid_element = computed_data['elements'][global_mid_idx]

        # O nó final do elemento (terceiro nó) fica no índice 1
        mid_node = mid_element[1]

    return mid_node, global_mid_idx

edge = 'top'
mid_node_number, mid_element_idx = get_midnode_edge(computed_data, edge)
print(f"Mid node number for '{edge}' edge: {mid_node_number}")
print(f"Mid element number for '{edge}' edge: {mid_element_idx}")

def create_input_file(filename, computed_data, input_params,time_function):
    """
    Creates an input file for the Fortran BEM program (quadplet) in the time domain
    based on computed mesh data and boundary conditions.
    """
    with open(filename, 'w') as f:
        # 1. Title Line
        f.write("INFINITE STRIP UNDER UNIFORM TRACTION\n\n")

        # Extract data from computed_data and input_params
        nodes = computed_data['coordinates']
        elements = computed_data['elements']
        num_elements = elements.shape[0]

        # 2. Basic Parameter Line
        ntfun = input_params.get('ntfunctions', 1)
        ipr = 1 # Printing interval
        nstep = 20 # Number of time steps
        tstep = float(input_params.get('tstep', 2.0E-4))

        # G = E / 2(1+v)
        G = float(input_params['E'] / (2 * (1 + input_params['nu'])))
        rho = float(input_params['rho'])
        nu = float(input_params['nu'])

        # Format float numbers to use uppercase E for Fortran compatibility
        f.write(f" {num_elements},{ntfun}   {ipr},{nstep},{tstep:.4E}   {G:.4E},{nu:.4E},{rho:.4E}\n\n")

        # 3. Boundary Nodes Coordinates Lines
        # Sequência geométrica desejada: no1 e no_meio de cada elemento
        node_coords_list = []
        for elem in elements:
            no1_idx = elem[0]
            no_meio_idx = elem[2]
            x1, y1 = nodes[no1_idx][0], nodes[no1_idx][1]
            xm, ym = nodes[no_meio_idx][0], nodes[no_meio_idx][1]
            node_coords_list.append(f"{float(x1):.4E} {float(y1):.4E}")
            node_coords_list.append(f"{float(xm):.4E} {float(ym):.4E}")

        f.write(" ")
        for i in range(0, len(node_coords_list), 8): # Wrap at 8 nodes per line
            f.write("  ".join(node_coords_list[i:i+8]) + "\n ")
        f.write("\n")

        # 4. Time Function lines
        #  scale 1.0
        scale = 1.
        n_time_points = time_function.shape[0]
        f.write(f" {n_time_points},{scale}\n ")
        
        # Write all points on a single line
        time_points_str = []
        for i in range(n_time_points):
            time_points_str.append(f"{time_function[i,0]},{time_function[i,1]}")
        f.write("  , ".join(time_points_str) + "\n\n")

        # 5. Boundary Conditions Lines
        segments = computed_data['segments']
        bc_info = computed_data['bc_info']

        for i in range(num_elements):
            segment_index = segments[i]
            # Find the boundary condition data corresponding to the current segment
            for bc_name, bc_data in bc_info.items():
                if bc_data['segment'] == segment_index:
                    tp_x = bc_data['type_x']
                    # Gracefully handle both lists and scalar floats from input_data
                    vl_x = float(bc_data['value_x'][0] if isinstance(bc_data['value_x'], list) else bc_data['value_x'])
                    tf_x = bc_data.get('tfun', 1)

                    tp_y = bc_data['type_y']
                    vl_y = float(bc_data['value_y'][0] if isinstance(bc_data['value_y'], list) else bc_data['value_y'])
                    tf_y = bc_data.get('tfun', 1)
                    break

            line_parts = []
            for n in range(3):
                line_parts.append(f"{tp_x},{vl_x:.4E},{tf_x} {tp_y},{vl_y:.4E},{tf_y}")
            f.write(f" {'  '.join(line_parts)}\n")

# Call the function to generate str12.dat
create_input_file(inp_data_filename, computed_data, global_input_data, time_function)


