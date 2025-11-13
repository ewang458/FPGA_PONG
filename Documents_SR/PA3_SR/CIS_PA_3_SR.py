import numpy as np
import os

# ======================== CONFIGURATION ==========================
# Path to the dataset folder (contains all PA3 debug and unknown files)
data_path = "2025 PA345 Student Data"
output_path = "OUTPUT"


# ===================== 3D GEOMETRY UTILITIES =====================
def pos_3D(x, y, z):
    """Return a 3D vector as NumPy array."""
    return np.array([x, y, z])


def rotate(R, P):
    """Apply rotation R to point/vector P."""
    P = np.asarray(P).reshape(3,)
    return R @ P


def translate(P, d):
    """Apply translation d to point P."""
    P, d = np.asarray(P).reshape(3,), np.asarray(d).reshape(3,)
    return P + d


def frame_transformation(R, Pin, d):
    """Perform full rigid-body transform: output = R * Pin + d."""
    Pin, d = np.asarray(Pin).reshape(3,), np.asarray(d).reshape(3,)
    R = np.asarray(R).reshape(3, 3)
    return R @ Pin + d


def inverse_transformation(R, p):
    """Compute inverse transformation: R_inv, p_inv st inv(R, p) * x = R_inv * x + p_inv"""
    R_inv = R.T
    p_inv = -R_inv @ p
    return R_inv, p_inv


# ==================== 3D POINT REGISTRATION ======================
def point2point_3Dregistration(A, B):
    """
    Compute rotation (R) and translation (p) aligning A -> B using SVD.
    Solves min_{R,p} ||B - (R A + p)||^2.
    """
    A, B = np.asarray(A, float), np.asarray(B, float)
    if A.shape != B.shape or A.shape[1] != 3:
        raise ValueError("A and B must have same shape (N,3)")

    # Compute centroids of both sets
    cA, cB = A.mean(axis=0), B.mean(axis=0)
    A0, B0 = A - cA, B - cB

    # Covariance matrix and SVD decomposition
    H = A0.T @ B0
    U, S, Vt = np.linalg.svd(H)

    # Compute rotation matrix
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:  # Ensure right-handed rotation
        Vt[2, :] *= -1
        R = Vt.T @ U.T

    # Translation vector
    p = cB - R @ cA
    return R, p


# MAIN FUNCTIONS
def find_closest_point_on_triangle(point, p, q, r):
    """
    Find the closest point on triangle (p, q, r) to the given point.
    |_ use least square approach
    |_ return the closest point on the triangle to the given point
    """
    point = np.asarray(point).reshape(3,)
    p = np.asarray(p).reshape(3,)
    q = np.asarray(q).reshape(3,)
    r = np.asarray(r).reshape(3,)
    
    #express point in terms of triangle basis: point - p = lambda*(q-p) + mu*(r-p)
    right1 = q - p  # edge vector 1
    right2 = r - p  # edge vector 2
    left = point - p  # vector from p to point
    
    # solve via least square
    A_matrix = np.array([right1, right2]).T  # Shape (3, 2)
    coeffs = np.linalg.lstsq(A_matrix, left, rcond=None)[0]
    lamda, mu = coeffs[0], coeffs[1]
    
    # get initial estimate
    c = p + lamda * right1 + mu * right2
    
    #BOUNDARY CASES
    # if the point is outside triangle, we need to project onto edges or vertices
    if lamda < 0 and lamda + mu > 1:
        # closest to vertex r
        c = r
    elif lamda + mu > 1 and mu < 0:
        # closest to vertex q
        c = q
    elif lamda < 0 and mu < 0:
        #closest to vertex p
        c = p
    elif lamda < 0:
        # project onto edge r-p
        c = project_on_segment(c, r, p)
    elif mu < 0:
        # project onto edge p-q
        c = project_on_segment(c, p, q)
    elif lamda + mu > 1:
        # project onto edge q-r
        c = project_on_segment(c, q, r)
    
    return c


def project_on_segment(point, p, q):
    """
    Project a point onto the line segment between p and q.
    """
    point = np.asarray(point).reshape(3,)
    p = np.asarray(p).reshape(3,)
    q = np.asarray(q).reshape(3,)
    
    #calculate t -> proj = p + t*(q-p)
    # Constrain to [0, 1]
    pq = q - p
    pp = point - p
    
    # Project pp onto pq
    t = np.dot(pp, pq) / np.dot(pq, pq)
    t = max(0, min(1, t)) 
    
    return p + t * pq


# ==================== CLOSEST POINT ON MESH ======================
def find_closest_point_on_mesh_linear(point, vertices, indices):
    """
    Find the closest point on a mesh to a given point using linear search.
    |_ iterate through all triangles to find the min distance.
    |_ return the closest point on the mesh to the given point
    """
    point = np.asarray(point).reshape(3,)
    vertices = np.asarray(vertices)
    indices = np.asarray(indices, dtype=int)
    
    # initialize first triangle
    triangle = get_triangle_coordinates(0, vertices, indices)
    closest_point = find_closest_point_on_triangle(point, triangle[0], triangle[1], triangle[2])
    min_dist = np.linalg.norm(point - closest_point)
    
    #iterate through all triangles
    for i in range(1, len(indices)):
        triangle = get_triangle_coordinates(i, vertices, indices)
        curr_point = find_closest_point_on_triangle(point, triangle[0], triangle[1], triangle[2])
        curr_dist = np.linalg.norm(point - curr_point)
        
        if curr_dist < min_dist:
            min_dist = curr_dist
            closest_point = curr_point
    
    return closest_point


def get_triangle_coordinates(triangle_idx, vertices, indices):
    """
    Get the coordinates of the three vertices of a triangle.
    """
    idx = indices[triangle_idx]
    return vertices[idx[0]], vertices[idx[1]], vertices[idx[2]]


#FILE PARSER
def parse_body_file(path):
    """
    Parse body definition file
    """
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    
    # line 1 = N_markers, filename
    first_line = lines[0].split()
    N_markers = int(first_line[0])
    
    # get marker coordinates
    markers = []
    for i in range(1, 1 + N_markers):
        coords = [float(x) for x in lines[i].split()]
        markers.append(coords)
    
    # tip coordinates in last line
    tip = [float(x) for x in lines[-1].split()]
    
    return N_markers, np.array(markers), np.array(tip)


def parse_mesh_file(path):
    """
    Parse surface mesh file, return the vertex coordinates and the triangle vertex indices
    """
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    
    #  line 1 = N_vertices
    N_vertices = int(lines[0])
    
    # read vertex coord.
    vertices = []
    for i in range(1, 1 + N_vertices):
        coords = [float(x) for x in lines[i].split()]
        vertices.append(coords)
    
    #N_triangles
    N_triangles = int(lines[1 + N_vertices])
    
    indices = []
    for i in range(1 + N_vertices + 1, 1 + N_vertices + 1 + N_triangles):
        idx_line = lines[i].split()
        i1, i2, i3 = int(idx_line[0]), int(idx_line[1]), int(idx_line[2])
        indices.append([i1, i2, i3])
    
    return np.array(vertices), np.array(indices)


def parse_sample_readings(path, N_A_markers, N_B_markers):
    """
    Parse sample readings file, get number of sample grames, A and B marker coordinates in tracker frame
    """
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    
    # line 1 = N_readings, N_samples, filename
    first_line_parts = lines[0].split(",")
    N_readings = int(first_line_parts[0].strip())
    N_samples = int(first_line_parts[1].strip())
    
    A_readings = []
    B_readings = []
    
    # parse each sample frame
    idx = 1
    for sample in range(N_samples):
        # Read A markers
        A_frame = []
        for i in range(N_A_markers):
            coords = [float(x.strip()) for x in lines[idx].split(",")]
            A_frame.append(coords)
            idx += 1
        
        # Read B markers
        B_frame = []
        for i in range(N_B_markers):
            coords = [float(x.strip()) for x in lines[idx].split(",")]
            B_frame.append(coords)
            idx += 1
        
        # Skip dummy markers
        idx += (N_readings - N_A_markers - N_B_markers)
        
        A_readings.append(np.array(A_frame))
        B_readings.append(np.array(B_frame))
    
    return N_samples, A_readings, B_readings


# FRAME PROCESSING
def process_frame(A_markers_body, B_markers_body, A_tip_body, A_readings_frame, B_readings_frame, 
                  vertices, indices, F_reg=None):
    """
    Process a single frame to compute pointer tip position and closest mesh point.
    
    For Problem 3:
    |_ compute FA,k: transformation from A markers (body frame) to tracker frame
    |_ compute FB,k: transformation from B markers (body frame) to tracker frame
    |_ compute dk = FB,k^-1 * FA,k * A_tip (pointer tip in B body frame)
    |_ compute sk = F_reg * dk (pointer tip in CT frame, for Problem 3, F_reg = I)
    |_ find closest point ck on mesh to sk

    Return pointer tip position in B body frame (d) and closest point on mesh to transformed tip position (c)
    """
    # FA,k: A_markers_body -> A_readings_frame
    R_A, p_A = point2point_3Dregistration(A_markers_body, A_readings_frame)
    
    # FT2B: B_readings_frame -> B_markers_body
    R_T2B, p_T2B = point2point_3Dregistration(B_readings_frame, B_markers_body)
    
    # Transform A_tip to tracker frame
    A_tip_tracker = frame_transformation(R_A, A_tip_body, p_A)
    
    # Transform A_tip_tracker from tracker frame to B body frame
    d = frame_transformation(R_T2B, A_tip_tracker, p_T2B)
    
    #registration transformation (F_reg = I, so sk = dk)
    if F_reg is None:
        # Identity transformation
        s = d
    else:
        R_reg, p_reg = F_reg
        s = frame_transformation(R_reg, d, p_reg)
    
    #find closest point on mesh to s
    c = find_closest_point_on_mesh_linear(s, vertices, indices)
    
    return d, c


# WRITE OUTPUT FILE
def write_output_file(output_path, problem_num, letter, N_samples, ds, cs):
    """
    FORMAT:
    |_ first line: N_samples, filename
    |_ for each sample: dx, dy, dz, cx, cy, cz, distance
    """
    os.makedirs(output_path, exist_ok=True)
    
    if letter.upper() <= 'F':
        filename = f"PA{problem_num}-{letter.upper()}-Debug-Output.txt"
    else:
        filename = f"PA{problem_num}-{letter.upper()}-Unknown-Output.txt"
    
    filepath = os.path.join(output_path, filename)
    
    with open(filepath, "w") as f:
        f.write(f"{N_samples}, {filename}\n")
        
        for i in range(N_samples):
            d = ds[i]
            c = cs[i]
            distance = np.linalg.norm(d - c)
            f.write(f"{d[0]:.3f}, {d[1]:.3f}, {d[2]:.3f}, {c[0]:.3f}, {c[1]:.3f}, {c[2]:.3f}, {distance:.3f}\n")
    
    return filepath


# DATA PROCESSING AND MAIN FUNCTION
def process_dataset(data_prefix, problem_num=3):
    """
    Process a given dataset
    """
    print(f"Processing {data_prefix}...")
    
    try:
        # 1. Parse body definition files
        body_A_path = os.path.join(data_path, f"Problem{problem_num}-BodyA.txt")
        body_B_path = os.path.join(data_path, f"Problem{problem_num}-BodyB.txt")
        
        if not os.path.exists(body_A_path):
            raise FileNotFoundError(f"Body A file not found: {body_A_path}")
        if not os.path.exists(body_B_path):
            raise FileNotFoundError(f"Body B file not found: {body_B_path}")
        
        N_A_markers, A_markers_body, A_tip_body = parse_body_file(body_A_path)
        N_B_markers, B_markers_body, B_tip_body = parse_body_file(body_B_path)
        
        print(f"  Loaded body A: {N_A_markers} markers")
        print(f"  Loaded body B: {N_B_markers} markers")
        
        # parse mesh file
        mesh_path = os.path.join(data_path, f"Problem{problem_num}MeshFile.sur")
        if not os.path.exists(mesh_path):
            mesh_path = os.path.join(data_path, f"Problem{problem_num}Mesh.sur")
        if not os.path.exists(mesh_path):
            raise FileNotFoundError(f"Mesh file not found: {mesh_path}")
        vertices, indices = parse_mesh_file(mesh_path)
        print(f"Loaded mesh: {len(vertices)} vertices, {len(indices)} triangles")
        
        # parse sample readings
        sample_readings_path = os.path.join(data_path, f"{data_prefix}-SampleReadingsTest.txt")
        if not os.path.exists(sample_readings_path):
            raise FileNotFoundError(f"Sample readings file not found: {sample_readings_path}")
        N_samples, A_readings, B_readings = parse_sample_readings(
            sample_readings_path, N_A_markers, N_B_markers
        )
        print(f"  Loaded {N_samples} sample frames")
        
        # process each frame
        # for PA3, F_reg = identity (so no registration needed)
        F_reg = None  # Identity
        
        ds = []
        cs = []
        
        for frame_idx in range(N_samples):
            d, c = process_frame(
                A_markers_body, B_markers_body, A_tip_body,
                A_readings[frame_idx], B_readings[frame_idx],
                vertices, indices, F_reg
            )
            ds.append(d)
            cs.append(c)
        
        ds = np.array(ds)
        cs = np.array(cs)
        
        # get letter / extract output
        letter = data_prefix.split("-")[1].lower()  
        output_file = write_output_file(output_path, problem_num, letter, N_samples, ds, cs)
        
        print(f"  Completed {data_prefix}")
        print(f"    Output: {output_file}")
        
    except Exception as e:
        print(f"  Error processing {data_prefix}: {e}")
        import traceback
        traceback.print_exc()


def main():
    """Main function to process all PA3 datasets."""
    problem_num = 3
    
    datasets = []
    for letter in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'J']:
        if letter <= 'F':
            datasets.append(f"PA{problem_num}-{letter}-Debug")
        else:
            datasets.append(f"PA{problem_num}-{letter}-Unknown")
    
    for data_prefix in datasets:
        sample_file = os.path.join(data_path, f"{data_prefix}-SampleReadingsTest.txt")
        if os.path.exists(sample_file):
            process_dataset(data_prefix, problem_num)
        else:
            print(f"Skipping {data_prefix} (file not found)")


if __name__ == "__main__":
    main()

