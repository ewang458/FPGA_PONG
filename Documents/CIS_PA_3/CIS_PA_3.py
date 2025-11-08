import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize
import os

# ======================== CONFIGURATION ==========================
# Path to the dataset folder (contains all PA1 debug and unknown files)
extract_path = r"C:\Users\aarus\Documents\CIS_PA_3\DATA\\"
datasets = os.listdir(extract_path)


# ===================== 3D GEOMETRY UTILITIES =====================
def pos_3D(x, y, z):
    """Return a 3D vector as NumPy array."""
    return np.array([x, y, z])

# --- Rotation matrices around X, Y, Z axes respectively ---
def Rot_X(theta):
    """Rotation matrix about X-axis by angle theta (radians)."""
    return np.array([
        [1, 0, 0],
        [0, np.cos(theta), -np.sin(theta)],
        [0, np.sin(theta),  np.cos(theta)]
    ])

def Rot_Y(theta):
    """Rotation matrix about Y-axis by angle theta (radians)."""
    return np.array([
        [ np.cos(theta), 0, np.sin(theta)],
        [ 0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])

def Rot_Z(theta):
    """Rotation matrix about Z-axis by angle theta (radians)."""
    return np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta),  np.cos(theta), 0],
        [0, 0, 1]
    ])

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


# ===================== PIVOT CALIBRATION =========================
def pivot_calibration(frames):
    """
    Compute stationary pivot point given multiple frames of probe markers.
    Each frame: G_k = R_k * g_local + t_k.
    Returns p_dimple (fixed pivot point).
    """
    if len(frames) < 2:
        raise ValueError("Need at least two frames.")

    ref = np.asarray(frames[0], float)
    G0 = ref.mean(axis=0)
    g_local = ref - G0  # Local marker coordinates

    A_blocks, b_list = [], []
    for F in frames:
        F = np.asarray(F, float)
        # Get rotation & translation of this frame
        Rk, tk = point2point_3Dregistration(g_local, F - G0)
        # Build least-squares blocks
        A_blocks.append(np.hstack([Rk, -np.eye(3)]))
        b_list.append(-(G0 + tk))

    # Stack all equations and solve A x = b
    A, b = np.vstack(A_blocks), np.concatenate(b_list)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    p_dimple = sol[3:]  # Extract pivot position (in global frame)
    return p_dimple


# ======================== FILE PARSERS ===========================
def parseBodyA(path):
    """Read BodyA file for A markers coordinates and A tip coordinates."""
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    Na= int(lines[0].split()[0])
    idx = 1
    A = [[float(x) for x in lines[idx+i].split()] for i in range (Na)]
    A_tip = [[float(x) for x in lines[idx + Na].split()]]
    return Na, np.array(A), np.array(A_tip)

def parseBodyB(path):
    """Read BodyB file for B markers coordinates and B tip coordinates."""
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    Nb= int(lines[0].split()[0])
    idx = 1
    B = [[float(x) for x in lines[idx+i].split()] for i in range (Nb)]
    B_tip = [[float(x) for x in lines[idx + Nb].split()]]
    return Nb, np.array(B), np.array(B_tip)

def parseMesh(path):
    """Parse Mesh file for vertices and indices"""
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    Nv= int(lines[0].split()[0])
    V = [[float(x) for x in lines[1 + i].split()] for i in range (Nv)]
    Nt = int(lines[1 + Nv].split()[0])
    Indices = [[int(x) for x in lines[Nv + i + 2].split()[:3]] for i in range (Nt)]
    return Nv, np.array(V), Nt, np.array(Indices)

def parseSampleReadings(path, Na, Nb):
    """Parse sample readings for measurements of A, B, and D markers in tracker coordinates"""
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    Nsum, Ns = map(int, lines[0].split(",")[:2])
    Nd = Nsum - Na - Nb
    idx = 1
    A_frames = []
    B_frames = []
    D_frames = []
    for _ in range(Ns):
        A = np.array([[float(x) for x in lines[idx + i].split(",")] for i in range(Na)], float); idx += Na
        B = np.array([[float(x) for x in lines[idx + i].split(",")] for i in range(Nb)], float); idx += Nb
        D = np.array([[float(x) for x in lines[idx + i].split(",")] for i in range(Nd)], float); idx += Nd
        A_frames.append(A), B_frames.append(B), D_frames.append(D)
    return A_frames, B_frames, D_frames, Nd, Ns

# ==================== CLOSEST POINT SEARCH  =====================

def compute_bounds(V, Indices):
    """
    V: (nv,3)
    Indices: (nt,3) integers indexing into V
    Returns: lower (nt,3), upper (nt,3), tris (nt,3,3)
    """
    V = np.asarray(V, dtype=float)
    Indices = np.asarray(Indices, dtype=int)
    tris = V[Indices]           # shape (nt, 3, 3)
    lower = np.min(tris, axis=1)
    upper = np.max(tris, axis=1)
    return lower, upper, tris



def project_on_segment(a, p, q):
    ap = a - p
    qp = q - p
    denom = np.dot(qp, qp)
    if denom == 0.0:
        c = p.copy()
    else: 
        lambd = np.dot(ap, qp) / np.dot(qp, qp)
        lam_seg = np.maximum(0, np.minimum(lambd, 1))  # Fixed function names
        c = p + lam_seg * qp
    return c, np.sum((a - c) ** 2)


def closest_point_on_triangle(a, p, q, r):
    """
    Barycentric-based closest point on triangle pqr to point a.
    Returns (closest_point (3,), squared_distance)
    """
    a = np.asarray(a, dtype=float)
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    r = np.asarray(r, dtype=float)

    pq = q - p
    pr = r - p
    pa = a - p

    # Dot-products for barycentric solve 
    A = np.dot(pq, pq)
    B = np.dot(pq, pr)
    C = np.dot(pr, pr)
    D = np.dot(pq, pa)
    E = np.dot(pr, pa)

    det = A * C - B * B
    # Degenerate triangle (collinear / numerical): fallback to nearest vertex
    if abs(det) < 1e-12:
        # fallback: check distances to vertices and edges
        c_p, d2_p = p, np.sum((a - p) ** 2)
        c_q, d2_q = q, np.sum((a - q) ** 2)
        c_r, d2_r = r, np.sum((a - r) ** 2)
        # also check edges robustly
        c1, d21 = project_on_segment(a, p, q)
        c2, d22 = project_on_segment(a, q, r)
        c3, d23 = project_on_segment(a, r, p)
        all_cs = [c_p, c_q, c_r, c1, c2, c3]
        all_d2 = [d2_p, d2_q, d2_r, d21, d22, d23]
        idx = int(np.argmin(all_d2))
        return all_cs[idx], all_d2[idx]

    lam = (C * D - B * E) / det
    mu  = (A * E - B * D) / det

    # If point projects inside the triangle
    if lam >= 0.0 and mu >= 0.0 and (lam + mu) <= 1.0:
        c = p + lam * pq + mu * pr
        return c, np.sum((a - c) ** 2)

    # Otherwise, the closest point is on one of the edges
    if (lam < 0):
        c, d = project_on_segment(a, r, p)
    elif (mu < 0):
        c, d = project_on_segment(a, p, q)
    else: 
        c, d = project_on_segment(a, q, r)
    return c, float(d)

'''
def search_with_boxes(a, lower, upper, tris):
    """
    Bounding-box accelerated search for closest point on mesh.
    a: (3,)
    lower, upper: (nt,3)
    tris: (nt,3,3)
    Returns: closest_point (3,), squared_distance
    """
    a = np.asarray(a, dtype=float)
    Nt = tris.shape[0]
    mag = float("inf")
    best_c = None

    # Precompute sqrt of best_d2 when needed; use squared comparisons to avoid sqrt calls:
    for i in range(Nt):
        vmin = lower[i]
        vmax = upper[i]
        # Compute squared distance from point to AABB quickly (vectorized scalar)
        # dist_box = max(0, vmin - a, a - vmax) componentwise
        diff_lower = vmin - a
        diff_upper = a - vmax
        dist_comp = np.maximum(0.0, np.maximum(diff_lower, diff_upper))
        d2_box = np.dot(dist_comp, dist_comp)
        if d2_box > best_d2:
            continue  # cannot improve

        p, q, r = tris[i]
        c, d2 = closest_point_on_triangle(a, p, q, r)
        if d2 < best_d2:
            best_d2 = d2
            best_c = c

    return best_c, best_d2
'''
def search_with_boxes(a, lower, upper, tris):
    a = np.asarray(a)
    Nt = tris.shape[0]
    bound2 = np.inf
    closest = None

    for i in range(Nt):
        if np.any((a < (lower[i] - np.sqrt(bound2))) | (a > (upper[i] + np.sqrt(bound2)))):
            continue

        p, q, r = tris[i]
        c, d2 = closest_point_on_triangle(a, p, q, r)
        if d2 < bound2:
            bound2 = d2
            closest = c

    return closest, bound2

def write_output_files(filename, d, c, Ns):
    """Write output file 1 with em and optical dimples and C expected,
        and write output file 2 with CT tip positions in required format for each dataset."""
    #Write Output File 1
    with open(filename, "w") as f:
        f.write(f"{Ns} {filename}\n")
        for s in range(Ns):
            magdiff = np.sum((d[s] - c[s])**2)
            diff = magdiff**(0.5)
            f.write(f"{d[s][0]:.2f}\t{d[s][1]:.2f}\t{d[s][2]:.2f}\t\t{c[s][0]:.2f}\t{c[s][1]:.2f}\t{c[s][2]:.2f}\t{diff:.3f}\n")
        




# ======================= DATASET PIPELINE ========================
def process_dataset(data_prefix):
    """Run full calibration pipeline for a single dataset prefix."""
    print(f"Processing {data_prefix}...")

    try:
        # --- 1. Parse input files ---
        Na, A, A_tip = parseBodyA(os.path.join(extract_path,"Problem3-BodyA.txt"))
        Nb, B, B_tip = parseBodyB(os.path.join(extract_path,"Problem3-BodyB.txt"))
        Nv , V, Nt, Indices = parseMesh(os.path.join(extract_path,"Problem3MeshFile.sur"))
        A_frames, B_frames, D_frames, Nd, Ns = parseSampleReadings(os.path.join(extract_path, f"{data_prefix}-SampleReadingsTest.txt"), Na, Nb)
    except Exception as e:
        print(f"  Error processing {data_prefix}: {e}")

    lower, upper, tris = compute_bounds(V, Indices)
    c = []
    d = [] 

    A_tip = np.asarray(A_tip).reshape(3,)

    for k in range(Ns): 
        A_meas = np.asarray(A_frames[k])
        B_meas = np.asarray(B_frames[k])

        R_ak, p_ak = point2point_3Dregistration(A, A_meas)
        R_bk, p_bk = point2point_3Dregistration(B, B_meas)
        a_tracker = frame_transformation(R_ak, A_tip, p_ak)
        R_bk_inv = R_bk.T
        p_bk_inv = -R_bk.T@p_bk
        d_k = frame_transformation(R_bk_inv , a_tracker, p_bk_inv)
        c_k, _ = search_with_boxes(d_k, lower, upper, tris)

        c.append(c_k)
        d.append(d_k)

    out_name = f"{data_prefix}-myOutput.txt"
    write_output_files(out_name, d , c, Ns)





def main():
    """Main loop — process all datasets in directory."""
    # Ask the user for the dataset path
    '''
    extract_path = input("Please enter the path to the dataset folder: ").strip()  # Get user input
    
    # Optionally, ensure the path is valid
    while not os.path.isdir(extract_path):
        print("The path you entered is invalid. Please try again.")
        extract_path = input("Please enter the path to the dataset folder: ").strip()

    datasets = os.listdir(extract_path)
    '''

    for fname in datasets:
        if fname.endswith("-SampleReadingsTest.txt"):
            data_prefix = fname[:-len("-SampleReadingsTest.txt")]
            process_dataset(data_prefix)
        else:
            print(f"Skipping {fname}")


if __name__ == "__main__":
    main()