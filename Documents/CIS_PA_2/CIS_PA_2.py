import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize
from scipy.special import comb
import os

# ======================== CONFIGURATION ==========================
# Path to the dataset folder (contains all PA2 debug and unknown files)
extract_path = r"C:\Users\aarus\Documents\CIS_PA_2\DATA\\"
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
    p_tip = sol[:3] 
    return p_tip, g_local

# ===================== BERNSTEIN DEWARPING ====================+==
def bernstein_polynomial(n, k, v):
    """Compute Bernstein polynomial"""
    return comb(n, k) * (v** k) * ((1 - v) ** (n - k))

def create_bernstein_tensors(q, degree=3, q_min=None, q_max=None):
    """
    Create tensor-product Bernstein basis for points q.
    If q_min/q_max are provided, normalize accordingly (used for applying corrections).
    """
    q = np.asarray(q)
    N = q.shape[0]

    if q_min is None:
        q_min = q.min(axis=0)
    if q_max is None:
        q_max = q.max(axis=0)
    
    scale = np.maximum(q_max - q_min, 1e-8)
    u = (q - q_min) / scale  # normalize to [0,1]

    # Compute 1D Bernstein basis for each axis
    Bx = np.array([bernstein_polynomial(degree, i, u[:,0]) for i in range(degree+1)]).T
    By = np.array([bernstein_polynomial(degree, j, u[:,1]) for j in range(degree+1)]).T
    Bz = np.array([bernstein_polynomial(degree, k, u[:,2]) for k in range(degree+1)]).T

    # Compute tensor-product 3D basis
    tensors = np.einsum('pi,pj,pk->pijk', Bx, By, Bz).reshape(N, -1)
    
    return tensors, q_min, q_max

def bernstein_distortion_correction(C_measured_frames, C_expected_frames, degree=3, lambda_reg=0.0):
    """
    Fit Bernstein polynomial distortion model.
    lambda_reg > 0 adds Tikhonov regularization to prevent overcorrection.
    """
    C_meas_all = np.vstack(C_measured_frames)
    C_exp_all = np.vstack(C_expected_frames)

    F, q_min, q_max = create_bernstein_tensors(C_meas_all, degree)

    correction_models = {
        'degree': degree,
        'points_min': q_min,
        'points_max': q_max,
        'coefficients': []
    }

    # Fit each coordinate
    for coord in range(3):
        P = C_exp_all[:, coord]
        if lambda_reg > 0:
            # Regularized least squares (ridge)
            C = np.linalg.solve(F.T @ F + lambda_reg * np.eye(F.shape[1]), F.T @ P)
        else:
            C, _, _, _ = np.linalg.lstsq(F, P, rcond=None)
        correction_models['coefficients'].append(C)

    return correction_models

def applyBernsteincorrection(q, correction_models):
    """
    Apply the fitted Bernstein correction to new points q.
    """
    q = np.asarray(q)
    if q.ndim == 1:
        q = q.reshape(1, -1)

    F, _, _ = create_bernstein_tensors(
        q,
        degree=correction_models['degree'],
        q_min=correction_models['points_min'],
        q_max=correction_models['points_max']
    )

    corrected = np.zeros_like(q)
    for coord in range(3):
        C = correction_models['coefficients'][coord]
        corrected[:, coord] = F @ C

    return corrected

# ======================== FILE PARSERS ===========================
def parseCalbody(path):
    """Read calbody file: known d, a, c marker coordinates."""
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    Nd, Na, Nc = map(int, lines[0].split(",")[:3])
    d = [[float(x) for x in lines[1 + i].split(",")] for i in range(Nd)]
    a = [[float(x) for x in lines[1 + Nd + i].split(",")] for i in range(Na)]
    c = [[float(x) for x in lines[1 + Nd + Na + i].split(",")] for i in range(Nc)]
    return np.array(d), np.array(a), np.array(c)

def parseCalReadings(path):
    """Read per-frame measured D, A, C marker coordinates."""
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    Nd, Na, Nc, Nf = map(int, lines[0].split(",")[:4])

    D_frames, A_frames, C_frames = [], [], []
    idx = 1
    for _ in range(Nf):
        D = np.array([[float(x) for x in lines[idx + i].split(",")] for i in range(Nd)], float); idx += Nd
        A = np.array([[float(x) for x in lines[idx + i].split(",")] for i in range(Na)], float); idx += Na
        C = np.array([[float(x) for x in lines[idx + i].split(",")] for i in range(Nc)], float); idx += Nc
        D_frames.append(D); A_frames.append(A); C_frames.append(C)
    return D_frames, A_frames, C_frames, Nd, Na, Nc, Nf

def parseEmpivot(path):
    """Read EM pivot file: G marker frames."""
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    Ng, Nf = map(int, lines[0].split(",")[:2])
    frames, idx = [], 1
    for _ in range(Nf):
        G = np.array([[float(x) for x in lines[idx + i].split(",")] for i in range(Ng)], float)
        idx += Ng
        frames.append(G)
    return frames, Ng, Nf

def parseOptpivot(path):
    """Read optical pivot file: D (base) and H (probe) frames."""
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    Nd, Nh, Nf = map(int, lines[0].split(",")[:3])
    D_frames, H_frames, idx = [], [], 1
    for _ in range(Nf):
        D = np.array([[float(x) for x in lines[idx + i].split(",")] for i in range(Nd)], float); idx += Nd
        H = np.array([[float(x) for x in lines[idx + i].split(",")] for i in range(Nh)], float); idx += Nh
        D_frames.append(D); H_frames.append(H)
    return D_frames, H_frames, Nd, Nh, Nf

def parseCTFiducials(path):
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    Nb = int(lines[0].split(",")[0])
    b = [[float(x) for x in lines[1 + i].split(",")] for i in range(Nb)]
    return np.array(b), Nb

def parseEMFiducials(path):
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    Ng, Nb = map(int, lines[0].split(",")[:2])
    G_frames, idx = [], 1
    for _ in range (Nb):
        G = np.array([[float(x) for x in lines[idx + i].split(",")] for i in range(Ng)], float); idx += Ng
        G_frames.append(G)
    return G_frames, Ng, Nb
def parseEMnav(path):
    with open(path, "r") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    Ng, Nf = map(int, lines[0].split(",")[:2])
    G_frames, idx = [], 1
    for _ in range (Nf):
        G = np.array([[float(x) for x in lines[idx + i].split(",")] for i in range(Ng)], float); idx += Ng
        G_frames.append(G)
    return G_frames, Ng, Nf

    
# ==================== CORE COMPUTATION (PA2) =====================
def compute_expected_C(d, a, c, D_frames, A_frames):
    """
    Compute expected positions of C markers for each frame j using:
    F_D^-1 * F_A = [R_d.T @ R_a | R_d.T @ p_a - R_d.T @ p_d].
    """
    C_expected_frames = []
    for Dj, Aj in zip(D_frames, A_frames):
        R_d, p_d = point2point_3Dregistration(d, Dj)
        R_a, p_a = point2point_3Dregistration(a, Aj)
        R_d_inv, p_d_inv = R_d.T, -R_d.T @ p_d
        R_comb, p_comb = R_d_inv @ R_a, R_d_inv @ p_a + p_d_inv
        Cj = np.array([frame_transformation(R_comb, ci, p_comb) for ci in c], float)
        C_expected_frames.append(Cj)
    return C_expected_frames


def improved_pivot_calibration(frames, correction_models):
    print(f"  Raw frames shape: {np.array(frames).shape}")
    
    corrected_frames = [applyBernsteincorrection(F, correction_models) for F in frames]
    
    # Check the magnitude of correction
    raw_points = np.vstack(frames)
    corrected_points = np.vstack(corrected_frames)
    correction_magnitude = np.linalg.norm(corrected_points - raw_points, axis=1)
    print(f"  Average correction magnitude: {np.mean(correction_magnitude):.3f} mm")
    print(f"  Max correction magnitude: {np.max(correction_magnitude):.3f} mm")
    
    # Check if corrections are reasonable (should be < 10mm typically)
    if np.max(correction_magnitude) > 50:
        print("  WARNING: Excessive distortion correction detected!")
    
    return pivot_calibration(corrected_frames)


def compute_fiducial_positions(em_fiducial_frames, correction_models, em_pivot, local_markers):
    corrected_frames = [applyBernsteincorrection(frame, correction_models) for frame in em_fiducial_frames]

    B_frames = []
    for frame in corrected_frames:
        R, t = point2point_3Dregistration(local_markers, frame)
        B = frame_transformation(R,em_pivot , t)
        B_frames.append(B)
    return np.array(B_frames)


def process_navigation_data(G_nav_frames, correction_models, em_pivot, R_reg, p_reg, local_markers):
    """
    Process navigation frames to compute tip positions in CT coordinates.
    """
    corrected_frames = [applyBernsteincorrection(frame, correction_models) for frame in G_nav_frames]

    tip_positions_ct = []
    for frame in corrected_frames:
        # Apply distortion correction
        R,t = point2point_3Dregistration(local_markers, frame)
        tip_em = frame_transformation(R, em_pivot, t)
        # Transform to CT coordinates
        tip_ct = frame_transformation(R_reg, tip_em, p_reg)
        tip_positions_ct.append(tip_ct)
    
    return np.array(tip_positions_ct)

def write_output_file(filename, tip_positions_ct):
    """Write output file in required format for each dataset."""
    with open(filename, "w") as f:
        f.write(f"{len(tip_positions_ct)},{filename}\n")
        for tip in tip_positions_ct:
            f.write(f"{tip[0]:.2f},{tip[1]:.2f},{tip[2]:.2f}\n")


# ======================= DATASET PIPELINE ========================
def process_dataset(data_prefix):
    """Run full calibration pipeline for a single dataset prefix."""
    print(f"Processing {data_prefix}...")

    try:
        # --- 1. Parse input files ---
        d, a, c = parseCalbody(os.path.join(extract_path, f"{data_prefix}-calbody.txt"))
        D_frames, A_frames, C_frames, Nd, Na, Nc, Nf_cal = parseCalReadings(
            os.path.join(extract_path, f"{data_prefix}-calreadings.txt"))
        G_frames, Ng, Nf_em = parseEmpivot(os.path.join(extract_path, f"{data_prefix}-empivot.txt"))
        D_frames_opt, H_frames, Nd_opt, Nh, Nf_opt = parseOptpivot(
            os.path.join(extract_path, f"{data_prefix}-optpivot.txt"))
        b, Nb = parseCTFiducials(os.path.join(extract_path, f"{data_prefix}-ct-fiducials.txt"))
        G_emfiducial, Ng_em, Nb_em = parseEMFiducials(os.path.join(extract_path, f"{data_prefix}-em-fiducialss.txt"))
        G_emNav, Ng_nav, Nf_nav = parseEMnav(os.path.join(extract_path, f"{data_prefix}-EM-nav.txt"))

        # --- 2. Compute expected C positions ---
        print("  Computing expected C positions...")
        C_expected_frames = compute_expected_C(d, a, c, D_frames, A_frames)
        print(f"C_Expected for {data_prefix}")

        print("  Computing expected C positions...")
        C_expected_frames = compute_expected_C(d, a, c, D_frames, A_frames)
        
        
        # --- 3. Distortion Model ---
    
        print("  Building distortion correction model...")
        distortion_models= bernstein_distortion_correction(C_frames, C_expected_frames, degree =5) 

        # --- 4. Pivot Calibration ---
        print("  Performing em pivot calibration...")
        em_pivot, local_markers = improved_pivot_calibration(G_frames, distortion_models)
        #em_pivot = improved_pivot_calibration(G_frames, distortion_models)
        #pivot_diff = np.linalg.norm(em_pivot - em_pivot_notfixed)
        #print(f"EM Pivot not fixed - {em_pivot_notfixed}")
        #print(f"Fixed EM Pivot {em_pivot}")
        #print(f" Pivt difference:{pivot_diff:.2f} mm")

        # --- 5. Fiducial Positions ---
        B_em = compute_fiducial_positions(G_emfiducial, distortion_models, em_pivot, local_markers)

        # --- 6. Registration ---
        R_reg, p_reg = point2point_3Dregistration(B_em, b)

        # --- 7. Navigation ---
        tip_positions_ct = process_navigation_data(G_emNav, distortion_models, em_pivot, R_reg, p_reg, local_markers)

        # --- 8. Write results ---
        out_name = f"{data_prefix}-myoutput2.txt"
        write_output_file(out_name, tip_positions_ct)

        print(f"   Completed {data_prefix}")
        print(f"    Output: {out_name}")

    except Exception as e:
        print(f"  Error processing {data_prefix}: {e}")


def main():
    """Main loop — process all datasets in directory."""
    for fname in datasets:
        if fname.endswith("-calbody.txt"):
            data_prefix = fname[:-len("-calbody.txt")]
            process_dataset(data_prefix)
        else:
            print(f"Skipping {fname}")

    

if __name__ == "__main__":
    main()