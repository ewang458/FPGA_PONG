import numpy as np

# =================== DEBUGGING HELPERS ===================
def debug_array(name, arr, verbose=True):
    """Prints info about a numpy array for debugging."""
    arr = np.asarray(arr)
    if verbose:
        print(f"{name}: shape={arr.shape}, dtype={arr.dtype}, min={arr.min():.3f}, max={arr.max():.3f}, mean={arr.mean():.3f}")
        if arr.ndim == 2 and arr.shape[1] == 3:
            norms = np.linalg.norm(arr, axis=1)
            print(f"  Norms: min={norms.min():.3f}, max={norms.max():.3f}, mean={norms.mean():.3f}")
    return arr

def check_rotation_matrix(R, tol=1e-6):
    """Check if a matrix is a valid rotation matrix."""
    RtR = R.T @ R
    I = np.eye(3)
    det = np.linalg.det(R)
    valid = np.allclose(RtR, I, atol=tol) and np.isclose(det, 1.0, atol=tol)
    print(f"Rotation check: det={det:.6f}, orthogonal={np.allclose(RtR,I,tol)} -> valid={valid}")
    return valid

def check_pivot_solution(pivot, name="Pivot"):
    """Sanity check for pivot coordinates."""
    pivot = np.asarray(pivot)
    print(f"{name}: shape={pivot.shape}, values={pivot}")
    if np.any(np.abs(pivot) > 1e4):
        print(f"  WARNING: Pivot coordinates unusually large!")

# =================== UNIT TEST FUNCTIONS ===================
def test_point2point_3Dregistration():
    """Test registration on a known rotation/translation."""
    A = np.random.rand(10, 3)
    R_true = Rot_X(np.pi/6) @ Rot_Y(np.pi/4)
    t_true = np.array([1.0, -2.0, 0.5])
    B = (R_true @ A.T).T + t_true

    R_est, t_est = point2point_3Dregistration(A, B)
    debug_array("A", A)
    debug_array("B", B)
    debug_array("t_est", t_est)
    check_rotation_matrix(R_est)

    # Check errors
    B_est = (R_est @ A.T).T + t_est
    error = np.linalg.norm(B - B_est, axis=1)
    print(f"Registration error: mean={error.mean():.6f}, max={error.max():.6f}")

def test_pivot_calibration():
    """Test pivot calibration with synthetic data."""
    pivot_true = np.array([5.0, -3.0, 2.0])
    g_local = np.random.rand(4,3)
    frames = []
    for _ in range(5):
        R = Rot_X(np.pi/8) @ Rot_Y(np.pi/12)
        t = np.random.rand(3)
        frame = (R @ g_local.T).T + pivot_true + t
        frames.append(frame)
    
    p_tip, g_local_est = pivot_calibration(frames)
    print("Pivot tip estimate:", p_tip)
    print("Local markers shape:", g_local_est.shape)

def test_applyBernsteincorrection():
    """Test Bernstein correction for known points."""
    points = np.random.rand(5,3)
    C_measured = points + 0.1*np.random.randn(*points.shape)
    C_expected = points
    model = bernstein_distortion_correction([C_measured], [C_expected], degree=3)
    corrected = applyBernsteincorrection(C_measured, model)
    debug_array("Corrected points", corrected)
    print("Expected vs corrected difference:", np.linalg.norm(corrected - C_expected, axis=1))

def test_file_parsers():
    """Test parser functions with dummy small files."""
    # Create dummy calbody content
    d = np.random.rand(2,3)
    a = np.random.rand(2,3)
    c = np.random.rand(2,3)
    print("Calbody arrays shapes:", d.shape, a.shape, c.shape)

    # For other parsers you could simulate small arrays
    # and verify shapes, number of frames, etc.
    # For example:
    D_frames = [d for _ in range(3)]
    A_frames = [a for _ in range(3)]
    C_frames = [c for _ in range(3)]
    Nd, Na, Nc, Nf = d.shape[0], a.shape[0], c.shape[0], 3
    debug_array("D_frames", np.vstack(D_frames))
    debug_array("A_frames", np.vstack(A_frames))
    debug_array("C_frames", np.vstack(C_frames))


# ===================== SYNTHETIC DATA GENERATORS =====================
def generate_synthetic_calbody():
    """Generate minimal synthetic calbody data (d, a, c markers)."""
    d = np.random.rand(2, 3)  # 2 markers for D
    a = np.random.rand(2, 3)  # 2 markers for A
    c = np.random.rand(2, 3)  # 2 markers for C
    return d, a, c

def generate_synthetic_frames(base_points, n_frames=3, noise=0.0):
    """Generate multiple frames by applying random rotations/translations."""
    frames = []
    for _ in range(n_frames):
        R = Rot_X(np.random.rand()*0.1) @ Rot_Y(np.random.rand()*0.1) @ Rot_Z(np.random.rand()*0.1)
        t = np.random.rand(3)
        frame = (R @ base_points.T).T + t + noise*np.random.randn(*base_points.shape)
        frames.append(frame)
    return frames

def generate_synthetic_fiducials(n_fiducials=3):
    """Generate synthetic EM or CT fiducials."""
    return np.random.rand(n_fiducials, 3)

# ===================== DRY-RUN PIPELINE =====================
def dry_run_pipeline(dataset_prefix="synthetic"):
    print(f"\n=== Dry-run pipeline for dataset '{dataset_prefix}' ===")
    try:
        # 1. Synthetic calbody
        d, a, c = generate_synthetic_calbody()
        debug_array("d", d)
        debug_array("a", a)
        debug_array("c", c)

        # 2. Synthetic measured frames
        D_frames = generate_synthetic_frames(d)
        A_frames = generate_synthetic_frames(a)
        C_frames = generate_synthetic_frames(c)

        # 3. Expected C positions
        C_expected_frames = compute_expected_C(d, a, c, D_frames, A_frames)
        debug_array("C_expected_frames", np.vstack(C_expected_frames))

        # 4. Distortion model
        distortion_models = bernstein_distortion_correction(C_frames, C_expected_frames, degree=3)
        print("Distortion model coefficients shapes:",
              [c.shape for c in distortion_models['coefficients']])

        # 5. EM pivot calibration
        G_frames = generate_synthetic_frames(np.random.rand(3, 3))
        em_pivot, local_markers = improved_pivot_calibration(G_frames, distortion_models)
        check_pivot_solution(em_pivot, "EM Pivot")
        debug_array("local_markers", local_markers)

        # 6. Fiducial positions
        G_emfiducial = generate_synthetic_fiducials()
        B_em = compute_fiducial_positions([G_emfiducial], distortion_models, em_pivot, local_markers)
        debug_array("B_em", B_em)

        # 7. Registration (CT)
        b = generate_synthetic_fiducials(B_em.shape[1])
        R_reg, p_reg = point2point_3Dregistration(B_em[0], b)
        check_rotation_matrix(R_reg)
        debug_array("p_reg", p_reg)

        # 8. Navigation (tip positions)
        G_emNav = generate_synthetic_frames(np.random.rand(2, 3))
        tip_positions_ct = process_navigation_data(G_emNav, distortion_models, em_pivot, R_reg, p_reg, local_markers)
        debug_array("tip_positions_ct", tip_positions_ct)

        print(f"Dry-run pipeline '{dataset_prefix}' completed successfully!\n")

    except Exception as e:
        print(f"Error during dry-run for dataset '{dataset_prefix}': {e}")



# =================== RUN ALL TESTS ===================
if __name__ == "__main__":
    print("=== Testing point2point_3Dregistration ===")
    test_point2point_3Dregistration()
    
    print("\n=== Testing pivot_calibration ===")
    test_pivot_calibration()
    
    print("\n=== Testing applyBernsteincorrection ===")
    test_applyBernsteincorrection()
    
    print("\n=== Testing file parsers / arrays ===")
    test_file_parsers()

    for prefix in ["a", "b"]:  # simulate two datasets
    dry_run_pipeline(prefix)

    #NOT RUNNABLE JUST ROUGH TESTS