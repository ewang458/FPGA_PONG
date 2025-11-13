import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

# ===================== Utilities / Geometry =====================

#Function Definitions 
def compute_bounds(V, Indices):
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
        lam_seg = np.maximum(0, np.minimum(lambd, 1))  
        c = p + lam_seg * qp
    return c, np.sum((a - c) ** 2)


def closest_point_on_triangle(a, p, q, r):
    """
    Barycentric-based closest point on triangle pqr to point a.
    """
    a = np.asarray(a, dtype=float)
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    r = np.asarray(r, dtype=float)


    M = np.column_stack((q - p, r - p))  # 3x2 matrix
    rhs = a - p

    # Least-squares solve for λ and μ
    lam_mu, _, _, _ = np.linalg.lstsq(M, rhs, rcond=None)
    lam, mu = lam_mu
    nu = 1 - lam - mu
 
    # If point projects inside the triangle
    if lam >= 0.0 and mu >= 0.0 and (lam + mu) <= 1.0:
        c = lam * q + mu * r + nu * p
        return c, np.sum((a - c) ** 2)

    # Otherwise, the closest point is on one of the edges
    if (lam < 0):
        c, d = project_on_segment(a, r, p)
    elif (mu < 0):
        c, d = project_on_segment(a, p, q)
    else: 
        c, d = project_on_segment(a, q, r)
    return c, float(d)

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

# ================= Visualization =================================

def visualize_triangle_query(V, Indices, query, closest):
    """
    3D visualization of triangle(s), query point, and closest point.
    """
    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection="3d")
    tris = np.asarray(V)[np.asarray(Indices)]
    for tri in tris:
        loop = np.vstack((tri, tri[0]))  # close loop
        ax.plot(loop[:, 0], loop[:, 1], loop[:, 2], 'k-', linewidth=1.2)

    # Query point
    ax.scatter(query[0], query[1], query[2], color='r', label='Query', s=50)
    if closest is not None:
        ax.scatter(closest[0], closest[1], closest[2], color='g', label='Closest', s=50)
        ax.plot([query[0], closest[0]], [query[1], closest[1]], [query[2], closest[2]], 'r--', alpha=0.6)
    ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
    ax.legend()
    ax.set_title("Closest Point Visualization")
    ax.set_box_aspect([1, 1, 0.5])
    plt.show()


# =================== Tests  ===================

def print_test_result(name, point, expected, result):
    print(f"\n--- {name} ---")
    print(f"Query point: {point}")
    print(f"Expected closest: {expected}")
    print(f"Computed closest: {result}")
    if result is None:
        print("Computed closest is None (mesh empty?)")
    else:
        print(f"Distance error: {np.linalg.norm(result - expected):.6e}")
    print("-" * 40)


def brute_force_mesh_search(a, V, Indices):
    """Brute force over all triangles for correctness checking."""
    V = np.asarray(V)
    Indices = np.asarray(Indices, dtype=int)
    best_c = None
    best_d2 = float("inf")
    for tri in Indices:
        p, q, r = V[tri]
        c, d2 = closest_point_on_triangle(a, p, q, r)
        if d2 < best_d2:
            best_d2 = d2
            best_c = c
    return best_c, best_d2


def test_project_on_segment_cases(visualize=True):
    V = np.array([
        [0, 0, 0],   # p
        [1, 0, 0],   # q
        [0, 1, 0]    # r
    ])
    Indices = np.array([[0, 1, 2]])
    lower, upper, tris = compute_bounds(V, Indices)

    tests = [
        ("Edge PQ", np.array([0.5, -0.3, 0]), np.array([0.5, 0.0, 0.0])),
        ("Edge PR", np.array([-0.2, 0.3, 0]), np.array([0.0, 0.3, 0.0])),
        ("Edge QR", np.array([0.6, 0.6, 0]), np.array([0.5, 0.5, 0.0])),
        ("Vertex P", np.array([-0.3, -0.3, 0]), np.array([0.0, 0.0, 0.0])),
        ("Vertex Q", np.array([1.3, -0.2, 0]), np.array([1.0, 0.0, 0.0])),
        ("Vertex R", np.array([-0.1, 1.4, 0]), np.array([0.0, 1.0, 0.0]))
    ]

    for name, point, expected in tests:
        res, d2 = search_with_boxes(point, lower, upper, tris)
        res_bf, d2_bf = brute_force_mesh_search(point, V, Indices)
        # consistency check
        if res is None:
            raise RuntimeError("search_with_boxes returned None unexpectedly")
        assert np.allclose(res, res_bf, atol=1e-8), f"Bounding-box vs brute mismatch for {name}"
        print_test_result(name, point, expected, res)
        if visualize:
            visualize_triangle_query(V, Indices, point, res)


def test_inside_and_above_cases(visualize=True):
    V = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0]
    ])
    Indices = np.array([[0, 1, 2]])
    lower, upper, tris = compute_bounds(V, Indices)

    tests = [
        ("Inside", np.array([0.25, 0.25, 1]), np.array([0.25, 0.25, 0])),
        ("Below", np.array([0.3, 0.3, -0.8]), np.array([0.3, 0.3, 0])),
    ]

    for name, point, expected in tests:
        res, d2 = search_with_boxes(point, lower, upper, tris)
        res_bf, d2_bf = brute_force_mesh_search(point, V, Indices)
        assert np.allclose(res, res_bf, atol=1e-8)
        print_test_result(name, point, expected, res)
        if visualize:
            visualize_triangle_query(V, Indices, point, res)


def test_multiple_triangles_box_filtering(visualize=True):
    V = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [5, 5, 0],
        [6, 5, 0],
        [5, 6, 0],
    ])
    Indices = np.array([
        [0, 1, 2],  # near origin
        [3, 4, 5]   # far away
    ])
    lower, upper, tris = compute_bounds(V, Indices)

    query1 = np.array([0.2, 0.3, 0.5])
    res1, _ = search_with_boxes(query1, lower, upper, tris)
    res1_bf, _ = brute_force_mesh_search(query1, V, Indices)
    assert np.allclose(res1, res1_bf, atol=1e-8)
    print_test_result("Box Filter Test (near)", query1, np.array([0.2, 0.3, 0]), res1)
    if visualize:
        visualize_triangle_query(V, Indices, query1, res1)

    query2 = np.array([5.2, 5.3, 0.5])
    res2, _ = search_with_boxes(query2, lower, upper, tris)
    res2_bf, _ = brute_force_mesh_search(query2, V, Indices)
    assert np.allclose(res2, res2_bf, atol=1e-8)
    print_test_result("Box Filter Test (far)", query2, np.array([5.2, 5.3, 0]), res2)
    if visualize:
        visualize_triangle_query(V, Indices, query2, res2)


def test_random_points():
    """
    Quick randomized consistency checks (no visualization).
    """
    np.random.seed(42)
    V = np.random.rand(50, 3) * 2.0
    Indices = np.random.randint(0, len(V), (120, 3))
    lower, upper, tris = compute_bounds(V, Indices)
    for i in range(10):
        a = np.random.rand(3) * 2.0
        res, d2 = search_with_boxes(a, lower, upper, tris)
        res_bf, d2_bf = brute_force_mesh_search(a, V, Indices)
        visualize_triangle_query(V, Indices, a, res)
        print(f"Expected : {res_bf}")
        print(f"Calculated : {res}")
        if res is None:
            raise RuntimeError("search_with_boxes returned None unexpectedly")
        if not np.allclose(res, res_bf, atol=1e-7):
            print("Mismatch on random test", i)
            print("res:", res, "res_bf:", res_bf)
            raise AssertionError("Bounding-box vs brute mismatch on random test")
    print("Randomized consistency tests passed.")



# =================== Main runner =====================

if __name__ == "__main__":
    print("Running tests with visualization (close plot windows to continue)...")
    test_inside_and_above_cases(visualize=True)
    test_project_on_segment_cases(visualize=True)
    test_multiple_triangles_box_filtering(visualize=True)
    test_random_points()
    print("\nAll tests passed.")
