import numpy as np

def fibonacci_sphere(center, r, n):
    center = np.asarray(center, dtype=float)
    i = np.arange(n)
    phi = np.pi * (3.0 - np.sqrt(5.0))        # golden angle
    y = 1 - 2 * (i + 0.5) / n                  # y from 1 to -1
    radius = np.sqrt(1 - y * y)
    theta = phi * i
    x = np.cos(theta) * radius
    z = np.sin(theta) * radius
    return center + r * np.column_stack([x, y, z])

def random_points_in_sphere(center, radius, n):
    center = np.array(center, dtype=float)
    points = []
    while len(points) < n:
        # sample in a cube, reject points outside the sphere
        p = np.random.uniform(-radius, radius, size=(3,))
        if np.linalg.norm(p) <= radius:
            points.append(center + p)
    return np.array(points)