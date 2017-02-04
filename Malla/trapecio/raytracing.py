import numpy as np
import matplotlib.pyplot as plt
import sys

w = 400
h = 300

def normalize(x):
    x /= np.linalg.norm(x)
    return x

def intersect_plane(O, D, P, N):
    # Return the distance from O to the intersection of the ray (O, D) with the 
    # plane (P, N), or +inf if there is no intersection.
    # O and P are 3D points, D and N (normal) are normalized vectors.
    denom = np.dot(D, N)
    if np.abs(denom) < 1e-6:
        return np.inf
    d = np.dot(P - O, N) / denom
    if d < 0:
        return np.inf
    return d

def intersect_sphere(O, D, S, R):
    # Return the distance from O to the intersection of the ray (O, D) with the 
    # sphere (S, R), or +inf if there is no intersection.
    # O and S are 3D points, D (direction) is a normalized vector, R is a scalar.
    a = np.dot(D, D)
    OS = O - S
    b = 2 * np.dot(D, OS)
    c = np.dot(OS, OS) - R * R
    disc = b * b - 4 * a * c
    if disc > 0:
        distSqrt = np.sqrt(disc)
        q = (-b - distSqrt) / 2.0 if b < 0 else (-b + distSqrt) / 2.0
        t0 = q / a
        t1 = c / q
        t0, t1 = min(t0, t1), max(t0, t1)
        if t1 >= 0:
            return t1 if t0 < 0 else t0
    return np.inf
    
def intersect_triangle(O, D, (v0, v1, v2)):
	# http://research.omicsgroup.org/index.php/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
	
	epsilon = sys.float_info.epsilon
	
	# Find vectors for two edges sharing V1
	e1 = v1 - v0
	e2 = v2 - v0
	
	# Begin calculating determinant - also used to calculate u parameter
	P = np.cross(D, e2)
	
	# if determinant is near zero, ray lies in plane of triangle
	det = np.dot(e1,P)
	
	# NOT CULLING
	if det > -epsilon and det < epsilon:
		return np.inf
		
	inv_det = 1. / det
	
	# calculate distance from V1 to ray origin
	T = O - v0
	
	# calculate u parameter and test bound
	u = np.dot(T,P) * inv_det
	
	# the intersection lies outside of the triangle
	if u < 0. or u > 1.:
		return np.inf
	
	# prepare to test v parameter
	Q = np.cross(T, e1)
	
	# calculate V parameter and test bound
	v = np.dot(D,Q) * inv_det
	
	# the intersection lies outside of the triangle
	if v < 0. or u + v > 1.:
		return np.inf

	t = np.dot(e2,Q) * inv_det
	
	if t < epsilon:
		return np.inf

	return t

def intersect(O, D, obj):
	if obj['type'] == 'plane':
		return intersect_plane(O, D, obj['position'], obj['normal'])
	elif obj['type'] == 'sphere':
		return intersect_sphere(O, D, obj['position'], obj['radius'])
	elif obj['type'] == 'triangle':
		return intersect_triangle(O, D, obj['position'])

def get_normal(obj, M):
    # Find normal.
    if obj['type'] == 'sphere':
        N = normalize(M - obj['position'])
    elif obj['type'] == 'plane' or obj['type'] == 'triangle':
        N = obj['normal']
    return N
    
def get_color(obj, M):
    color = obj['color']
    if not hasattr(color, '__len__'):
        color = color(M)
    return color

def trace_ray(rayO, rayD):
    # Find first point of intersection with the scene.
    t = np.inf
    for i, obj in enumerate(scene):
        t_obj = intersect(rayO, rayD, obj)
        if t_obj < t:
            t, obj_idx = t_obj, i
    # Return None if the ray does not intersect any object.
    if t == np.inf:
        return
    # Find the object.
    obj = scene[obj_idx]
    # Find the point of intersection on the object.
    M = rayO + rayD * t
    # Find properties of the object.
    N = get_normal(obj, M)
    color = get_color(obj, M)
    toO = normalize(O - M)
    # Start computing the color.
    col_ray = ambient
    for tup in light_position_and_color:
		L, color_light = tup
		toL = normalize(L - M)
		# Shadow: find if the point is shadowed or not.
		l = [intersect(M + N * .0001, toL, obj_sh) 
				for k, obj_sh in enumerate(scene) if k != obj_idx]
		if l and min(l) >= np.inf:
			# Lambert shading (diffuse).
			col_ray += obj.get('diffuse_c', diffuse_c) * max(np.dot(N, toL), 0) * color
			# Blinn-Phong shading (specular).
			col_ray += obj.get('specular_c', specular_c) * max(np.dot(N, normalize(toL + toO)), 0) ** specular_k * color_light
    return obj, M, N, col_ray


def add_sphere(position, radius, color):
    return dict(type='sphere', position=np.array(position), 
        radius=np.array(radius), color=np.array(color), reflection=.5)
    
def add_plane(position, normal):
    return dict(type='plane', position=np.array(position), 
        normal=np.array(normal),
        color=lambda M: (color_plane0 
            if (int(M[0] * 2) % 2) == (int(M[2] * 2) % 2) else color_plane1),
        diffuse_c=.75, specular_c=.5, reflection=.25)
        
        
#TRIANGLES
def add_triangle(v0, v1, v2, color):
    n = np.cross(np.subtract(v1,v0),np.subtract(v2,v0))

    return dict(type='triangle',
		position=np.array([v0,v1,v2]),
		normal=n,
		color=np.array(color),
		reflection=.5)

def mesh_triangles(Points, color):
    meshtriangles = [add_triangle(Points[0], Points[1], Points[2], color)]
    
    for i in range(3,len(Points)):
        if i%2 != 0:
            meshtriangles = meshtriangles + [add_triangle(Points[i], Points[i-1], Points[i-2], color)]
        else:
            meshtriangles = meshtriangles + [add_triangle(Points[i], Points[i-2], Points[i-1], color)]
    
    return meshtriangles
    

# List of objects.
color_plane0 = 1. * np.ones(3)
color_plane1 = 0. * np.ones(3)
scene = [add_sphere([.75, .1, 1.], .6, [0., 0., 1.]),
         add_sphere([-.75, .1, 2.25], .6, [.5, .223, .5]),
         add_sphere([-2.75, .1, 3.5], .6, [1., .572, .184]),
         add_plane([0., -.5, 0.], [0., 1., 0.]),
	]


triangles = mesh_triangles([[-1,-0.5,0.5],[-0.5,0.5,0.5],[0.,-0.5,0.5],[0.,0.5,0.5],[0.5,-0.5,0.5]], [0.,1.,1.])

for i in triangles:
    scene += [i]
 
 
# Light position.
light_position1 = np.array([5., 5., -10.])
light_position2 = np.array([-8.,7.,-10.])
light_position3 = np.array([1.,15.,-3.])

# Light color.
color_light1 = np.ones(3)
color_light2 = np.array([1., 0., 0.])
color_light3 = np.array([0., 1., 0.])

light_position_and_color = np.array([[light_position1,color_light1],[light_position2,color_light2],[light_position3,color_light3]])

# Default light and material parameters.
ambient = .05
diffuse_c = 1.
specular_c = 1.
specular_k = 50

depth_max = 5  # Maximum number of light reflections.
col = np.zeros(3)  # Current color.
O = np.array([0., 0.35, -1.])  # Camera.
Q = np.array([0., 0., 0.])  # Camera pointing to.
img = np.zeros((h, w, 3))

r = float(w) / h
# Screen coordinates: x0, y0, x1, y1.
S = (-1., -1. / r + .25, 1., 1. / r + .25)

# Loop through all pixels.
for i, x in enumerate(np.linspace(S[0], S[2], w)):
    if i % 10 == 0:
        print i / float(w) * 100, "%"
    for j, y in enumerate(np.linspace(S[1], S[3], h)):
        col[:] = 0
        Q[:2] = (x, y)
        D = normalize(Q - O)
        depth = 0
        rayO, rayD = O, D
        reflection = 1.
        # Loop through initial and secondary rays.
        while depth < depth_max:
            traced = trace_ray(rayO, rayD)
            if not traced:
                break
            obj, M, N, col_ray = traced
            # Reflection: create a new ray.
            rayO, rayD = M + N * .0001, normalize(rayD - 2 * np.dot(rayD, N) * N)
            depth += 1
            col += reflection * col_ray
            reflection *= obj.get('reflection', 1.)
        img[h - j - 1, i, :] = np.clip(col, 0, 1)

plt.imsave('fig.png', img)
