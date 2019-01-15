## An optical ray tracer

During my sophomore year studying physics at Imperial College London, I did a Python project tracing a single and a bundle of rays propagating through different medium (refraction).


Loading packages
```python
import numpy
from numpy import *
import scipy as sp
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
```

Defining function here based on Snell's law
```python
def snell(k1, n, n1, n2):
    #kn is the incident vector, nn is the normal vector, and kn and nn are both normalised
    kn = k1/linalg.norm(k1)
    nn = n/linalg.norm(n)
    cos_theta1 = inner(nn, kn)     
    #theta_1 is the incident angle, inner() means it's finding the dot product
    
    sin_theta1 = (1-(cos_theta1)**2)**0.5
    if sin_theta1 <= n2/n1:
        sin_theta2 = (n1/n2)*sin_theta1          
        #theta_2 is the refractive angle
        cos_theta2 = (1. - sin_theta2**2)**0.5
    
        coe_1 = sin_theta1/sin_theta2
        if cos_theta1 >= 0:
            coe_2 = cos_theta2-a*cos_theta1
        else: coe_2 = -cos_theta2-coe_1*cos_theta1
        return coe_1*kn+coe_2*nn
        #coe_1 and coe_2 are just coefficents
```

Generating a bundle of rays
```python
def bun(ray,n,rmax,m):
    bundle1 = []
    for r,t in uniform(n, rmax, m):
        a = Ray(ray.r+array([r *cos(t), r*sin(t),0]),ray.dir)
        bundle1.append(a)
        
    return bundle1
```

        
Generating a sequence of radii that are uniformly distributed over a disk
```python
def uniform(n, rmax, m):       
    #rmax = the max radius of disk (float), n = number of rings of points (integer), m changes the number of points in each ring (integer)
    rstep = rmax/n
    rset = arange(0,rmax+rstep,rstep)
    ni = 1
    i = 0
    for r in rset:
        theta = 0
        tstep = 2.*pi/ni
        for t in arange(0., 2.*pi,tstep):
            yield r,t
        i += 1
        ni = m*i
```

Defining the ray
```python
class Ray:
    def __init__(self, r = 0, dir = 0): #r is position and dir is direction
        self.r = numpy.array(r)
        self.dir = numpy.array(dir)
        self.listofpositions = []
        self.listofdirections = []
        self.listofpositions.append(self.r)
        self.listofdirections.append(self.dir)
        
        
    def p(self):
        return self.r
        # position vector
        
    def k(self):
        return self.dir/linalg.norm(self.dir)
        # direction vector is normalised
       
    def append(self, p, k):
        self.listofpositions.append(p)
        self.listofdirections.append(k)
        self.r = p
        self.dir = k
        
    def vertices(self):
        return self.listofpositions
```

Propagating the ray through the optical element
```python
class OpticalElement:
    def propagate_ray(self, Ray):
```

Defining a spherical surface for the refraction to happen
```python
class SphericalRefraction(OpticalElement):
    
    def __init__(self, z0=0, curvature=0, n1=0, n2=0, aperture_radius=0):
        self.z0 = z0 
        self.curv = curvature
        self.n1 = n1
        self.n2 = n2
        self.ap_r = aperture_radius
        
    #then we define a method to find the first valid intercept of the ray"""  
    
    def intercept(self,Ray):
        #to find the intercept of the ray and the surface
        #inside_squareroot decides if the incident ray intersects with the object of refraction
        
        #first, define the quantities we need
        
        #define the centre of curvature
        centre_of_curvature = array([0.,0.,self.z0+(1/self.curv)])
        
        
        #define the displacement from the centre of curve to a position of a ray
        displacement_vector_from_centre = Ray.p() - centre_of_curvature
      
        dot_product = inner(displacement_vector_from_centre,Ray.k())
        
        if dot_product >= 0:
            return "the ray is travelling away from the spherical surface"
            
        if self.curv != 0:
            centre_of_curvature = array([0.,0.,self.z0+(1/self.curv)])
            displacement_vector_from_centre = Ray.p() - centre_of_curvature
            
            length_r = linalg.norm(displacement_vector_from_centre)
            
            inside_squareroot = (dot_product)**2 - ((length_r)**2 - (1/self.curv)**2)
            
            if inside_squareroot >= 0:
                l1 = -(dot_product) + (inside_squareroot)**0.5
                if l1 < length_r:
                    intercept1_vector = l1*Ray.k() + displacement_vector_from_centre
                    if abs(intercept1_vector[0]) > self.ap_r:
                        return "no intercept case 1"
                    elif abs(intercept1_vector[1]) > self.ap_r:
                        return "no intercept case 1"
                    else: return intercept1_vector+centre_of_curvature
                     
        
                else:
                    l2 = -(dot_product) - (inside_squareroot)**0.5
                    intercept2_vector = l2*Ray.k() + displacement_vector_from_centre
                    if abs(intercept2_vector[0]) > self.ap_r:
                        return  "no intercept case 2"
                    elif abs(intercept2_vector[1]) > self.ap_r:
                        return  "no intercept case 2"
                    else: return intercept2_vector+centre_of_curvature
            
            else: return "no intercept"
            
        else: 
            a = (self.z0 - (Ray.p()[2]))/(Ray.k()[2])
            return Ray.p() + a*Ray.k()
             
    def propagate_ray(self,Ray):
        # finding the new position and direction vectors of the ray after refraction 
        new_position = self.intercept(Ray)    
        new_displacement_vector_from_centre = new_position - array([0,0,self.z0+(1/self.curv)])
        normal = new_displacement_vector_from_centre/linalg.norm(new_displacement_vector_from_centre)
        new_direction = snell(normal,Ray.k(),self.n1,self.n2)
        Ray.append(new_position,new_direction)
        return Ray
```

Defining OutputPlane, a derived class from OpticalElement propagating rays to an ending point on the z axis
```python
class OutputPlane(OpticalElement):
    def __init__(self,z0=0):
        self.z0 = z0
        
    def intercept(self,Ray):
        # locates the intercept of a ray with the OutputPlane 
        a = (self.z0 - Ray.p()[2])/Ray.k()[2]
        return Ray.p() + a*Ray.k()
            
    def propagate_ray(self,Ray):
        #find the new position when it hits the outputplane 
        Ray.append(self.intercept(Ray),Ray.k())
        return Ray
```

Plotting a single ray
```python
P1 = Ray(r = [0.001,0,0], dir = [0,0,1])
A0 = SphericalRefraction(100.,0.03,1,1.5,100)
A0.intercept(P1)
A0.propagate_ray(P1)
P1.k()
P1.p()
A1 = OutputPlane(123)
A1.intercept(P1)
A1.propagate_ray(P1)
P1.p()
P1.k()
P1.vertices()


z = []
x = []

# drawing a 2D graph for z, x axis to find the paraxial focus
   
for i in P1.vertices():
    z.append(i[2])
    x.append(i[0])

fig=plt.figure(figsize=(20,10))    
plt.plot(z,x)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)
```

<img src="https://github.com/FlorisWu/optical-ray-tracer/blob/master/single_ray.jpg?raw=true" width="900"/>

Plotting a bundle of rays
```python
D1 = SphericalRefraction(105.,0.02,1,1.5168,100)
E1 = OutputPlane(138.555)
H1 = Ray(r = [0.001,0,80], dir = [0,0,1])
B1 = bun(H1,4,10,8)

for i in B1:
    D1.propagate_ray(i)
    E1.propagate_ray(i)
    x=[]
    z=[]
    
    for q in i.vertices():
        x.append(q[0])
        z.append(q[2])
        plt.plot(z,x)
        plt.xticks(fontsize=25)
        plt.yticks(fontsize=25)
```


<img src="https://github.com/FlorisWu/optical-ray-tracer/blob/master/bundle_of_rays.png?raw=true" width="900"/>
