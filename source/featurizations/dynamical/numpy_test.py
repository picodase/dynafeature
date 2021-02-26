# %%
import numpy as np
# %%
dmin = 0
dmax = 10
steps = 12j

a = np.mgrid[dmin:dmax:steps,dmin:dmax:steps,dmin:dmax:steps]
a.shape
# %%
for x in np.nditer(a,itershape=(3,4,4,4)):
    print(x)

# %%
mg = np.meshgrid(np.arange(0,10,1), np.arange(0,10,1),np.arange(0,10,1))
# %%
dmin = 0
dmax = 10
steps = 3j

a = np.mgrid[dmin:dmax:steps,dmin:dmax:steps,dmin:dmax:steps]
a
# %%
fx1y1 = 

a00 = fx1y1
a10 = fx2y1 - fx1y1
a01 = fx1y2 - fx1y1
a11 = fx1y1 - fx2y1 - fx1y2 + fx2y2

b00 = gx1y1
b10 = gx2y1 - gx1y1
b01 = gx1y2 - gx1y1
b11 = gx1y1 - gx2y1 - gx1y2 + gx2y2

c00 = a11*b00 - a00*b11
c10 = a11*b10 - a10*b11
c01 = a11*b01 - a01*b11

# (-a11*c10)s^2 + (-a11*c00 - a01*c10 + a10*c01) s + (a00*c01 - a01*c00)=0
a = (-a11*c10)
b = (-a11*c00 - a01*c10 + a10*c01)
c = (a00*c01 - a01*c00)

sl = (-b - sqrt(b^2 - 4*a*c))/(2*a)
sh = (-b + sqrt(b^2 - 4*a*c))/(2*a)

t = (-c00/c01) - (c10/(c01*s))

