# -*- coding: utf-8 -*-
"""
Reaction-Diffusion in 1-dimension
u_t = d_u*Lapace u + .1 - u + u*u*v
v_t = d_v*Lapace v + 1 -  u*u*v
Neumann BC
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import random as rnd
#from matplotlib import colors
from fun_VanGor import Laplace1d, Nl_Laplace1d, Nr_Laplace1d # cmap, norm, bounds


d_u = .0002
d_v = .02
uast = 1.1
vast = pow(1/uast,2)


gs = 200 #changing to 400 does not chenge the plot much. What is the relation to the choice of dh?


dh = 5*pow(10,-3) #pow(1/(gs*100),1) #at least .5*10^6  steps
dt = pow(dh,2)
steps = math.floor(pow(10,4)/dt)
"""  10 000 = steps *dt
dt = .0001
dh = dt#1/gs

dh = pow(1/(gs*100),1/2) # gs=400 dh = 5*10^(-4) then #at least .5 milion steps
dt = pow(dh,2)
"""
#xi_u = rnd.uniform(0,1)*pow(10,-2)*3
#xi_v = rnd.uniform(0,1)*pow(10,-2)*3
"""
dh = 1/gs #1/gs   dt/dh^2 < 1/0.0004 < 1/.04    dt < dh^2/0.0004   dt = dh^2/0.0008
dt = pow(dh,2)/0.0008#.0001. calculate the point in time using D dt / dh^2 = .5
"""

Grid_u = [uast+rnd.uniform(0,1)*pow(10,-3) for i in range(0,gs)]  #np.array([uast+xi_u]*(gs*gs)).reshape(gs,gs)
Grid_v = [vast+rnd.uniform(0,1)*pow(10,-3) for i in range(0,gs)]  #np.array([vast+xi_v]*(gs*gs)).reshape(gs,gs)
System = {'u': Grid_u, 'v': Grid_v }
Next_u = [uast]*gs 
Next_v = [vast]*gs 


steps_plot = pow(10,5)
for t in range(0,steps):
        printtime = t/steps_plot
        if printtime.is_integer() == True:
            plt.plot([i/gs for i in range(len(System['u']))],System['u'])
            plt.title('u: t = ' + str(round(t *dt,2)))
            plt.show()
        for i in range(1,gs-1):
                Next_u[i] = System['u'][i] + (pow(dh,-2)*d_u*Laplace1d(System['u'],i)+.1-System['u'][i] +pow(System['u'][i] ,2)*System['v'][i] )*dt
                Next_v[i] = System['v'][i]  + (pow(dh,-2)*d_v*Laplace1d(System['v'],i)+1-pow(System['u'][i] ,2)*System['v'][i] )*dt
        for i in [0]:
                Next_u[i] = System['u'][i] + (pow(dh,-2)*d_u*Nl_Laplace1d(System['u'],i)+.1-System['u'][i] +pow(System['u'][i] ,2)*System['v'][i] )*dt
                Next_v[i] = System['v'][i]  + (pow(dh,-2)*d_v*Nl_Laplace1d(System['v'],i)+1-pow(System['u'][i] ,2)*System['v'][i] )*dt
        for i in [gs-1]:
                Next_u[i] = System['u'][i] + (pow(dh,-2)*d_u*Nr_Laplace1d(System['u'],i)+.1-System['u'][i] +pow(System['u'][i] ,2)*System['v'][i] )*dt
                Next_v[i] = System['v'][i]  + (pow(dh,-2)*d_v*Nr_Laplace1d(System['v'],i)+1-pow(System['u'][i] ,2)*System['v'][i] )*dt

        System['u'], Next_u = Next_u, System['u'] #swap lists
        System['v'], Next_v = Next_v, System['v'] #swap lists

"""
ax1 = plt.axes(frameon=False)
plt.imshow([bounds], cmap=cmap, norm=norm)
plt.title('Values')
plt.xticks(list(range(len(bounds))), [round(bounds[i],1) for i in range(len(bounds))])
ax1.axes.get_yaxis().set_visible(False)
"""

""" #Discrete:
    #PDE:   
  
        
  Check the new system steady state and add colors number graphs.
    
   v=1/u^2
   u =1.1
   v =pow(1/u,2)  
   
   Van's
d_u = .0002
d_v = .02
uast = 1.1
vast = pow(1/uast,2)
xi_u = rnd.uniform(0,1)*pow(10,-1)
xi_v = rnd.uniform(0,1)*pow(10,-1)

Grid_u = [[uast+xi_u for i in range(0,gs)]
         for i in range(0,int(gs))]
Grid_v = [[vast + xi_v for i in range(0,gs)]
        for i in range(0,int(gs))]


 
    """

   