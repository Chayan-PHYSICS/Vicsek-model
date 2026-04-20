import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import random
import math

def gen_data(N, L):
    """N = number of particle
       L = Box length"""
    x = np.zeros(shape=N)
    y = np.zeros(shape=N)
    com_x = np.zeros(shape=N)
    com_y = np.zeros(shape=N)

    for i in range(N):
        x[i] = random.random() * L # X position
        y[i] = random.random() * L # Y position
        com_x[i] = math.cos(random.random()*2*np.pi) # x component of vlocity
        com_y[i] = math.sin(random.random()*2*math.pi) # y component of velocity
    return x, y ,com_x, com_y

def move_all(x,y,com_x, com_y,d,L):
    # update all position
    # d is speed 
    for i in range(len(x)):
        x[i] = x[i] + com_x[i] * d
        y[i] = y[i] +  com_y[i] * d
        # Boundary condition for x
        if x[i] > L:
            x[i] -= L
        if x[i] < 0:
            x[i] = L-x[i]
        # Boundary condition for y
        if y[i] > L:
            y[i] -= L
        if y[i] < 0:
            y[i] = L-y[i]     

def neighbour_average_direction(x,y,com_x,com_y,k,R):
    # k = k_th neighbour
    sum_com_x = 0
    sum_com_y = 0
    for i in range(len(x)):
        if (x[i] -x[k])**2 + (y[i]-y[k])**2 <= R**2:
            sum_com_x += com_x[i]
            sum_com_y += com_y[i]
    length = math.sqrt(sum_com_x**2 + sum_com_y**2)
    return sum_com_x/length, sum_com_y/length 

L = 100      
eta = 0.2 # noise paremeter
R = 0.10*L # neighbour's radius
V = 0.01*L # speed
x,y,com_x,com_y = gen_data(300,L)

# Plotting the data
fig = plt.figure()
particles = plt.quiver(x,y,com_x,com_y,np.random.uniform(-np.pi,np.pi,size=len(x)))

def animate(i):
    move_all(x,y,com_x,com_y,V,L)
    ax = np.zeros(shape=len(x))
    ay = np.zeros(shape=len(x))

    for  k in range(len(x)):
        ax[k], ay[k] = neighbour_average_direction(x,y,com_x,com_y,k,R)
        # Add noise 
        t = math.atan2(ay[k],ax[k]) +  eta*np.random.uniform(-np.pi, np.pi) #(eta - 2*eta*random.random())
        ax[k] = math.cos(t)
        ay[k] = math.sin(t)
    com_x[:] = ax
    com_y[:] = ay

    pos = np.column_stack((x,y))
    particles.set_offsets(pos)
    particles.set_UVC(com_x, com_y)
    # ax1.axis('equal')
    return particles, 

anim = animation.FuncAnimation(fig,animate,np.arange(1, 200),interval=1)

plt.show()
    

    
    
      