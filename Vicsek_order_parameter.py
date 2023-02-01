import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import random
import math

def gen_data(N, L):
    """N = number of particle
       L = length of box """
    x = np.zeros(shape=N)
    y = np.zeros(shape=N)
    x_com = np.zeros(shape=N)
    y_com = np.zeros(shape=N)

    for i in range(N):
        x[i] = random.random() * L # X position
        y[i] = random.random() * L # Y position
        x_com[i] = math.cos(random.random()*2*np.pi) # x component of vlocity
        y_com[i] = math.sin(random.random()*2*math.pi) # y component of velocity
    return x, y ,x_com, y_com

def move_all(x,y,x_com,y_com,v0,L):
    # update all position
    # d is speed 
    for i in range(len(x)):
        x[i] = x[i] + x_com[i] * v0
        y[i] = y[i] + y_com[i] * v0
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

def neighbour_average_direction(x,y,x_com,y_com,k,R):
    # k = k_th neighbour
    sum_x_com = 0
    sum_y_com = 0
    for i in range(len(x)):
        if (x[i] -x[k])**2 + (y[i]-y[k])**2 <= R**2:
            sum_x_com += x_com[i]
            sum_y_com += y_com[i]
    length = math.sqrt(sum_x_com**2 + sum_y_com**2)
    return sum_x_com/length, sum_y_com/length 

l = 100      
eta = 0
R = 0.10*l # neighbour's radius
V = 0.02*l # speed
x,y,x_com,y_com = gen_data(200,L=l )

# Plotting the data
plt.rc('font', **{'family': 'serif', 'size':10})
fig,ax = plt.subplots()
ax.set_xlim(0,1000)
ax.set_ylim(-0.1,1.2)
plt.suptitle(r'$\eta = 0, N = {200}$', fontsize=20)
ax.set_xlabel("Time",fontsize=18)
ax.set_ylabel("Order Parameter",fontsize=18)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)

# variables to calculate order parameter 
time = []
order_parameter = []
order_param, = ax.plot(0,0)

# variables to calculate average noise parameter
sum = 0
no = 0

def animate(i):
    move_all(x,y,x_com,y_com,V,L=l)
    ax = np.zeros(shape=len(x))
    ay = np.zeros(shape=len(x))

    for  k in range(len(x)):
        ax[k], ay[k] = neighbour_average_direction(x,y,x_com,y_com,k,R)
        t = math.atan2(ay[k],ax[k]) + eta*np.random.uniform(-np.pi, np.pi) #(eta - 2*eta*random.random())
        ax[k] = math.cos(t)
        ay[k] = math.sin(t)
    x_com[:] = ax
    y_com[:] = ay

    # ploting order parameter
    order_parameter.append(math.sqrt(np.sum(x_com)**2 + np.sum(y_com)**2)/200)
    time.append(i)
    order_param.set_xdata(time)
    order_param.set_ydata(order_parameter)
    # print(order_parameter)

    # average order parameter
    global sum,no
    if i <= 1200:
        if i >=400:
            sum += order_parameter[i]
            no += 1

    return order_param, 

anim = animation.FuncAnimation(fig,animate,interval=1)
plt
plt.show()

print('Average order parameter = ',sum/no)

