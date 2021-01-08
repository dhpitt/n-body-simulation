# David Pitt
# Nov 6th, 2020
# First try - a proof-of-concept simulation of the 3-body problem.
# units are in SI, time is in years.

import numpy as np
from math import pi
import time
from scipy import integrate
import plotly.graph_objects as go
import matplotlib.pylab as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as anim

from scipy.integrate import odeint
import plotly.graph_objects #as go
from three_body_planet import pointMass

# set up the simulation

planet_1 = pointMass(0,1e7,0,0,1.0)
planet_2 = pointMass(0,-3e7,9e4,0,0.012)
planet_3 = pointMass(0.0,-7e6,-4e4,1e5,2.1)
planet_4 = pointMass(-8e6,-3e6,1e4,8e5,1.25)

planets_2b = []
planets_2b.append(planet_1)
planets_2b.append(planet_2)

planets_3b = []
planets_3b.append(planet_1)
planets_3b.append(planet_2)
planets_3b.append(planet_3)


planets_4b = []
planets_4b.append(planet_1)
planets_4b.append(planet_2)
planets_4b.append(planet_3)
planets_4b.append(planet_4)


def RK4_solve_3body(t_i,t_f,N,planets): 

    # simulates gravitation and orbits of 3 celestial bodies

    Ts = np.linspace(1,((t_f - t_i) * N - 1),((t_f - t_i) * N - 1))
    Ts = Ts.astype(int)

    dt = (t_f - t_i) / N

    for p in planets:
            hPlace = p.history_[0]
            p.history_ = np.zeros([N*(t_f-t_i),2])
            p.history_[0] = hPlace
            vPlace = p.velHist_[0]
            p.velHist_ = np.zeros([N*(t_f-t_i),2])
            p.velHist_[0] = vPlace

    for t in Ts:
        for p in planets:
            RK4_step(t,dt,p,planets)
            p.history_[t] = p.pos_
            p.velHist_[t] = p.vel_

    return planet_1.history_,planet_2.history_,planet_3.history_,Ts

def RK4_solve_nBody(t_i,t_f,N,planets): 

    Ts = np.linspace(1,((t_f - t_i) * N - 1),((t_f - t_i) * N - 1))
    Ts = Ts.astype(int)
    #print(Ts)
    #print(len(Ts))
    numSteps = (t_f - t_i) / N 
    dt = numSteps #*86400
    #print("dt = " + str(dt))
    for p in planets:
            hPlace = p.history_[0]
            p.history_ = np.zeros([N*(t_f-t_i),2])
            p.history_[0] = hPlace
            vPlace = p.velHist_[0]
            p.velHist_ = np.zeros([N*(t_f-t_i),2])
            p.velHist_[0] = vPlace

    for t in Ts:
        for p in planets:
            RK4_step(t,dt,p,planets)
            p.history_[t] = p.pos_
            p.velHist_[t] = p.vel_

    #return planet_1.history_,planet_2.history_,planet_3.history_,planet_4.history_,Ts
    return [p.history_ for p in planets],Ts

def RK4_step(t,dt,body,scene):
    # An adaptation of Runge-Kutta 4th Order
    # that outputs the next state of a planet (position and velocity)
    # given current state, history, and the scene.

    y_old = np.array([body.history_[t-1],body.velHist_[t-1]])

    k1 = np.multiply(body.RK4_helper(dt, [0,0], scene), dt) 
    k2 = np.multiply(body.RK4_helper(dt, np.multiply(0.5,k1), scene), dt)
    k3 = np.multiply(body.RK4_helper(dt, np.multiply(0.5,k2), scene), dt)
    k4 = np.multiply(body.RK4_helper(dt, k3, scene), dt)
    
    y_new = y_old + np.multiply((1/6.),(k1 + np.multiply(k2,2) + np.multiply(k3,2) + k4))

    body.pos_ = y_new[0]
    body.vel_ = y_new[1]

#a,b,c,tspace = RK4_solve_3body(0,20,100,planets_3b) # 3 body system output

#results,tspace = RK4_solve_4body(0,20,50,planets_4b) # 4 body system output
#results,tspace = RK4_solve_nBody(0,20,50,planets_2b) # 4 body system output
results,tspace = RK4_solve_nBody(0,20,50,planets_3b) # 4 body system output


a = results[0]
b = results[1]
c = results[2]
'''

d = results[3]
'''


#print(a)
#print(tspace)

# Animate the simulation

plt.style.use('seaborn-pastel')
fig, ax = plt.subplots()
ax.set(xlim=(-1e8,1e8),ylim=(-1e8,1e8))

d1, = plt.plot(a[0][0],a[0][1],'ro',color='blue')
d2, = plt.plot(b[0][0],b[0][1],'ro',color='green')
d3, = plt.plot(c[0][0],c[0][1],'ro')
#d4, = plt.plot(d[0][0],d[0][1],'ro',color='black')


def animate3b(i):
    d1.set_data(a[i][0],a[i][1])
    d2.set_data(b[i][0],b[i][1])
    d3.set_data(c[i][0],c[i][1])
    return d1,d2,d3,

'''
def animate4b(i):
    d1.set_data(a[i][0],a[i][1])
    d2.set_data(b[i][0],b[i][1])
    d3.set_data(c[i][0],c[i][1])
    d4.set_data(d[i][0],d[i][1])

    return d1,d2,d3,d4,
'''   

anim = anim.FuncAnimation(fig, animate3b, frames=tspace, \
                                      interval=50, blit=False, repeat=True)
#plt.draw()
plt.draw()
anim.save('3body1.mp4', fps=60, extra_args=['-vcodec', 'libx264'])
