# David Pitt
# Nov 6th, 2020
# First try - a proof-of-concept simulation of the 3-body problem.
# point mass class

import numpy as np

# All values in SI units
real_G_const = 6.67e-11 # N*m^3/kg^2
M_Earth = 5.972e26 #3e10  #3e10 

# Note: add a measure to scale distance dynamically
# or mass

class pointMass:
    # basic point mass object. 
    # a point mass is instantiated with initial conditions for
    # x,y,vx and vy.

    def __init__(self,x,y,v_naught_x,v_naught_y,m):
        self.pos_ = np.array([x,y])
        self.vel_ = np.array([v_naught_x,v_naught_y])
        self.m_ = m*M_Earth
        self.p_ = np.multiply(self.vel_, self.m_)
        self.accel_ = np.array([0.,0.])
        self.uTot_ = 0
        self.history_ = np.array([self.pos_])
        self.velHist_ = np.array([self.vel_])

    def gravitationOf(self,p2,offset):
        # method to calculate gravitational force from p2
        
        '''
        print("self:")
        print(self.pos_)
        print("p2:")
        print(p2.pos_)
        print("offset:")
        print(offset)
        '''
        #r = [(p2.pos_[0] - (self.pos_[0]+offset[0])),(p2.pos_[1] - (self.pos_[1]+offset[1]))]
        r = p2.pos_-self.pos_+offset ##THIS IS RIGHT
        #print("r = ")
        #print(r)
        r_hat = np.divide(r,np.linalg.norm(r))
        #print("rhat = ")
        #print(r_hat)
        fGrav = np.multiply(r_hat,(real_G_const*self.m_ * p2.m_) / (np.linalg.norm(r)**2))
        print("Grav =: ")
        print(fGrav)
        return [fGrav[0], fGrav[1]]

    def RK4_helper(self,dt,offset,scene): 
        # helper function for Runge-Kutta 4 that gives an incremental
        # acceleration and velocity to the planet 
        
        #y = np.array([self.vel_,self.accel_])

        self.accel_ = [0.,0.]
        dt = float(dt)
        for p in scene:
            if self != p:
                self.accel_ += np.true_divide(self.gravitationOf(p,offset),self.m_)
        print("a =:")
        print(self.accel_)
        #y[0] += np.multiply(self.vel_,dt) +  0.5*np.multiply(self.accel_,dt**2)
        velStep = 0.5*np.multiply(self.accel_,dt**2) # Update velocity
        
        return np.array([self.vel_+velStep,self.accel_])
    
    def calculateKE(self):
        # method to calculate kinetic energy
        #return (0.5*self.m_*(self.vel_[0]^2+self.vel_[1]^2))

        return (0.5*self.m_*np.linalg.norm(self.vel_)**2)

    def calculateGPE(self,scene):
        # calculate the gravitational potential
        potential = 0
        for p in scene:
            if (p.pos_[0] != self.pos_[0] and p.pos_[1] != self.pos_[1]):
                potential += (real_G_const*self.m_ * p.m_) / np.sqrt((self.pos_[1] - p.pos_[1])**2+(self.pos_[0] - p.pos_[0])**2)
        return -(potential)
    
    def calculateUTotal(self,scene):
        # Unused.
        return self.calculateKE() + self.calculateGPE(scene)