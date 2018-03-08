#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import time
from shocktubecalc import sod

#equation of state
def eos(e, V, g):
    return (g-1) * e / V

#function to look at values of quantities over box and make nice plots
def inspect(x,u,e,p,V,rho0,q,t):
    ##numerical values to see plots for, uncomment to see
    plt.plot(1/V, label = "Density")
    plt.plot(p, label = "Pressure")
    plt.plot(u, label = "Velocity")
    #plt.plot(e, label = "Energy")

    pos, reg, val = sod.solve(left_state=(pl,rhol,0), right_state=(pr,rhor,0),
                    geometry=(0,L,L/2),t=t,gamma=g,npts=500)
    
    ##analytical values to show plots for, uncomment to see 
    #plt.plot(val['x']*n/L,val['rho'],label = "Analytical Density")
    #plt.plot(val['x']*n/L,val['u'],label = "Analytical Velocity")
    #plt.plot(val['x']*n/L,val['p'],label = "Analytical Pressure")
    print(u[300], val['u'][300])

    plt.title("Sod Shock Tube at t = %s s, %d Particles"
            % ((np.around(t,2).astype(str)[:4]),n))
    plt.xlabel("Particle Number")
    plt.ylabel("Value (dimensionless)")
    plt.grid(True)
    plt.legend()
    plt.show()

#actual functional part of the code
#performs one iteration in time of all of the equations
#recalculates array one timestep in the future
def timestep(xi,ui,ei,pi,rho0,Vi,qi,dt):
    
    #make arrays for final variables
    xf = np.zeros(n)
    uf = np.zeros(n)
    ef = np.zeros(n)
    pf = np.zeros(n)
    Vf = np.zeros(n)
    qf = np.zeros(n)

    #update velocity value
    for j in range(n-1):
        uf[j] = ui[j]-(2*dt/(rho0[j]+rho0[j+1]))*(p[j+1]-p[j]+q[j+1]-q[j])/dx

        #set to 0 at boundaries
        if uf[j] <= 0 or uf[j] >= L:
            uf[j] = 0
    
    #reset end value to prevent blowing up (for every variable)
    uf[-1] = uf[-2]

    #recalculate x
    for j in range(n-1):
        xf[j] = xi[j]+(uf[j]*dt)

    xf[-1] = xf[-2]

    #recalculate V
    for j in range(n-1):
        Vf[j] = (1/rho0[j+1])*(xf[j+1]-xf[j])/dx

    Vf[-1] = Vf[-2]

    #recalculate q with appropriate conditions
    for j in range(n-1):
        if uf[j+1]-uf[j] < 0:
            qf[j+1] = .5*a*a*((uf[j+1]-uf[j])**2)*((1/Vf[j])+(1/Vi[j]))
        else:
            qf[j+1]=0

    #recalculate e
    for j in range(n-1):
        ef[j] = ei[j] - ((pi[j+1]+qf[j+1])*(Vf[j]-Vi[j]))

    ef[-1] = ef[-2]

    #recalculate p with equation of state
    for j in range(n-1):
        pf[j+1] = eos(ef[j], Vf[j], g)

    pf[-1] = pf[-2]
    
    #recalculate speed of sound at each point, determine appropriate dt
    cs = np.sqrt(g*pf*Vf)
    dt = b*np.nanmin(dx/(np.abs(uf)+cs))

    return xf,uf,ef,pf,Vf,qf,dt

if __name__ == "__main__":
    #set parameters for simulation
    cs = 1
    b = .5
    L = 2
    n = 500
    a = 1
    g = 7/5
    dt = .001
    t = 0

    #set initial conditions
    barrier = n//2
    pl = 1
    rhol = 1
    pr = .1
    rhor = .125

    #create arrays for variables
    x = np.linspace(0,L,n)
    u = np.zeros(n)
    e = np.zeros(n)
    p = np.zeros(n)
    V = np.zeros(n)
    rho0 = np.zeros(n)
    q = np.zeros(n)
    dx = x[1]-x[0]
    
    #assign initial conditions for arrays
    p[:barrier] = pl
    p[barrier:] = pr
    V[:barrier] = 1/rhol
    V[barrier:] = 1/rhor
    rho0[:barrier] = rhol
    rho0[barrier:] = rhor
    e = p * V / (g-1)

    #see starting conditions for arrays
    inspect(x,u,e,p,V,rho0,q,t)

    #iterate through many timesteps and occasionally look at result
    
    for i in range(10000):
        ti = t
        x,u,e,p,V,q,dt = timestep(x,u,e,p,rho0,V,q,dt)
        t += dt
        if t%.05 < ti%.05:
            inspect(x,u,e,p,V,rho0,q,t)
    '''
    ns = [50,100,500,1000,5000]
    for i in range(5):
        n = ns[i]
        barrier = n//2
        #create arrays for variables
        x = np.linspace(0,L,n)
        u = np.zeros(n)
        e = np.zeros(n)
        p = np.zeros(n)
        V = np.zeros(n)
        rho0 = np.zeros(n)
        q = np.zeros(n)
        dx = x[1]-x[0]
        
        #assign initial conditions for arrays
        p[:barrier] = pl
        p[barrier:] = pr
        V[:barrier] = 1/rhol
        V[barrier:] = 1/rhor
        rho0[:barrier] = rhol
        rho0[barrier:] = rhor
        e = p * V / (g-1)


        t=0
        while t<.2:
           
            x,u,e,p,V,q,dt = timestep(x,u,e,p,rho0,V,q,dt)
            t+=dt
        
        plt.plot(x/n,p)#,label = ns[i])
    #plt.xlim((0,L))
    #    plt.show()
    plt.legend()
    plt.show()
    ''' 
