# OLG_krusell_smith_dyn.py
"""
Created on January 24, 2022

@author: Burkhard Heer

purpose: solves the stochastic Overlapping Generations Model 
        in Chapter 10.2.2 of Heer/Maussner, 
        'Dynamic General Equilibrium Modeling: Computational Methods
         and Applications', 3rd edition (scheduled for 2022)

        part 2: computes the dynamics
        
        FIRST RUN 'OLG_krusell_smith_ss.py' WHICH COMPUTES
        THE STEADY STATE (INPUT INTO THIS PROGRAM)!!!!!!!
        
method: Krusell-Smith Algorithm         

"""



# Part 1: import libraries
import pandas as pd
import numpy as np
import scipy.optimize 
import quantecon as qe
from scipy.stats import norm
from scipy import interpolate
import time
import math
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.iolib.summary2 import summary_col


# part 2: definition of the functions
#
# wage function
def wage(k,l,z):
    return (1-alpha) * z* k**alpha * l**(-alpha)

# interest rate function
def interest(k,l,z):
    return alpha * z* k**(alpha - 1) * l**(1-alpha)-delta

# production function
def production(k,l,z):
    return z*k**(alpha) * l**(1-alpha)

# utility function 
def u(c,l): 
    if eta==1:
        y =  gamma*log(c) +(1-gamma)*log(1-l) 
    else:
        y = (c**(gamma*(1-eta)) *(1-l)**((1-gamma)*(1-eta))) / (1-eta)
    return y


# marginal utility from consumption
def uc(c,l):
	y = gamma*c**(gamma*(1-eta)-1) * (1-l)**((1-gamma)*(1-eta))
	return y

# marginal disutility from labor
def ulabor(c,l):
	y = (1-gamma)*c**(gamma*(1-eta)) * (1-l)**((1-gamma)*(1-eta)-1) 
	return y

# computes the gini for distribution where x has measure g(x) 
def ginid(x,gy,ng):
    x =  np.where(x<=0, 0, x)
    xmean = np.dot(x,gy)
#    y = np.c_[x,gy]
    # sorting according to first column still has to be implemented
#    y0 = np.sort(y,axis=0)
#    x = y0[:,0]
#    g = y0[:,1]
    g = gy
    f = np.zeros(ng)  # accumulated frequency */
    f[0] = g[0]*x[0]/xmean
    gini = 1- f[0] * g[0]
    for i in range(1,ng):
        f[i] = f[i-1] + g[i]*x[i]/xmean
        gini = gini - (f[i]+f[i-1])*g[i]
	
    return gini


#  prod1(x)
#  function to compute individual productivity levels
#
#  two conditions:
#  1. average productivity equal to one
#  2. log variance equal to vartheta
#
def prod1(x,*args):
    z1=x[0]
    z2=x[1]
    y=np.zeros(2)
    y[0]=z1+z2-2.0   # normalization to one 
    # variance of log productivity 
    varth=args[0]
    varz=( (log(z1))**2 + (log(z2))**2 ) / 2.0  
    y[1]=varz-varth
    return y

# getpolicy()
# computes the optimal policy functions in steady state
def getpolicy():
    apolicy = np.zeros((nperm,ntheta,na,nage,nk,nz))  # next-period wealth at age iage,
                                        # wealth ia, permanent prod iperm, 
                                        # stochastic prod itheta
    cpolicy = np.zeros((nperm,ntheta,na,nage,nk,nz))  # consumption
    lpolicy = np.zeros((nperm,ntheta,na,nage,nk,nz))  # labor supply
     
    for ik in range(nk):
        k0 = kgrid[ik]
        for iz in range(nz):
            z0 = zz[iz]
            x0 = np.array([1.0,iz,log(k0),iz*log(k0)])
            k1 = np.dot(d,x0.T)
            k1 = exp**k1
            n0 = np.dot(f,x0.T) # dot product of two vectors f and x0
            n0 = exp**n0
            trbar0 = np.dot(dtr,x0.T)
            w0 = wage(k0,n0,z0)
            r0 = interest(k0,n0,z0)
            pen0=replacement_ratio*(1-taulbar)*lbar0*w0
            
            
            # last period nage of life
            for ia in range(na):
                a0 = agrid[ia]
                for iperm in range(nperm):
                    if pentype==0:
                        pene=pen0
                    else:
                        pene = pen0*perm[iperm]
                
                    for itheta in range(ntheta):
                        c= (1+(1-tauk)*r0)*a0+pene+trbar0 
                        c = c/(1+tauc)
                        apolicy[iperm,itheta,ia,nage-1,ik,iz] = 0
                        cpolicy[iperm,itheta,ia,nage-1,ik,iz] = c
                        lpolicy[iperm,itheta,ia,nage-1,ik,iz] = 0
            
    
    # computation of the decision rules for the retired 
    # age nage, nage-1,--,nw+1
    
    for iage in range(nage-2,nw-1,-1):
        print(q,iage)
        for ik in range(nk):
            k0 = kgrid[ik]
            for iz in range(nz):
                z0 = zz[iz]
                x0 = np.array([1.0,iz,log(k0),iz*log(k0)])
                k1 = np.dot(d,x0.T)
                k1 = exp**k1
                n0 = np.dot(f,x0.T) # dot product of two vectors f and x0
                n0 = exp**n0
                trbar0 = np.dot(dtr,x0.T)
                w0 = wage(k0,n0,z0)
                r0 = interest(k0,n0,z0)
                pen0=replacement_ratio*(1-taulbar)*lbar0*w0
        
                for ia in range(na):
                    a0=agrid[ia]
                    for iperm in range(nperm):
                        if pentype==0:
                            pene=pen
                        else:
                            pene = pen*perm[iperm]
                
                        for itheta in range(ntheta):
                            args1 = (a0,iage,iperm,itheta,pene,trbar0,r0,k1,iz,cpolicy)
                            x1=rfold(0,args1)
                                                  
                            if x1>0:     # corner solution?
                                aopt0=0
                            else:
                                amax2=pene+(1+(1-tauk)*r0)*a0+trbar0-psic   # maximum capital stock for c=0; */
                                x1=amax2
                                x1=min(x1,amax)
                                
                                if ia==0:
                                    agrid1 = np.linspace(0, 0.5, na)
                                else:
                                    x0=apolicy[iperm,itheta,ia-1,iage,ik,iz]
                                    
                                    agrid1 = np.linspace(x0,x1,na)
                                # search for initial value such that rfold >0
                                
                                ia0=-1
                                y0=-1
                                while ia0<na-1 and y0<0:
                                    ia0=ia0+1
                                    x0 = agrid1[ia0]
                                    y0 = rfold(x0,args1)
                                    
                                
                                if x0>=agrid[na-2]:
                                    if rfold(amax,args1)<=0:
                                        aopt0=amax # optimal next-period asset lies outside agrid
                                    else:
                                        aopt0 = scipy.optimize.fsolve(rfold1,x0,args=args1)
                                else:
                                    aopt0 = scipy.optimize.fsolve(rfold1,x0,args=args1)
                                   
                                if aopt0<amax and abs(rfold(aopt0,args1))>0.001:
                                    print("accuracy of solution in rfold:")
                                    print(abs(rfold(aopt0,args1)))
                                    print("itheta,iperm,ia,iage,ik,iz")
                                    print(itheta,iperm,ia,iage,ik,iz)
                                    x0 = input("Press Enter to continue: ")
                     
    
                     
    
                            c=pene+trbar0+(1+(1-tauk)*r0)*a0-ygrowth*aopt0
                            c=c/(1+tauc)
                    
                        
                            apolicy[iperm,itheta,ia,iage,ik,iz] = aopt0
                            cpolicy[iperm,itheta,ia,iage,ik,iz] = c
                            lpolicy[iperm,itheta,ia,iage,ik,iz] = 0    
    
                            # check for monotonocity of c(.) and a'(.)
                            if ia>0: 
                                if aopt0<apolicy[iperm,itheta,ia-1,iage,ik,iz]:
                                    print("apolicy not monotone")
                                    print(itheta,iperm,ia,iage,ik,iz)
                                    x0 = input("Press Enter to continue: ")
                                              
                                if c<cpolicy[iperm,itheta,ia-1,iage,ik,iz]:
                                    print("cpolicy not monotone")
                                    print(c,aopt0)
                                    print(itheta,iperm,ia,iage,ik,iz)
                                    x0 = input("Press Enter to continue: ")

    # computation of the decision rules for the worker 
    # age nw, nw-1,--,1 corresponding to variable iage=nw-1,..0
    
    for iage in range(nw-1,-1,-1):
        print(q,iage)
        for ik in range(nk):
            k0 = kgrid[ik]
            for iz in range(nz):
                z0 = zz[iz]
                x0 = np.array([1.0,iz,log(k0),iz*log(k0)])
                k1 = np.dot(d,x0.T)
                k1 = exp**k1
                n0 = np.dot(f,x0.T) # dot product of two vectors f and x0
                n0 = exp**n0
                trbar0 = np.dot(dtr,x0.T)
                w0 = wage(k0,n0,z0)
                r0 = interest(k0,n0,z0)
                pen0=replacement_ratio*(1-taulbar)*lbar0*w0
        
                for ia in range(na):
                    a0=agrid[ia]
                    for iperm in range(nperm):
                        if pentype==0:
                            pene=pen
                        else:
                            pene = pen*perm[iperm]
                
                        for itheta in range(ntheta):                            
                            args1 = (a0,iage,iperm,itheta,pene,trbar0,r0,w0,k1,iz,cpolicy)
                            x1=rfyoung(0,args1)
                          
                            weff = (1-taulbar)*perm[iperm]*theta[itheta]*ef[iage]*w0         
                            if x1>0:     # corner solution?
                                aopt0=0
                            else:
                                amax2=weff*lmax+(1+(1-tauk)*r0)*a0+trbar0-psic   # maximum capital stock for c=0; */
                                x1=min(amax2,amax)
                                
                                if ia==0:
                                    agrid1 = np.linspace(0, 0.5, na)
                                else:
                                    x0=apolicy[iperm,itheta,ia-1,iage,ik,iz]
                                    agrid1 = np.linspace(x0,x1,na)
                                # search for initial value such that rfold >0
                                ia0=-1
                                y0=-1
                                while ia0<=na-2 and y0<0:
                                    ia0=ia0+1
                                    x0 = agrid1[ia0]
                                    y0 = rfyoung(x0,args1)
                                    
                                
                                if x0>=agrid[na-2]:
                                    if rfyoung(amax,args1)<=0:
                                        aopt0=amax # optimal next-period asset lies outside agrid
                                    else:
                                        aopt0 = scipy.optimize.fsolve(rfyoung1,x0,args=args1)
                                else:
                                    aopt0 = scipy.optimize.fsolve(rfyoung1,x0,args=args1)
                                   
                                if aopt0<amax and abs(rfyoung(aopt0,args1))>0.001:
                                    print("accuracy of solution in rfyoung:")
                                    print(abs(rfyoung(aopt0,args1)))
                                    print("itheta,iperm,ia,iage,ik,iz")
                                    print(itheta,iperm,ia,iage,ik,iz)
                                    x0 = input("Press Enter to continue: ")
                     
    
                     
                            l0 = gamma-(1-gamma)*(trbar0+(1+(1-tauk)*r0)*a0-ygrowth*aopt0)/weff
    
                            if l0<0:
                                l0=0
                            elif l0>lmax:
                                l0=lmax
                                
                            c = weff*l0+trbar0+(1+(1-tauk)*r0)*a0-ygrowth*aopt0
                            c = c/(1+tauc)
                                
                            if iperm==0:
                                if itheta==0:
                                    if ia==0:
                                        if iage==0:
                                            if ik==3:
                                                if iz==1:
                                                    print(weff)
                                                    print(l0)
                                                    print(trbar0)
                                                    print(tauk)
                                                    print(aopt0)
                                                    print(c)
                                                    print(a0)
                                                    print(r0)
                        
                            apolicy[iperm,itheta,ia,iage,ik,iz] = aopt0
                            cpolicy[iperm,itheta,ia,iage,ik,iz] = c
                            lpolicy[iperm,itheta,ia,iage,ik,iz] = l0    
                          
    return apolicy, cpolicy, lpolicy


# rfyoung: 
# 
# purpose: computes optimal next-period asset of the young
#
# first order condition for young agent with l>0 
# input: next-period asset
# output: solution to Euler equation
def rfyoung(x,*args):
    a1 = x
    y = args[0]
    a0 = y[0]
    iage = y[1]
    iperm = y[2]
    itheta = y[3]
    pen0 = y[4]
    trbar0 = y[5]
    rbar0 = y[6]
    wbar0 = y[7]
    k1 = y[8]
    iz0 = y[9]
    cpolicy0 = y[10]
    
    # optimal labor from first-order condition w.r.t. labor 
    weff = (1-taulbar)*theta[itheta]*perm[iperm]*ef[iage]*wbar0 # effective wage
    l0 = gamma-(1-gamma)*(trbar0+(1+(1-tauk)*rbar0)*a0-ygrowth*a1)/weff
    
    if l0<0:
        l0=0
    elif l0>lmax:
        l0=lmax
    
    c0=weff*l0+trbar0+(1+(1-tauk)*rbar0)*a0-ygrowth*a1
    c0 = c0/(1+tauc)
    
    if c0<=psic:
        return (1+(c0-psic)**2)*1e5 # penalty function

    # Euler equation
    rf=uc(c0,l0)
    for iz1 in range(nz):
        zz1 = zz[iz1] 
        # next-period aggregate labor and transfers given expected K_1 and shock s[is1]
        x1 = np.array([1.0,iz1,log(k1),iz1*log(k1)])
        n1 = np.dot(f,x1.T) # dot product of two vectors f and x1 
        n1 = exp**n1 
        r1 = interest(k1,n1,zz1)    
        w1 = wage(k1,n1,zz1)
        factor=sp1[iage]*beta1*(1+(1-tauk)*r1)    
        for itheta1 in range(ntheta):
            # bi-linear interpolation of consumption at age iage+1
            ctemp = cpolicy0[iperm,itheta1,:,iage+1,:,iz1]
            c1 = bilint(agrid,kgrid,ctemp,a1,k1)
            # labor supply at age iage+1:
            if iage<nw-1:                
                weff1 = (1-taulbar)*theta[itheta1]*perm[iperm]*ef[iage+1]*w1 # effective wage
                l1=1-(1-gamma)/gamma*c1*(1+tauc)/weff1
                if l1<0:
                    l1=0
                elif l1>lmax:
                    l1=lmax
            else:
                l1=0
            
            if c1<=psic:
                rf = rf - factor*probth[itheta,itheta1]*probz[iz0,iz1]*(1+(c1-psic)**2)*1e5
            else:
                rf = rf - factor*probth[itheta,itheta1]*probz[iz0,iz1]*uc(c1,l1)
    
    return rf
 
# linear interpolation routine 
# that also extrapolates
def lininter(xd,yd,x,nx):
    if x <= xd[0]:
        y = yd[0] + (x-xd[0])*(yd[1]-yd[0])/(xd[1]-xd[0]) # extrapolation below xd[0]                
            
    elif x >= xd[nx-1]:
        y = yd[na-1] + (yd[na-1]-yd[na-2])*(x-xd[na-1])/(xd[na-1]-xd[na-2])
    else:
        j = sum(xd<=x) # x lies between xd[j-1] and xd[j]
        y = yd[j-1]+(yd[j]-yd[j-1])*(x-xd[j-1])/(xd[j]-xd[j-1])
    return y


# rfyoung: same as rfoungss, but as input into fsolve() => handling of tuple 'args'
# 
# purpose: computes optimal next-period asset of the young
#
# first order condition for young agent with l>0 
# input: next-period asset
# output: solution to Euler equation
def rfyoung1(x,*args):
    a1 = x
    a0 = args[0]
    iage = args[1]
    iperm = args[2]
    itheta = args[3]
    pen0 = args[4]
    trbar0 = args[5]
    rbar0 = args[6]
    wbar0 = args[7]
    k1 = args[8]
    iz0 = args[9]
    cpolicy0 = args[10]
 
    
 
    
    # optimal labor from first-order condition w.r.t. labor 
    weff = (1-taulbar)*theta[itheta]*perm[iperm]*ef[iage]*wbar0 # effective wage
    l0 = gamma-(1-gamma)*(trbar0+(1+(1-tauk)*rbar0)*a0-ygrowth*a1)/weff
    
    if l0<0:
        l0=0
    elif l0>lmax:
        l0=lmax
    
    
    c0=weff*l0+trbar0+(1+(1-tauk)*rbar0)*a0-ygrowth*a1
    c0 = c0/(1+tauc)
    
    if c0<=0:
        return (1+c0**2)*1e5 # penalty function

    # Euler equation
    rf=uc(c0,l0)
    for iz1 in range(nz):
        zz1 = zz[iz1] 
        # next-period aggregate labor and transfers given expected K_1 and shock s[is1]
        x1 = np.array([1.0,iz1,log(k1),iz1*log(k1)])
        n1 = np.dot(f,x1.T) # dot product of two vectors f and x1 
        n1 = exp**n1 
        r1 = interest(k1,n1,zz1)    
        w1 = wage(k1,n1,zz1)
        factor=sp1[iage]*beta1*(1+(1-tauk)*r1)    
        for itheta1 in range(ntheta):
            # bi-linear interpolation of consumption at age iage+1
            ctemp = cpolicy0[iperm,itheta1,:,iage+1,:,iz1]
            c1 = bilint(agrid,kgrid,ctemp,a1,k1)
            # labor supply at age iage+1:
            if iage<nw-1:                
                weff1 = (1-taulbar)*theta[itheta1]*perm[iperm]*ef[iage+1]*w1 # effective wage
                l1=1-(1-gamma)/gamma*c1*(1+tauc)/weff1
                if l1<0:
                    l1=0
                elif l1>lmax:
                    l1=lmax
            else:
                l1=0
            
            if c1<=psic:
                rf = rf - factor*probth[itheta,itheta1]*probz[iz0,iz1]*(1+(c1-psic)**2)*1e5
            else:
                rf = rf - factor*probth[itheta,itheta1]*probz[iz0,iz1]*uc(c1,l1)

    return rf

# rfold()
# first-order condition for retired agent
# returns residual function (=0 in the optimum)
# input: a' (next-period asset)
# output: residual of the Euler equation
def rfold(x,*args):
    a1 = x
    y = args[0]
    a0 = y[0]
    iage = y[1]
    iperm = y[2]
    itheta = y[3]
    pen0 = y[4]
    trbar0 = y[5]
    rbar0 = y[6]
    k1 = y[7]
    iz = y[8]
    cpolicy0 = y[9]
    
    
    c0 = trbar0+pen0+(1+(1-tauk)*rbar0)*a0-ygrowth*a1
    c0=c0/(1+tauc)
    
    if c0<=psic:
        return (1+(c0-psic)**2)*1e5 # penalty function 
    
    rf=uc(c0,0.0)
 
    # here, we assume that the idiosyncratic component z[iz] remains
    # stochastic during retirement --> innocuous assumption because
    # pension does not depend on z[iz]
    for iz1 in range(nz):
        zz1 = zz[iz1] 
        # next-period aggregate labor and transfers given expected K_1 and shock s[is1]
        x1 = np.array([1.0,iz1,log(k1),iz1*log(k1)])
        n1 = np.dot(f,x1.T) # dot product of two vectors f and x1 
        n1 = exp**n1 
        r1 = interest(k1,n1,zz1)    
        factor=sp1[iage]*beta1*(1+(1-tauk)*r1)    
        for itheta1 in range(ntheta):
            # bi-linear interpolation of consumption at age iage+1
            ctemp = cpolicy0[iperm,itheta1,:,iage+1,:,iz1]
            c1 = bilint(agrid,kgrid,ctemp,a1,k1)
            if c1<=psic:
                rf = rf - factor*probth[itheta,itheta1]*probz[iz,iz1]*(1+(c1-psic)**2)*1e5
            else:
                rf = rf - factor*probth[itheta,itheta1]*probz[iz,iz1]*uc(c1,0.0)
    
    return rf

# rfold1()
# first-order condition for retired agent
# returns residual function (=0 in the optimum)
# input: a' (next-period asset)
# output: residual of the Euler equation
def rfold1(x,*args):
    a1 = x
    a0 = args[0]
    iage = args[1]
    iperm = args[2]
    itheta = args[3]
    pen0 = args[4]
    trbar0 = args[5]
    rbar0 = args[6]
    k1 = args[7]
    iz = args[8]
    cpolicy0 = args[9]
    
    
    c0 = trbar0+pen0+(1+(1-tauk)*rbar0)*a0-ygrowth*a1
    c0=c0/(1+tauc)
#    print("iperm: " +str(iperm))
#    print("itheta: " +str(itheta))
#    print("a0: " +str(a0))
#    print("a1: " +str(a1))
#    print("c0: " +str(c0))
    if c0<=psic:
        return (1+(c0-psic)**2)*1e5 # penalty function 
    
    rf=uc(c0,0.0)
 
    # here, we assume that the idiosyncratic component z[iz] remains
    # stochastic during retirement --> innocuous assumption because
    # pension does not depend on z[iz]
    for iz1 in range(nz):
        zz1 = zz[iz1] 
        # next-period aggregate labor and transfers given expected K_1 and shock s[is1]
        x1 = np.array([1.0,iz1,log(k1),iz1*log(k1)])
        n1 = np.dot(f,x1.T) # dot product of two vectors f and x1 
        n1 = exp**n1 
        r1 = interest(k1,n1,zz1)    
        factor=sp1[iage]*beta1*(1+(1-tauk)*r1)    
        for itheta1 in range(ntheta):
            # bi-linear interpolation of consumption at age iage+1
            ctemp = cpolicy0[iperm,itheta1,:,iage+1,:,iz1]
            c1 = bilint(agrid,kgrid,ctemp,a1,k1)
            if c1<=psic:
                rf = rf - factor*probth[itheta,itheta1]*probz[iz,iz1]*(1+(c1-psic)**2)*1e5
            else:
                rf = rf - factor*probth[itheta,itheta1]*probz[iz,iz1]*uc(c1,0.0)
            
    
    return rf




# bi-linear interpolation: if x or y are outside the grid, extrapolation 
# input:
# - xd: column vector: dimension na 
# - yd: column vector: dimension nk
# - zd: matrix with rows xd and columns yd: (na x nk)
# - (x,y): interpolation points
# output: bi-linear interpolation
def bilint(xd,yd,zd,x,y):
    i1 = sum(xd<=x)
    i2 = sum(yd<=y)
    # if x/y is outside the vector agrid/kgrid
    # set i1/i2 to the corner point 
    i1 = max(1,i1)
    i2 = max(1,i2)
    i1 = min(i1,na-1)
    i2 = min(i2,nk-1)
    z1=(x-xd[i1-1])/(xd[i1]-xd[i1-1])
    z2=(y-yd[i2-1])/(yd[i2]-yd[i2-1])    
    
    z=(1-z1)*(1-z2)*zd[i1-1,i2-1]
    z=z+z1*(1-z2)*zd[i1,i2-1]
    z=z+z1*z2*zd[i1,i2]
    z=z+(1-z1)*z2*zd[i1-1,i2]
    return z




start_time = time.time()

# abbreviations
exp = np.e
log = math.log


# ---------------------------------------------------------------------------
#
# Step 1: Parameterization  
#
# -------------------------------------------------------------------------- 

# Step 1.1: Import data
data = pd.read_excel (r'survival_probs.xlsx') 
df = pd.DataFrame(data, columns= ['sp1','ef'])
print(df)
arr = np.array(df)
print(arr)
sp1 = arr[:,0]
ef = arr[0:45,1]
print(ef)


pentype = 1         # 0 -- lump-sum pensions for all
                    # 1 -- proporational to efficiency types

firstrun = 0        # 0 -- normal run, 1 -- loads aopt, copt, lopt during
                    # coding process

# Step 1.2: Parameterization of the model
# demographics
nage=70                # maximum age                  
nw=45                  # number of working years        
Rage=46                # first period of retirement 
nr=nage-Rage+1         # number of retirement years
popgrowth = 0.0075400000 # population growth rate

# preferences
beta1=1.011             # discount factor 
gamma=0.29		        # weight of consumption in utility
lbar=0.30		        # steady-state labor supply
eta=2.0		           # 1/IES

# production
ygrowth=1.02		    # annual growth factor
alpha=0.35             # production elasticity of capital
delta=0.083            # depreciation rate
rbbar=1.04		        # initial annual real interest rate on bonds

# fiscal policy and social security
taulbar=0.28	        # both taun+taup=0.28!!, see Mendoza, Razin (1994), p. 311
tauk=0.36              # capital income tax rate
tauc=0.05              # consumption tax rate
taup=0.124             # initial guess for social security contribution rate
replacement_ratio=0.494	# net replacement ratio
gybar=0.18					# government consumption-output ratio


#  Step 1.3: computational parameters 
kind_tol = 'linear'     # interpolation of policy functions:
                        # 'linear' or 'Cubic'
nq = 30        # number of outer iterations
psic=0.00001        # small constant 
lmax=0.9            # maximum labor
maxit=100           # maximum number of iterations over steady state k */
tolk=1e-5        # tolerance for equilibrium capital stock */
tol = 1e-5      # tolerance interpolation
psi0=0.7            # updating parameter: dynamics of K', L, Tr
kritold=100

# individual asset grid 
amin=0.0
amax=10.0
na=50     # grid on assets for value function 
agrid = np.linspace(amin, amax, na)   # asset grid policy function
# net income: net wage income plus net interest income 
ynetmin=0
ynetmax=2
ny=500
ynetgrid = np.linspace(ynetmin, ynetmax, ny)   # asset grid net income
# total income: wage income plus interest income plus pensions 
ytotmin=0
ytotmax=2
ytotgrid = np.linspace(ytotmin, ytotmax, ny)   # asset grid gross income
gtot = np.zeros(na)     # distribution total income



# productivity of workers
nperm=2;        # number of permanent productivities
perm = np.zeros(2)
perm[0]=0.57
perm[1]=1.43

ntheta=2            # number of stochastic prod types
varth=0.08;    # variance of individual productivity
probth = np.array([[0.9800,0.0200],[0.0200,0.9800]]) # transition matrix
theta = np.zeros(2)                                            
# compute productivities with function prod1(theta)
thetainitial = [0.8,1.2]
theta = scipy.optimize.fsolve(prod1, thetainitial, args=varth)
print('Idiosyncratic productivities')
print(theta)


# aggregate shock Z
# a period is calibrated as one year; the expected duration of 
# a business cycle is equal to six years 
nz = 2
zz = [0.98,1.02]
probz=np.array([[2/3,1/3],[1/3,2/3]]) # transition matrix aggregate prod Z

# Step 1.4
# initialization of aggregate value
bigl=0.2        # aggregate effective labor supply 
lbar=0.3        # average hourly labor supply 
lbarold=lbar
lbar0 = 0.30
rbar=0.02
kbar=(alpha/(delta+rbar))**(1/(1-alpha))
print("kbar initializiation: " +str(kbar))
bigk=kbar*bigl
wbar=wage(bigk,bigl,1.0)
print("wbar: " + str(wbar))
trbar=0         # initializiation transfers
ybar=production(bigk,bigl,1.0)
pen=replacement_ratio*(1-taulbar)*lbar*wbar
print("rbar: " +str(interest(bigk, bigl, 1.0)))


# Step 1.5
# computation of cohort measure for the model with nage periods
mass = np.ones(nage)
for iage in range(nage-1):
	mass[iage+1]=mass[iage]*sp1[iage]/(1+popgrowth)

mass=mass/sum(mass)   # normalization of the total cohort mass to one


# ---------------------------------------------------------------------------
#
# Step 2: computation of the steady state 
#
# -------------------------------------------------------------------------- 

# initialization policy functions in steady state 
gwages = np.zeros((nw,nperm,ny,2)) 
ass = np.zeros((nperm,ntheta,na,nage))  # next-period wealth at age iage,
                                        # wealth ia, permanent prod iperm, 
                                        # stochastic prod itheta
css = np.zeros((nperm,ntheta,na,nage))  # consumption
lss = np.zeros((nperm,ntheta,na,nage))  # labor supply
ga = np.zeros((nperm,ntheta,na,nage))
bequests=0
gbar=0
bigkold=0


    
# load results
ass = np.load('ass.npy')
css = np.load('css.npy')
lss = np.load('lss.npy')
ga = np.load('ga.npy')
lbar0 = np.load('lbar0.npy')
bequests = np.load('bequests.npy')
gbar = np.load('gbar.npy')
bigk = np.load('bigk.npy')
bigl = np.load('bigl.npy')
trbar = np.load('trbar.npy')


# -------------------------------------------------------------------------------
#
# Step 3: Computation of the dynamics with algorithm of Krussell/Smith
#
# ------------------------------------------------------------------------------ 

# Step 3.1: choice of grids for aggregate variables
#           initialization of individual policy functions
#           as functions of aggregate capital K, aggregate technology Z, 
#           and individual variables eps, z and age s



# aggregate capital stock grid
kmin=0.8*bigk
kmax=1.2*bigk
nk=7     # grid for individual policy functions
kgrid = np.linspace(kmin, kmax, nk)

# simulation parameter 
ndiscard=20   # discard the first ndiscard periods 
nsim=180+ndiscard  # number of periods 
nsim1 = 280+ndiscard # number of periods after the first q=4 simulations
#nsim = 30 + ndiscard
psi1=0.9   # update of d,f, dtr 
psi2 = 0.8  # update after the first q=4 simulations
nq=30      # maximum number of iterations over dynamics K', N' 

tolbeta1=0.001      # criterion to stop iteration over dynamics
kritbeta1=1+tolbeta1

# initialization:
# policy functions over business cycle 
aopt = np.zeros((nperm,ntheta,na,nage,nk,nz))    # optimal next period wealth a'(K,z,eps,a,age)
copt = np.zeros((nperm,ntheta,na,nage,nk,nz))    # optimal consumption c
lopt = np.zeros((nperm,ntheta,na,nage,nk,nz))   # optimal labor supply l



# -------------------------------------------------------------------
#
# Step 3.2.: initialization of the dynamics for K' and N'
#
# log(K') = d' * ( 1 z==1 log(K) (z==1)*log(K) )'
#
# log(N) = f' * (1 z==1 log(K) (z==1)*log(K) )'
#
#  tr is close to zero, therefore:
# tr = dtr' * (1 z==1 K (z==1)*K )'
#
# see: Carroll, Porapokkarm and Young
# 
# new: transfers are a function of K and Z
#
#---------------------------------------------------------------------

d = np.zeros(4)
d[2]=0.9
d[0]=(1.0-d[2])*log(bigk)


f=np.zeros(4)
f[2]=0
f[0]=log(bigl)

dtr = np.zeros(4)
dtr[2]=0
dtr[0]=trbar


# initialization: coefficients of OLS regression of K', L and tr on K and Z
beta1nt = np.zeros((nq,4))
beta1kt = np.zeros((nq,4))
beta1trt = np.zeros((nq,4))


# -------------------------------------------------------------------------
#
# Step 3.3: assigning steady-state policy function values to
#          the optimal policy functions during the transition as initial guess
#
# -------------------------------------------------------------------------- 

for ik in range(nk):
    for iz in range(nz):
        aopt[:,:,:,:,ik,iz] = ass
        copt[:,:,:,:,ik,iz] = css
        lopt[:,:,:,:,ik,iz] = lss
        


q=-1
crit = 1+tol
# loop over dynamics d, f, dtr
while q<nq-1 and kritbeta1>tolbeta1:
    q=q+1
    print(q)
    
    if q==4:    # increase of simulation periods and update parameter 1-psi1
        nsim = nsim1
        psi1 = psi2
    
    # firstrun=1 is only used during the coding process
    # to test the loop over the dynamics d, f, dtr
    if firstrun==0:
        aopt,copt,lopt = getpolicy()        
        np.save('aopt',aopt)
        np.save('copt',copt)   
        np.save('lopt',lopt)
    else:
        aopt = np.load('aopt.npy')
        copt = np.load('copt.npy')
        lopt = np.load('lopt.npy')
        

    zshock = np.random.uniform(0.0,1.0,nsim)
    
    # time series of aggregate variables K, N, ...
    kt = np.zeros(nsim)     # predicted value of K_t in period t-1
    nt = np.zeros(nsim)     # predicted value of N_t at beginning of period t
    trt = np.zeros(nsim)    # predicted Tr_t at beginning of period t
    
    ztsim = np.zeros(nsim)  # realization of Z_t in period t 
    ktsim = np.zeros(nsim)  # simulated values 
    ntsim = np.zeros(nsim)  # aggregate labor
    ytsim = np.zeros(nsim)  # aggregate production
    trtsim= np.zeros(nsim)  # government transfers
    beqtsim = np.zeros(nsim)    # accidental bequests
    giniynetsim = np.zeros(nsim) # Gini coefficient net income
    giniytotsim = np.zeros(nsim) # Gini gross income
    sharey20sim = np.zeros(nsim) # income share bottom quintile
    sharey40sim = np.zeros(nsim)
    sharey60sim = np.zeros(nsim)
    sharey80sim = np.zeros(nsim)
    sharey95sim = np.zeros(nsim)
    sharey100sim = np.zeros(nsim)


    kt[0]=bigk 
    k0=bigk
    ktsim[0]=bigk
    nt[0]=bigl
    ntsim[0]=bigl
    trt[0]=trbar 
    trtsim[0]=trbar
    beqtsim[0] = bequests
    
    if zshock[0]<0.5:
        ztsim[0]=0 
        iz=0
    else:
        ztsim[0]=1 
        iz=1
    

    ga0=ga  # distribution in the first period: steady state distribution 

    isim=0 # loop over simulation periods
    while isim<nsim-1 and ((k0>kmin) and (k0<kmax)):
        isim = isim+1
        print(q,isim,k0)
        
        gynet = np.zeros(ny)
        gytot = np.zeros(ny)
        gwealth = np.zeros(na)
        k0 = ktsim[isim-1]
        iz = ztsim[isim-1]
        iz = iz.astype(int) # conversion to type integer so that 
                            # it can be used as an argument of copt, lopt, aopt

        # computation of aggregate labor and consumption in period isim
        # given optimal policy functions l(.) and c(.) for exogenous 
        # aggregates K=k0 and Z=zz[iz1]
        nsum = 0           # aggregate labor
        bigc0 = 0           # aggregate consumption
        for iage in range(nage):
            for ia in range(na):
                for iperm in range(nperm):
                    for itheta in range(ntheta):
                        mass0 = ga0[iperm,itheta,ia,iage]
                        if iage<=nw-1:
                            labor = lininter(kgrid,lopt[iperm,itheta,ia,iage,:,iz],k0,nk)
                            nsum = nsum + mass0 * labor*ef[iage]*perm[iperm]*theta[itheta]
                        
                        else:
                            labor = 0
                        
                        # optimal consumption in period isim
                        # with individual states (iperm,itheta,ia,iage)
                        # and aggregate states K=k0 and iz
                        c0 = lininter(kgrid,copt[iperm,itheta,ia,iage,:,iz],k0,nk)
                        bigc0 = bigc0 + mass0*c0
             
        print("bigl~nsum~bigc: ")
        print(bigl,nsum,bigc0)
        ntsim[isim-1]=nsum;
        
        
        #
        # prediction of k' given Z=s[is] and K=k0
        x0 = np.array([1.0,iz,log(k0),iz*log(k0)])
        k1 = np.dot(d,x0.T)
        k1 = exp**k1
        kt[isim] = k1
        
        if zshock[isim]<probz[iz,0]:
            iz1 = 0
        else:
            iz1 = 1
        
        ztsim[isim] = iz1
        print(iz1)
        
        # predicted aggregate labor L and transfers tr
        x1 = np.array([1.0,iz1,log(k1),iz1*log(k1)])
        n1 = np.dot(f,x1.T) # dot product of two vectors f and x0
        n1 = exp**n1
        nt[isim]=n1        # predicted labor in period isim
        n0 = nt[isim-1]    # actual (simulated) labor in period isim-1
    

        r0 = interest(k0,n0,zz[iz])
        w0 = wage(k0,n0,zz[iz])
        pen0 = replacement_ratio*(1-taulbar)*lbar0*w0

        # computation of aggregate transfers
        totalpen0 = sum(mass[nw:nage])*pen0
        # consolidated budget: fiscal + social security
        trbar0 = beqtsim[isim-1]+taulbar*w0*n0+tauk*r0*k0+tauc*bigc0-gbar-totalpen0
        trtsim[isim-1] = trbar0
        trbar1 =  np.dot(dtr,x1.T)
        trt[isim] = trbar1    # predicted transfers trbar in period isim

        #
        # dynamics of the distribution
        # 
        ga1 = np.zeros((nperm,ntheta,na,nage)) # next-period distribution: period isim 
        bigasum = 0          # aggregate wealth
        masstest = 0

        # initialization: assets at age 1 are 0 
        for iperm in range(nperm):
            for itheta in range(ntheta):
                # assumptions: all productivities have the same measure 
                ga1[iperm,itheta,0,0] = mass[0]/(nperm*ntheta)

        # compute the wealth in period isim of the households 
        # aged iage=1, .., nage-1
        # with wealth agrid[ia], permanent productivity perm[iperm] and 
        # idiosyncratic productivity theta[itheta]
        #        
        for iage in range(nage-1):
            # print(iage)
            for ia in range(na):
                for iperm in range(nperm):
                    if pentype==0:
                        pene=pen
                    else:
                        pene = pen0*perm[iperm]
                        
                    for itheta in range(ntheta):
                        # distribution in period isim-1
                        mass0 = ga0[iperm,itheta,ia,iage] 
                        # wealth distribution in period isim-1
                        gwealth[ia] = gwealth[ia]+mass0        
                        masstest = masstest + mass0
                        
                        # share of those who survive until the next period
                        mass1 = mass0*sp1[iage]/(1+popgrowth);   
                        # interpolation of optimal savings at K=k0
                        a1 = lininter(kgrid,aopt[iperm,itheta,ia,iage,:,iz],k0,nk)
                        # aggregateion accidental bequests
                        beqtsim[isim] = beqtsim[isim] + mass0*(1-sp1[iage])*a1*ygrowth
                        ynet0 = (1-tauk)*r0*agrid[ia]    # net income in period isim-1
                        ytot0 = r0*agrid[ia]+trbar0      # gross income in period isim-1
                        
                        if iage<=nw-1:   # labor supply in period isim-1 
                            labor =  lininter(kgrid,lopt[iperm,itheta,ia,iage,:,iz],k0,nk)
                            ynet0 = ynet0+(1-taulbar)*labor*ef[iage]*perm[iperm]*theta[itheta]*w0
                            ytot0 = ytot0+labor*ef[iage]*perm[iperm]*theta[itheta]*w0
                        else:
                            ynet0 = ynet0 + pene
                            ytot0 = ytot0 +pene

                        # distribution gross and net income in period isim-1
                        iy1 = sum(ynetgrid<ynet0)          
                        if iy1==0:
                            gynet[0]=gynet[0]+mass0
                        elif iy1==ny:
                            gynet[na-1]=gynet[na-1]+mass0
                        else:
                            iy0=(ynet0-ynetgrid[iy1-1])/(ynetgrid[iy1]-ynetgrid[iy1-1])
                            gynet[iy1-1]=gynet[iy1-1]+(1-iy0)*mass0
                            gynet[iy1]=gynet[iy1]+iy0*mass0
                    
                
                        iy1 = sum(ytotgrid<ytot0)
                        if iy1==0:
                            gytot[0]=gytot[0]+mass0
                        elif iy1==ny:
                            gytot[na-1]=gytot[na-1]+mass0
                        else:
                            iy0=(ytot0-ytotgrid[iy1-1])/(ytotgrid[iy1]-ytotgrid[iy1-1])
                            gytot[iy1-1]=gytot[iy1-1]+(1-iy0)*mass0
                            gytot[iy1]=gytot[iy1]+iy0*mass0
                

                        # aggregate capital in period isim
                        bigasum=bigasum+a1*mass1 
                        
                        if a1<=amin:
                            ia1=1
                            iz1=0
                            for itheta1 in range(ntheta):
                                ga1[iperm,itheta1,0,iage+1]=ga1[iperm,itheta1,0,iage+1]+probth[itheta,itheta1]*mass1    
                           
                        elif a1>=amax:
                            ia1=na-1
                            for itheta1 in range(ntheta):
                                ga1[iperm,itheta1,ia1,iage+1]=ga1[iperm,itheta1,ia1,iage+1]+probth[itheta,itheta1]*mass1    
                           
                        else:   # a1 lies between agrid[0] and agrid[na-1]
                            ia1=sum(agrid<a1) # a1 lies between agrid[ia-2] and agrid[ia-1]
                            is1=(a1-agrid[ia1-1])/(agrid[ia1]-agrid[ia1-1])
                            for itheta1 in range(ntheta):
                                ga1[iperm,itheta1,ia1,iage+1] = ga1[iperm,itheta1,ia1,iage+1]+is1*probth[itheta,itheta1]*mass1
                                ga1[iperm,itheta1,ia1-1,iage+1] = ga1[iperm,itheta1,ia1-1,iage+1]+(1-is1)*probth[itheta,itheta1]*mass1
                        
    
        iage=nage-1
        print(iage)
        for ia in range(na):
            for iperm in range(nperm):
                if pentype==0:
                    pene=pen0
                else:
                    pene = pen0*perm[iperm]
                        
                for itheta in range(ntheta):
                    # distribution in period isim-1
                    mass0 = ga0[iperm,itheta,ia,iage] 
                    # wealth distribution in period isim-1
                    gwealth[ia] = gwealth[ia]+mass0    
                    masstest = masstest + mass0
                    
                    a1=0   # no savings in the final period of life
                    # net and gross income in the last period of retirement
                    ynet0 = (1-tauk)*r0*agrid[ia]+pene
                    ytot0 = r0*agrid[ia]+trbar0+pene

                    # distribution gross and net income in period isim-1
                    iy1 = sum(ynetgrid<ynet0)          
                    if iy1==0:
                        gynet[0]=gynet[0]+mass0
                    elif iy1==ny:
                        gynet[na-1]=gynet[na-1]+mass0
                    else:
                        iy0=(ynet0-ynetgrid[iy1-1])/(ynetgrid[iy1]-ynetgrid[iy1-1])
                        gynet[iy1-1]=gynet[iy1-1]+(1-iy0)*mass0
                        gynet[iy1]=gynet[iy1]+iy0*mass0
                    
                
                    iy1 = sum(ytotgrid<ytot0)
                    if iy1==0:
                        gytot[0]=gytot[0]+mass0
                    elif iy1==ny:
                        gytot[na-1]=gytot[na-1]+mass0
                    else:
                        iy0=(ytot0-ytotgrid[iy1-1])/(ytotgrid[iy1]-ytotgrid[iy1-1])
                        gytot[iy1-1]=gytot[iy1-1]+(1-iy0)*mass0
                        gytot[iy1]=gytot[iy1]+iy0*mass0
           
        ktsim[isim] = bigasum
        ytsim[isim-1] = zz[iz]*ktsim[isim-1]**alpha*ntsim[isim-1]**(1-alpha)
        print("q~isim: ")
        print(q,isim)
        print("K'~N~bigk~bigl~tr: ")
        print(bigasum,nsum,bigk,bigl,trbar0)
        # Gini coefficients of distributions
        giniwealth=ginid(agrid,gwealth,na)
        print("Gini wealth: " +str(giniwealth))
     
                    
        y11=ginid(ynetgrid,gynet,ny)
        giniynetsim[isim-1]=y11
        y11=ginid(ytotgrid,gytot,ny)
        giniytotsim[isim-1]=y11
        
        gtotsum = np.cumsum(gytot)      
        totalincome = np.dot(gytot,ytotgrid) # total income
        gtotsum1 = gytot*ytotgrid
        gtotsum1 = np.cumsum(gtotsum1)
        nshare20 = sum(gtotsum<=0.20)-1 # first element in vector has index "0" 
        sharey20sim[isim-1]=gtotsum1[nshare20]/gtotsum1[ny-1]
        nshare40 = sum(gtotsum<=0.4)-1 
        sharey40sim[isim-1] = (gtotsum1[nshare40]-gtotsum1[nshare20])/gtotsum1[ny-1]
        nshare60 = sum(gtotsum<=0.6)-1
        sharey60sim[isim-1] = (gtotsum1[nshare60]-gtotsum1[nshare40])/gtotsum1[ny-1]
        nshare80 = sum(gtotsum<=0.8)-1
        sharey80sim[isim-1]=(gtotsum1[nshare80]-gtotsum1[nshare60])/gtotsum1[ny-1]
        nshare95 = sum(gtotsum<=0.95)-1
        sharey95sim[isim-1] = (gtotsum1[nshare95]-gtotsum1[nshare80])/gtotsum1[ny-1]
        sharey100sim[isim-1] = 1-gtotsum1[nshare95]/gtotsum1[ny-1]
        ga0=ga1


    # correlation of income shares with detrended output
    
    yt1 = np.log(ytsim[ndiscard:nsim-1])
    incomet1 = giniytotsim[ndiscard:nsim-1]
    incomet2 = sharey20sim[ndiscard:nsim-1]
    incomet3 = sharey40sim[ndiscard:nsim-1]
    incomet4 = sharey60sim[ndiscard:nsim-1]
    incomet5 = sharey80sim[ndiscard:nsim-1]
    incomet6 = sharey95sim[ndiscard:nsim-1]
    incomet7 = sharey100sim[ndiscard:nsim-1]
    cycle, trend = sm.tsa.filters.hpfilter(yt1, 100) # annual data, HP filter with 100
    
    print("standard deviation of y: " +str(np.std(cycle)))
    matrix = np.vstack((cycle,incomet1))
    matrix = np.vstack((matrix,incomet2))
    matrix = np.vstack((matrix,incomet3))
    matrix = np.vstack((matrix,incomet4))
    matrix = np.vstack((matrix,incomet5))
    matrix = np.vstack((matrix,incomet6))
    matrix = np.vstack((matrix,incomet7))
    corin = np.corrcoef(matrix)
    print(corin[:,0]) # correlations with cyclical component of income

    # update of dynamics for predicted K', N, tr
    #
    # update of coefficients d, f, dtr
    #
    
    # setting up the regressor xk
    xk = np.zeros((nsim-ndiscard-2,4))
    xk0 = np.ones(nsim-ndiscard-2)
    xk1= ztsim[ndiscard:nsim-2]
    xk2 = ktsim[ndiscard:nsim-2]
    xk2 = np.log(xk2)
    xk3 = xk1*xk2    
    xk[0:nsim-ndiscard-2,0] = xk0
    xk[0:nsim-ndiscard-2,1] = xk1
    xk[0:nsim-ndiscard-2,2] = xk2
    xk[0:nsim-ndiscard-2,3] = xk3
    # regression of ln (K') 
    yk = ktsim[ndiscard+1:nsim-1]
    yk = np.log(yk)
    betak_hat = np.linalg.solve(xk.T @ xk, xk.T @ yk)
    print("d: ")
    print(d)
    print("regression beta: ")
    print(betak_hat)

    # regresssion of N on [1~ln(K)]
    yn = ntsim[ndiscard:nsim-2]
    yn = np.log(yn)
    betan_hat = np.linalg.solve(xk.T @ xk, xk.T @ yn)
    print("f: ")
    print(f)
    print("regression beta: ")
    print(betan_hat)
    
    kritbetak1 = max(abs(d-betak_hat))
    kritbetan1 = max(abs(f-betan_hat))
    kritbeta1 = max(kritbetak1,kritbetan1)

    # regression of transfers on K and Z, not on ln K 
    xk_tr = np.zeros((nsim-ndiscard-2,4))
    xk0 = np.ones(nsim-ndiscard-2)
    xk1= ztsim[ndiscard:nsim-2]
    xk2 = ktsim[ndiscard:nsim-2]
    xk3 = xk1*xk2    
    xk_tr[0:nsim-ndiscard-2,0] = xk0
    xk_tr[0:nsim-ndiscard-2,1] = xk1
    xk_tr[0:nsim-ndiscard-2,2] = xk2
    xk_tr[0:nsim-ndiscard-2,3] = xk3
    # regression of tr on 1~K
    ytr = trtsim[ndiscard:nsim-2]
    betatr_hat = np.linalg.solve(xk_tr.T @ xk_tr, xk_tr.T @ ytr)
    print("dtr: ")
    print(dtr)
    print("regression beta: ")
    print(betatr_hat)
    
    beta1nt[q,:] = betan_hat
    beta1kt[q,:] = betak_hat
    beta1trt[q,:] = betatr_hat

    if q<3:
        # update d[2] and d[4] and choose d[1] and d[3]
        # so that mean capital stocks in simulation are equal 
        # to mean capital stock of dynamics */
        d[2]=psi1*d[2]+(1-psi1)*betak_hat[2]
        d[3]=psi1*d[3]+(1-psi1)*betak_hat[3]
        
        d[0] = np.log(bigk)*(1-d[2])
        d[1] = np.log(bigk)*(1-d[2]-d[3])-d[0]

    else:
        d = psi1*d + (1-psi1)*betatr_hat

    f = psi1*f+(1-psi1)*betan_hat
    dtr = psi1*dtr+(1-psi1)*betatr_hat

    print("new d,f,dtr")
    print(d,f,dtr)
    print("kritbeta1: " +str(kritbeta1))
    
    
    # print regression results at the end of the loop q
    if kritbeta1<tolbeta1 or q==nq-1:
        print("regression: employment")
        regn = sm.OLS(endog=yn, exog=xk, missing='drop')
        type(regn)
        resultsn = regn.fit()
        type(resultsn)
        print(resultsn.summary())
        
        print("regression: next-period capital")
        regk = sm.OLS(endog=yk, exog=xk, missing='drop')
        type(regk)
        resultsk = regk.fit()
        type(resultsk)
        print(resultsk.summary())

        print("regression: transfers")
        regtr = sm.OLS(endog=ytr, exog=xk_tr, missing='drop')
        type(regtr)
        resultstr = regtr.fit()
        type(resultstr)
        print(resultstr.summary())
        
        # plotting series over a sub-period        
        nsub = min(30,nsim-ndiscard-2)     
        
        plt.xlabel('period')
        plt.ylabel('production Y')
        plt.plot(range(1,nsub+1),ytsim[ndiscard+1:ndiscard+nsub+1])
        plt.show()          
  
        fig, ax = plt.subplots()
        label1 = 'simulated'
        label2 = 'predicted'
        plt.xlabel('period')
        plt.ylabel('capital stock K')
        ax.plot(range(1,nsub+1),ktsim[ndiscard+1:ndiscard+nsub+1], linewidth=2, label=label1)
        ax.plot(range(1,nsub+1),kt[ndiscard+1:ndiscard+nsub+1], linewidth=2, label=label2)
        ax.legend()
        plt.show()
        
          
        fig, ax = plt.subplots()
        label1 = 'simulated'
        label2 = 'predicted'
        plt.xlabel('period')
        plt.ylabel('employment L')
        ax.plot(range(1,nsub+1),ntsim[ndiscard+1:ndiscard+nsub+1], linewidth=2, label=label1)
        ax.plot(range(1,nsub+1),nt[ndiscard+1:ndiscard+nsub+1], linewidth=2, label=label2)
        ax.legend()
        plt.show()
        
        
          
        fig, ax = plt.subplots()
        label1 = 'simulated'
        label2 = 'predicted'
        plt.xlabel('period')
        plt.ylabel('transfers Tr')
        ax.plot(range(1,nsub+1),trtsim[ndiscard+1:ndiscard+nsub+1], linewidth=2, label=label1)
        ax.plot(range(1,nsub+1),trt[ndiscard+1:ndiscard+nsub+1], linewidth=2, label=label2)
        ax.legend()
        plt.show()
        
        
print("runtime: --- %s seconds ---" % (time.time() - start_time))    
sec = (time.time() - start_time)
ty_res = time.gmtime(sec)
res = time.strftime("%H : %M : %S", ty_res)
print(res)