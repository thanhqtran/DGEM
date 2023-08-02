# OLG_krusell_smith_ss.py
"""
Created on January 24, 2022

@author: Burkhard Heer

purpose: solves the stochastic Overlapping Generations Model 
        in Chapter 10.2.2 of Heer/Maussner, 
        'Dynamic General Equilibrium Modeling: Computational Methods
         and Applications', 3rd edition (scheduled for 2022)

        part 1: computes non-stochatsic steady state

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
def getpolicyss(pen,rbar,trbar,wbar):
    apolicy = np.zeros((nperm,ntheta,na,nage))  # next-period wealth at age iage,
                                        # wealth ia, permanent prod iperm, 
                                        # stochastic prod itheta
    cpolicy = np.zeros((nperm,ntheta,na,nage))  # consumption
    lpolicy = np.zeros((nperm,ntheta,na,nage))  # labor supply
     
    # last period nage of life
    for ia in range(na):
        a0 = agrid[ia]
        for iperm in range(nperm):
            if pentype==0:
                pene=pen
            else:
                pene = pen*perm[iperm]
                
            for itheta in range(ntheta):
                c= (1+(1-tauk)*rbar)*a0+pene+trbar 
                c = c/(1+tauc)
                apolicy[iperm,itheta,ia,nage-1] = 0
                cpolicy[iperm,itheta,ia,nage-1] = c
                lpolicy[iperm,itheta,ia,nage-1] = 0
            
    
    # computation of the decision rules for the retired 
    # age nage, nage-1,--,nw+1
    
    for iage in range(nage-2,nw-1,-1):
        print(q,iage)
        for ia in range(na):
            a0=agrid[ia]
            for iperm in range(nperm):
                if pentype==0:
                    pene=pen
                else:
                    pene = pen*perm[iperm]
                
                for itheta in range(ntheta):
                    args1 = (a0,iage,iperm,itheta,pene,trbar,rbar,cpolicy)
                    x1=rfoldss(0,args1)
                    if x1>0:     # corner solution? 
                        aopt0=0
                    else:
                        # find optimal starting value
                        # option 1: half the maximum possible a'
                        # option 2: a'=a
                        # option 3: a' for the same agent, one year older
                        # option 4: a' for the same agent, wealth one
                        #           gridpoint lower agrid[ia-1]
                        amax2=pene+trbar+(1+(1-tauk)*rbar)*a0-psic   # maximum capital stock for c=0; */
                        x0=amax2/2
                        y0=rfoldss(x0,args1)
                        y1=rfoldss(a0,args1)
                        if abs(y1)<abs(y0):
                            x0=a0
                            y0=y1
                        y1 = rfoldss(apolicy[iperm,itheta,ia,iage+1],args1)
                        if abs(y1)<abs(y0):
                            x0=apolicy[iperm,itheta,ia,iage+1]
                            y0=y1
                        if ia>0:
                            y1=rfoldss(apolicy[iperm,itheta,ia-1,iage],args1)
                            if abs(y1)<abs(y0):
                                x0=apolicy[iperm,itheta,ia-1,iage]
                        
                        aopt0 = scipy.optimize.fsolve(rfoldss1,x0,args=args1)
                        if abs(rfoldss(aopt0,args1))>0.001:
                            print("accuracy of solution in rfoldss")
                            print(abs(rfoldss(aopt0,args1)))
                            x0 = input("Press Enter to continue: ")
                     
                        
                    c=pene+trbar+(1+(1-tauk)*rbar)*a0-ygrowth*aopt0
                    c=c/(1+tauc)
                    
                        
                    apolicy[iperm,itheta,ia,iage] = aopt0
                    cpolicy[iperm,itheta,ia,iage] = c
                    lpolicy[iperm,itheta,ia,iage] = 0    


    # computation of the decision rules for the workers 
    # age nw, nw-1,.. 1 => index iage=nw-1,..,0
    for iage in range(nw-1,-1,-1):
        print(q,iage)
        for ia in range(na):
            a0=agrid[ia]
            for iperm in range(nperm):
                if iage==nw-1:  # next period retiree
                    if pentype==0:
                        pene=pen
                    else:
                        pene = pen*perm[iperm]
                
                for itheta in range(ntheta):
                    args1 = (a0,iage,iperm,itheta,pene,trbar,rbar,wbar,cpolicy)
                    x1=rfyoungss(0,args1)
                    weff=(1-taulbar)*theta[itheta]*perm[iperm]*ef[iage]*wbar 
                    
                    
                    if x1>0:     # corner solution for a'? 
                        aopt0=0
                    else:                       
                        # looking for a good initial value for a'
                        # 1. half the maximum possible a' for l=lmax
                        # 2. a'=a
                        # 3. optimal a' at age it+1
                        # 4. optimal a' at wealth a[ia-1] for ia>1
                        
                        amax2=(weff*lmax+trbar+(1+(1-tauk)*rbar)*a0-psic)/ygrowth   # c=0 
                        amax2 = amax2/2
                        x0=apolicy[iperm,itheta,ia,iage+1]
                        if abs(rfyoungss(amax2,args1))<abs(rfyoungss(x0,args1)):
                            x0=amax2
                        
                        if abs(rfyoungss(a0,args1))<abs(rfyoungss(x0,args1)):
                            x0 = a0
                            
                        if ia>1:
                            x1=apolicy[iperm,itheta,ia-1,iage+1]
                            if abs(rfyoungss(x1,args1))<abs(rfyoungss(x0,args1)):
                                x0=x1
                            
                        # find optimal next-period assets as solution to
                        # the Euler equation rfyoungss1
                        aopt0 = scipy.optimize.fsolve(rfyoungss1,x0,args=args1)   
                     

                    l = gamma-(1-gamma)*(trbar+(1+(1-tauk)*rbar)*a0-ygrowth*aopt0)/weff
                    if l<0:    # corner solution?
                        l=0
                    elif l>lmax:
                        l=lmax
                    
                    c = weff*l+trbar+(1+(1-tauk)*rbar)*a0-ygrowth*aopt0
                    c=c/(1+tauc)
                    
                    # check for monotonocity of a', l and c
                    # a' monoton increasing?
                    if ia>0:
                        if apolicy[iperm,itheta,ia-1,iage]>aopt0:
                            print("iage: " +str(iage))
                            print("ia: " +str(ia))
                            print("itheta: " +str(itheta))
                            print("iperm: " +str(iperm))
                            print("a' not monotone increasing in a")
                            x0 = input("Press Enter to continue: ")
                        if cpolicy[iperm,itheta,ia-1,iage]>c:
                            print("iage: " +str(iage))
                            print("ia: " +str(ia))
                            print("itheta: " +str(itheta))
                            print("iperm: " +str(iperm))
                            print("c not monotone increasing in a")
                            x0 = input("Press Enter to continue: ")
                        if lpolicy[iperm,itheta,ia-1,iage]<l:
                            print("iage: " +str(iage))
                            print("ia: " +str(ia))
                            print("itheta: " +str(itheta))
                            print("iperm: " +str(iperm))
                            print("l not monotone decreasing in a")
                            x0 = input("Press Enter to continue: ")
                            
                    apolicy[iperm,itheta,ia,iage] = aopt0
                    cpolicy[iperm,itheta,ia,iage] = c
                    lpolicy[iperm,itheta,ia,iage] = l    
                  
    return apolicy, cpolicy, lpolicy


# rfyoungss: 
# 
# purpose: computes optimal next-period asset of the young
#
# first order condition for young agent with l>0 
# input: next-period asset
# output: solution to Euler equation
def rfyoungss(x,*args):
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
    cpolicy0 = y[8]
    
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

    rf=uc(c0,l0)
    factor=sp1[iage]*beta1*(1+(1-tauk)*rbar0)
    
    
    for itheta1 in range(ntheta):
        # linear interpolation between two grid points of a' 
        # interpolation: next-period consumption
        
        # if a1<=amin:
        #    c1 = cpolicy0[iperm,itheta1,0,iage+1]
        # elif a1>=amax:
        #    c1 = cpolicy0[iperm,itheta1,na-1,iage+1]
        # else:
        #    c_polate = interpolate.interp1d(agrid,cpolicy0[iperm,itheta1,:,iage+1], kind=kind_tol)
        #    c1 = c_polate(a1)
         
        c1 = lininter(agrid,cpolicy0[iperm,itheta1,:,iage+1],a1)
        if iage<nw-1:   # next period: worker
            # optimal labor from first-order condition w.r.t. labor
            weff1=(1-taulbar)*theta[itheta1]*perm[iperm]*ef[iage+1]*wbar0
            l1=1-(1-gamma)/gamma*c1*(1+tauc)/weff1  # foc with respect to labor
            if l1<0:
                l1=0
            elif l1>lmax:
                l1=lmax
        else:
            l1=0
        
        rf = rf - factor*probth[itheta,itheta1]*uc(c1,l1)
    
    return rf
 
# linear interpolation routine
# that also extrapolates
def lininter(xd,yd,x):
    if x <= xd[0]:
        y = yd[0] + (x-xd[0])*(yd[1]-yd[0])/(xd[1]-xd[0]) # extrapolation below xd[0]                
            
    elif x >= xd[na-1]:
        y = yd[na-1] + (yd[na-1]-yd[na-2])*(x-xd[na-1])/(xd[na-1]-xd[na-2])
    else:
        j = sum(xd<=x) # x lies between xd[j-1] and xd[j]
        y = yd[j-1]+(yd[j]-yd[j-1])*(x-xd[j-1])/(xd[j]-xd[j-1])
    return y



# rfyoungss1: same as rfoungss, but as input into fsolve() => handling of tuple 'args'
# 
# purpose: computes optimal next-period asset of the young
#
# first order condition for young agent with l>0 
# input: next-period asset
# output: solution to Euler equation
def rfyoungss1(x,*args):
    a1 = x
    a0 = args[0]
    iage = args[1]
    iperm = args[2]
    itheta = args[3]
    pen0 = args[4]
    trbar0 = args[5]
    rbar0 = args[6]
    wbar0 = args[7]
    cpolicy0 = args[8]
    
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

    rf=uc(c0,l0)
    factor=sp1[iage]*beta1*(1+(1-tauk)*rbar0)
    
    
    for itheta1 in range(ntheta):
        # linear interpolation between two grid points of a' 
        # interpolation: next-period consumption
        # if a1<=amin:
        #     c1 = cpolicy0[iperm,itheta1,0,iage+1]
        # elif a1>=amax:
        #    c1 = cpolicy0[iperm,itheta1,na-1,iage+1]
        # else:
        #    c_polate = interpolate.interp1d(agrid,cpolicy0[iperm,itheta1,:,iage+1], kind=kind_tol)
        #    c1 = c_polate(a1)
        
        c1 = lininter(agrid,cpolicy0[iperm,itheta1,:,iage+1],a1)
        if iage<nw-1:   # next period: worker
            # optimal labor from first-order condition w.r.t. labor
            weff1=(1-taulbar)*theta[itheta1]*perm[iperm]*ef[iage+1]*wbar0
            l1=1-(1-gamma)/gamma*c1*(1+tauc)/weff1  # foc with respect to labor
            if l1<0:
                l1=0
            elif l1>lmax:
                l1=lmax
        else:
            l1=0
        
        rf = rf - factor*probth[itheta,itheta1]*uc(c1,l1)
    
    return rf

# rfoldss()
# first-order condition for retired agent: steady state 
# returns residual function (=0 in the optimum)
# input: k' (next-period asset)
# output: residual of the Euler equation
def rfoldss(x,*args):
    a1 = x
    y = args[0]
    a0 = y[0]
    iage = y[1]
    iperm = y[2]
    itheta = y[3]
    pen0 = y[4]
    trbar0 = y[5]
    rbar0 = y[6]
    cpolicy0 = y[7]
    
    c0 = trbar0+pen0+(1+(1-tauk)*rbar0)*a0-ygrowth*a1
    c0=c0/(1+tauc)
    
    if c0<=0:
        return (1+c0**2)*1e5 # penalty function 
    
    rf=uc(c0,0.0)
    factor=sp1[iage]*beta1*(1+(1-tauk)*rbar0)
 
    # here, we assume that the idiosyncratic component z[iz] remains
    # stochastic during retirement --> innocuous assumption because
    # pension does not depend on z[iz]
    for itheta1 in range(ntheta):
        # linear interpolation between two grid points of a' 
        # interpolation: next-period consumption
        # c_polate = interpolate.interp1d(agrid,cpolicy0[iperm,itheta1,:,iage+1], kind=kind_tol)
        # c1 = c_polate(a1)
        # print(c1)
        c1 = lininter(agrid,cpolicy0[iperm,itheta1,:,iage+1],a1)
        # check if interpolate too slow
        rf = rf - factor*probth[itheta,itheta1]*uc(c1,0.0)
    
    return rf



# rfoldss1(): same as rfoldss, however: modified for input into fsolve(.)
# first-order condition for retired agent: steady state 
# returns residual function (=0 in the optimum)
# input: k' (next-period asset)
# output: residual of the Euler equation
def rfoldss1(x,*args):
    a1 = x
    a0 = args[0]
    iage = args[1]
    iperm = args[2]
    itheta = args[3]
    pen0 = args[4]
    trbar0 = args[5]
    rbar0 = args[6]
    cpolicy0 = args[7]
    
    c0 = trbar0+pen0+(1+(1-tauk)*rbar0)*a0-ygrowth*a1
    c0=c0/(1+tauc)
    
    if c0<=0:
        return (1+c0**2)*1e5 # penalty function 
    
    rf=uc(c0,0.0)
    factor=sp1[iage]*beta1*(1+(1-tauk)*rbar0)
 
    # here, we assume that the idiosyncratic component z[iz] remains
    # stochastic during retirement --> innocuous assumption because
    # pension does not depend on z[iz]
    for itheta1 in range(ntheta):
        # linear interpolation between two grid points of a' 
        # interpolation: next-period consumption
        # if a1<=0:
        #    c1=cpolicy0[iperm,itheta1,0,iage+1]
        # else:
            # c_polate = interpolate.interp1d(agrid,cpolicy0[iperm,itheta1,:,iage+1], kind=kind_tol)
            # c1 = c_polate(a1)
        c1 = lininter(agrid,cpolicy0[iperm,itheta1,:,iage+1],a1)
        # check if interpolate too slow
        rf = rf - factor*probth[itheta,itheta1]*uc(c1,0.0)
    
    return rf

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
nq = 50        # number of outer iterations
#nq=1
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

bigkold=0


q=-1
crit = 1+tol
# loop over aggregate variables K, L, trbar
while q<nq-1 and crit>tol:
    q=q+1
    ass,css,lss = getpolicyss(pen,rbar,trbar,wbar)

    # computation of aggregate employment and savings 
    ga = np.zeros((nperm,ntheta,na,nage))    # distribution of assets (capital) 
    bigl=0  # aggregate employment 
    bigc=0  # aggregate consumption
    bigk=0  # aggregate assets k 
    lbar=0  # average labor supply 
    bequests=0     # accidental bequests
    totalmass=0
    
    
    # initialization: assets at age 1 are 0 
    for iperm in range(nperm):
        for itheta in range(ntheta):
            # assumptions: all productivities have the same measure 
            ga[iperm,itheta,0,0] = mass[0]/(nperm*ntheta)
            
    for iage in range(nage-1):
        print(iage)
        for ia in range(na):
            for iperm in range(nperm):
                for itheta in range(ntheta):
                    mass0 = ga[iperm,itheta,ia,iage]*sp1[iage]/(1+popgrowth)
                    a1 = ass[iperm,itheta,ia,iage]
                    if a1<=amin:
                        for itheta1 in range(ntheta):
                            ga[iperm,itheta1,0,iage+1] = ga[iperm,itheta1,0,iage+1]+probth[itheta,itheta1]*mass0
                        
                    elif a1>=amax:
                        for itheta1 in range(ntheta):
                            ga[iperm,itheta1,na-1,iage+1] = ga[iperm,itheta1,na-1,iage+1]+probth[itheta,itheta1]*mass0
                        
                    else:  # a1 lies between agrid[0] and agrid[na-1]
                        ia1=sum(agrid<a1) # a1 lies between agrid[ia-2] and agrid[ia-1]
                        is1=(a1-agrid[ia1-1])/(agrid[ia1]-agrid[ia1-1])
                        for itheta1 in range(ntheta):
                            ga[iperm,itheta1,ia1,iage+1] = ga[iperm,itheta1,ia1,iage+1]+is1*probth[itheta,itheta1]*mass0
                            ga[iperm,itheta1,ia1-1,iage+1] = ga[iperm,itheta1,ia1-1,iage+1]+(1-is1)*probth[itheta,itheta1]*mass0
                        
      
    grossysum=0
    ysum=0    
    gag= np.zeros(nage) # distribution assets between cohorts
    glg = np.zeros(nage) # distribution labor
    gyg = np.zeros(nage)    # distribution income 
    fa = np.zeros(na) # distribution ssets
    gytot = np.zeros(ny)    # distribution total income
    gynet = np.zeros(ny)    # and net income
    totalpen = sum(mass[nw:nage])*pen # total spending on pensions
    totalwage=0              

    for iage in range(nage):
        for ia in range(na):
            for iperm in range(nperm):
                for itheta in range(ntheta):
                    mass0 = ga[iperm,itheta,ia,iage]
                    totalmass = totalmass+mass0
                    a1 = ass[iperm,itheta,ia,iage]
                    c1 = css[iperm,itheta,ia,iage]
                    labor = lss[iperm,itheta,ia,iage]
                    if iage>nw-1:
                        if pentype==0:
                            pene=pen
                        else:
                            pene = pen*perm[iperm]
                        w0 = pene
                        wnet = w0 # net income, no taxes on pensions
                     
                    else:
                        w0= wbar*perm[iperm]*theta[itheta]*ef[iage]*labor
                        wnet=(1-taulbar)*w0
                    
                    ytot0 = rbar*agrid[ia]+w0+trbar
                    ynet0 = (1-tauk)*rbar*agrid[ia]+wnet
                
                    # distribution gross and net income
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
                
                    
                    
                    fa[ia] = fa[ia]+mass0
                    bigk = bigk+mass0*agrid[ia]
                    bigc = bigc + mass0*c1
                    bequests = bequests + mass0*(1-sp1[iage])*a1*ygrowth
                
                    if iage<=nw-1:
                        lbar = lbar+mass0*labor
                        bigl=bigl+mass0*labor*perm[iperm]*theta[itheta]*ef[iage]
                        totalwage=totalwage+mass0*labor*perm[iperm]*theta[itheta]*ef[iage]*wbar
                    
                
                    gag[iage] = gag[iage]+mass0*agrid[ia]
                    gyg[iage] = gyg[iage]+ytot0*mass0
                    
                    if iage<=nw-1:
                        glg[iage] = glg[iage]+mass0*labor
       


    print("bigk: " +str(bigk))
    crit=abs((bigkold-bigk)/bigk)
    bigkold=bigk
    
    
    # Gini coefficients of distributions
    giniwealth=ginid(agrid,fa,na)
    print("gini wealth: " +str(giniwealth))
    giniynet=ginid(ynetgrid,gynet,ny)
    print("gini net income: " +str(giniynet))
    giniygross=ginid(ytotgrid,gytot,ny)
    print("gini gross income: " +str(giniygross))
    
    # aggregation of L, tr, ... and update
    lbar=lbar/sum(mass[0:nw])
    lbar0=psi0*lbar0+(1-psi0)*lbar
    print("average labor supply: " +str(lbar))
    print("L: " +str(bigl))

    ybar = production(bigk,bigl,1.0)
    rbarnew = interest(bigk, bigl, 1.0)
    print("rbarnew: " +str(rbarnew))
    
    rbar=psi0*rbar+(1-psi0)*rbarnew
    wbarnew= wage(bigk,bigl,1.0)    
    print("wbarnew: " +str(wbarnew))
    wbar=psi0*wbar+(1-psi0)*wbarnew
    pen=replacement_ratio*(1-taulbar)*lbar0*wbar
    
    # update of transfers tr from the fiscal + social security budget 
    gbar = gybar*ybar
    trbarnew = bequests+taulbar*wbar*bigl+tauk*rbar*bigk+tauc*bigc-gbar-totalpen
    trbar=psi0*trbar+(1-psi0)*trbarnew
    print("trbar: " +str(trbar))    
    
    
# save results
np.save('ass',ass)
np.save('css',css)
np.save('lss',lss)
np.save('ga',ga)
np.save('lbar0',lbar0)
np.save('bequests',bequests)
np.save('gbar',gbar)



print("runtime: --- %s seconds ---" % (time.time() - start_time))    
sec = (time.time() - start_time)
ty_res = time.gmtime(sec)
res = time.strftime("%H : %M : %S", ty_res)
print(res)