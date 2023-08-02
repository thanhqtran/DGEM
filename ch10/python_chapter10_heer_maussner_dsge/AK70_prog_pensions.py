# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 14:36:08 2021

@author: heerburk
"""



# Part 1: import libraries
import pandas as pd
import numpy as np
#import scipy.optimize 
import quantecon as qe
from scipy.stats import norm
from scipy import interpolate
import time
import math
import matplotlib.pyplot as plt

# part 2: definition of the functions
#
# wage function
def wage(k,l):
    return (1-alpha) * k**alpha * l**(-alpha)

# interest rate function
def interest(k,l):
    return alpha * k**(alpha - 1) * l**(1-alpha)

# production function
def production(k,l):
    return k**(alpha) * l**(1-alpha)

# pension schedule
def pension(x,aincome,pen0):    # see Huggett and Parra (2010, JPE)   
    bendpoint1 = 0.20*aincome
    bendpoint2 = 1.24*aincome
    bmin = pen0
    
    if x<bendpoint1:
        y = bmin+0.9*x
    elif x>=bendpoint1 and x<=bendpoint2:
        y = bmin + 0.9*bendpoint1 + 0.32*(x-bendpoint1)
    else:
        y = bmin+ 0.9*bendpoint1 + 0.32*(bendpoint2-bendpoint1) + 0.15*(x-bendpoint2)
    
    return y



# utility function 
def u(c,l): 
#    c[c<=0.0] = psi1
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


# marginal utility from consumption next period given a',iperm,iy1,iage1
# in next period
#  uc1(k1,iperm,iy,iy1,ice,e1,iage+1,labor)
def uc1(a1,iperm0,iy0,iy01,ice0,e10,iage01,labor0):
    
    ia1 = sum(a<=a1)-1
    ie1 = sum(e<=e10)
    ie1 = max(1,ie1)
    ie1 = min(nce,ie1)
    
    if iage01<=nw-1:
    
        if e10<=0: 
            c1 = cwopt[ia1,0,iy01,iperm0,iage01]
            labor = lopt[ia1,0,iy01,iperm0,iage01]
        elif e10>=emax:
            c1 = cwopt[ia1,nce-1,iy01,iperm0,iage01]
            labor =  lopt[ia1,nce-1,iy01,iperm0,iage01]
        else:    
            c_polate = interpolate.interp1d(e,cwopt[ia1,:,iy01,iperm0,iage01], kind=kind_tol)
            c1 = c_polate(e10)
            l_polate = interpolate.interp1d(e,lopt[ia1,:,iy01,iperm0,iage01], kind=kind_tol)
            labor = l_polate(e10)
    else:
        
        labor = 0
        if e10<=0: 
            c1 = cropt[ia1,0,iage01-nw]
        elif e10>=emax:
            c1 = cropt[ia1,nce-1,iage01-nw]
        else:
            c_polate = interpolate.interp1d(e,cropt[ia1,:,iage01-nw], kind=kind_tol)
            c1 = c_polate(e10)
            
    y = uc(c1,labor)
    return y



def optimal_labor(a0,a1):
    w0=(1-taun-taup)*w*ef[iage]*perm[iperm]*ye1[iy]
    labor = gamma - ((1+(1-tauk)*(r-delta))*a0+trbar-ygrowth*a1)*(1-gamma)/w0    
    return labor

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



start_time = time.time()

# abbreviations
exp = np.e
log = math.log

# Step 1.1: Import data
data = pd.read_excel (r'survival_probs.xlsx') 
df = pd.DataFrame(data, columns= ['sp1','ef'])
print(df)
arr = np.array(df)
print(arr)
sp1 = arr[:,0]
ef = arr[0:45,1]
print(ef)


# Step 1.2: Parameterization of the model
# demographics
nage=70                # maximum age                  
nw=45                  # number of working years        
Rage=46                # first period of retirement 
nr=nage-Rage+1         # number of retirement years
popgrowth0 = 0.0075400000 # population growth rate


# preferences
beta1=1.011             # discount factor 
gamma=0.33		        # weight of consumption in utility
lbar=0.25		        # steady-state labor supply
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
replacement_ratio=0.352	# gross replacement ratio US
bybar=0.63					# debt-output ratio
gybar=0.18					# government consumption-output ratio
penminy=0.1242          # minimum pension as a percentage of average earnings

kind_tol = 'linear'     # interpolation of policy functions:
                        # 'linear' or 'Cubic'
phi=0.80       # update aggregate variables in outer iteration over K, L, tr, taup, taun
tol=0.001     # percentage deviation of final solution K and L
tol_final=tol
psi1 = 0.000001   # minimum consumption
nq = 50        # number of outer iterations
crit_counter=0  # is set to 1 if na is increased to na1=400
#nq=1


# productivity of workers

nperm=2;        # number of permanent productivities
perm = np.zeros(2)
perm[0]=0.57
perm[1]=1.43

lamb0=.96           # autoregressive parameter 
sigmay1=0.38        # variance for 20-year old, log earnings */
sigmae=0.045        # earnings disturbance term variance */
ny=5                # number of productivity types
m=1                 # width of the productivity grid
                    # -> calibrated to replicate Gini wages=0.37


# compute productivity grid
sy = np.sqrt(sigmae/(1-lamb0**2))    # earnings variance 
ye =  np.linspace(-m*sy, m*sy, ny)
print(ye)



# transition matrix using Tauchen's approximation
# return is a class 'Markov Chain'
mc = qe.markov.approximation.tauchen(lamb0, sigmae, 0.0,  m, ny)
# transition matrix is stored in object 'P
py = mc.P

# mass of the workers
muy = np.zeros(ny)
w = ye[1]-ye[0]

# first year mass distribution
muy[0] = norm.cdf( (ye[0]+w/2)/np.sqrt(sigmay1) ) 
muy[ny-1] = 1-norm.cdf( (ye[ny-1]-w/2)/np.sqrt(sigmay1))


for i in range(1,ny-1):
    muy[i] = ( norm.cdf( (ye[i]+w/2)/np.sqrt(sigmay1) ) - 
        norm.cdf( (ye[i]-w/2)/np.sqrt(sigmay1) ) )
    


#    muy[i,:] = muy[i-1,:] @ py * sp1[i] / (1+popgrowth0)

# transform ye so that mean exp(ye) = 1
# (mean efficiency equal to one)

ye1=np.exp(ye)
meane=np.dot(muy,ye1)
#ye = ye-log(meane)
print(ye)
ye1=np.exp(ye)

# asset grids
kmin=0              # inidividual wealth
kmax=20             # upper limit of capital grid 
na=200              # number of grid points over assets a in [kmin,kmax]
# na=400
na1=400             # number of grid points during final iterations over aggregates
# na1=600
a =  np.linspace(kmin, kmax, na)   # asset grid policy function


nce=10           # individual earnings grid point number
nce1=10
emin=0
emax = 3.0
e =  np.linspace(emin, emax, nce)   # grid over cumulated earnings

labormin=0
labormax=0.6        # maximum labor supply
nlabor = 30
nlabor1 = 60
lgrid =  np.linspace(labormin, labormax, nlabor)   # grid over labor

ymin=0
ymax=4.0
nin=100
incomegrid =  np.linspace(ymin, ymax, nin)   # grid over income

# measure of living persons
mass = np.ones(nage)

for i in range(nage-1):
    mass[i+1]=mass[i]*sp1[i]/(1+popgrowth0)

mass = mass / mass.sum()


# save measures
np.save('mass',mass)
np.save('muy',muy)  

# initial guesses
#
rbar=0.03   # interest rate
nbar=0.30    # aggregate efficienct labor L
nold=100    # initialization 
mean_labor=0.3   # average working hours
kbar=(alpha/(rbar+delta))**(1/(1-alpha))*nbar # capital stock 
kold=100
omega = kbar*1.2    # aggregate wealth
trbar=0.01          # transfers, initial guess 
w=wage(kbar,nbar) 
r=interest(kbar,nbar)
rb = (1-tauk)*(r-delta)
pen=replacement_ratio*(1-taulbar)*w*mean_labor*sum(mass)/sum(mass[0:nw])
taup=pen*sum(mass[nw:nage])/sum(mass)/(w*nbar)  # balanced budet social security
taun = taulbar-taup
bequests=0

#
# computation of the Gini coefficients of hourly wages    
# wage inequality
#
gwages = np.zeros((nw,nperm,ny,2))  # distribution of wages

# initialization of wage distribution at age 1
for iperm in range(nperm):
    for iy in range(ny):
        gwages[0,iperm,iy,0] = ye1[iy]*perm[iperm]*ef[0] # hourly wage at age 1
        gwages[0,iperm,iy,1] = 1/2*muy[iy]*mass[0]/sum(mass[0:nw]) # measure of households

for i in range(1,nw,1):
    print(i)
    for iperm in range(nperm):
        for iy in range(ny):
            gwages[i,iperm,iy,0] = ye1[iy]*perm[iperm]*ef[i]
            for iy1 in range(ny):
                gwages[i,iperm,iy1,1] = gwages[i,iperm,iy1,1]+sp1[i-1]/(1+popgrowth0)*py[iy,iy1]*gwages[i-1,iperm,iy,1]

# wages at age nw
for iperm in range(nperm):
    for iy in range(ny):
        gwages[nw-1,iperm,iy,0] = ye1[iy]*perm[iperm]*ef[nw-1] # hourly wage at age 1                

i0=-1
fwage=np.zeros((nw*nperm*ny,2))
for i in range(nw):
    for iperm in range(nperm):
        for iy in range(ny):
            i0=i0+1
            fwage[i0,0]=gwages[i,iperm,iy,0]
            fwage[i0,1]=gwages[i,iperm,iy,1]
            
#
# sorting according to first column 
# --> preparing for function ginid()
#
y0 = fwage[fwage[:,0].argsort()]
x = y0[:,0]
g = y0[:,1]
print("Gini_w= " + str(ginid(x,g,nw*nperm*ny)))


#
# outer iteraton over aggregates K, L, tr, taup: q = 0,..,nq
#
#


# saves aggregates during iteration: K, L, Omega, mean labor, taun, taup, tr, bequests
aggregateq = np.zeros((nq,8))

q=-1
crit = 1+tol
while q<nq-1 and crit>tol:
    q=q+1
    print("q: " + str(q))
    crit=max(abs((kbar-kold)/kbar),abs((nbar-nold)/nbar)) # percentage deviation of K and L below tol
    aggregateq[q,0] = kbar
    aggregateq[q,1] = nbar
    aggregateq[q,2] = omega
    aggregateq[q,3] = mean_labor
    aggregateq[q,4] = taun
    aggregateq[q,5] = taup
    aggregateq[q,6] = bequests
    aggregateq[q,7] = trbar
    
    w=wage(kbar,nbar)
    r=interest(kbar,nbar)   # marginal product of capital
    rb = (1-tauk)*(r-delta) # interest rate on bonds
    ybar=production(kbar,nbar)
    debt = bybar*ybar
    gbar = gybar*ybar
    meanincome=w*nbar/sum(mass[0:nw])  # average earnings of the workers
    penmin=penminy*meanincome             # minimum pension
    kold=kbar
    nold=nbar


    # 
    # computation of the policy function: value function iteration
    #
    
    
    policy_time = time.time()   # start timer for policy function computation
    
    vr=np.zeros((na,nce,nr))    # value function with lump-sum pensions: only depends on assets a
    aropt=np.zeros((na,nce,nr))  # optimal asset 
    cropt=np.zeros((na,nce,nr)) # optimal consumption 

    # workers' value function 
    vw = np.zeros((na,nce,ny,nperm,nw))
    awopt = np.zeros((na,nce,ny,nperm,nw))
    lopt = np.zeros((na,nce,ny,nperm,nw))
    cwopt = np.zeros((na,nce,ny,nperm,nw))

    for ia in range(na):
        for ice in range(nce):
            pen = pension(e[ice],meanincome,penmin)
            c = a[ia]*(1+(1-tauk)*(r-delta))+pen+trbar
            c = c/(1+tauc)
            vr[ia,ice,nr-1] = u(c,0)
            cropt[ia,ice,nr-1] = c
            aropt[ia,ice,nr-1] = 0

    for iage in range(nr-1,0,-1):
        print(iage)
        zd = vr[:,:,iage]
        for ia in range(na):
            for ice in range(nce):
                c = (1+rb)*a[ia]+pension(e[ice],meanincome,penmin)+trbar-ygrowth*a
                c = c/(1+tauc)
                c[c<=0] = psi1      # set all non-positive entries equal to psi1 
                y = u(c,np.zeros(na))+ygrowth**(gamma*(1-eta))*beta1*sp1[nw+iage-1]*zd[:,ice]
                indexopt = np.where(y == np.max(y))
                k1 = a[indexopt]
                aropt[ia,ice,iage-1] = k1
                vr[ia,ice,iage-1] = y[indexopt]
                cropt[ia,ice,iage-1]=( (1+rb)*a[ia]+pension(e[ice],meanincome,penmin)+trbar-ygrowth*k1 )/(1+tauc)


    #  compuate decsion rules of workers

    for i in range (nw,0,-1):
        print("Compute the policy function of the worker of age= " + str(i))
        print("K: " + str(kbar))    
        print("q: " + str(q))
        print("crit: " +str(crit))
        print("na: " + str(na))
        for iperm in range(nperm):
            if i==nw:           # next-period value function is first year of retirement
                y=vr[:,:,0]       # value function of retired at retirement period i+1
                #y.shape = (nce,na)  
                #y = np.transpose(y) # check if y correct (na x nce)-matrix 
            #else:
             #   for iy1 in range(ny):
             #      print(iy1)
              #      # the na x nce variables of the productivity type prod are stored in y 
               #     y=vw[:,:,iy1,iperm,i]    
                #    y.shape = (nce,na)  
                  #  # y = np.transpose(y) 
                 #   if iy1==0:
                   #     zd1 = y
            #        elif iy1==1:
             #           zd2 = y
             #       elif iy1==2:
             #           zd3 = y
             #       elif iy1==3:
             #           zd4 = y
             #       elif iy1==4:
             #           zd5 = y
             #       else:
             #           print("not implemented yet")
             #            # index
                         
            for iy in range(ny):  # productivity at age i
                ymax = w*ef[i-1]*perm[iperm]*ye1[iy]*labormax     # maximum labor income
                for ia in range(na):     # asset leval at age i and productivity iy
                    k0=a[ia]
                    for ice in range(nce):   #  cumulated labor earnings at age i                
                        e0=e[ice] 
                        # next-period cumulated earnings  
                        e1 = e0*(i-1)/i + w*ef[i-1]*perm[iperm]*ye1[iy]*lgrid/i
                        labor1 = lgrid*np.ones((na,nlabor))
                        zz = a*np.ones((nlabor,na)) 
                        zz = np.transpose(zz)
                        c0 = (1+rb)*k0 +(1-taup-taun)*w*ef[i-1]*perm[iperm]*ye1[iy]*labor1
                        c0 = c0 + trbar - ygrowth*zz
                        c0 = c0/(1+tauc)
                        c0[c0<0] = psi1
                        if i==nw:
                            # contruction of the next period value function 
                            for ilabor in range(nlabor):
                                e10=e1[ilabor]
                                # interpolation of value function at e10
                                # values of e below e10
                                i1 = np.count_nonzero([e<e10])
                                if i1<1:
                                    i1 = 1
                                elif i1 > nce-1:
                                    i1 = nce-1
                            
                                z1 = (e10-e[i1-1])/(e[i1]-e[i1-1])
                            
                                # stores the linear interpolated next-period value function for e1[ilabor]
                                v0 = (1-z1)*y[:,i1-1] + z1*y[:,i1]
                            
                                # np.c_[a,b] column stack, np.r_[a,b] row stack
                                if ilabor==0:
                                    vtemp = v0
                                else:
                                    vtemp = np.c_[vtemp,v0]

                            bellman = u(c0,labor1)
                            bellman = bellman+ygrowth**(gamma*(1-eta))*sp1[i-1]*beta1*vtemp
                            # index where bellman eq attains maximum
                            zmax = np.where(bellman == np.max(bellman))
                            iamax = zmax[0]
                            ilabormax = zmax[1]
                            k1 = a[iamax]
                            e1 = e1[ilabormax]
                            labor = lgrid[ilabormax]
                            c0 = (1+rb)*k0+(1-taup-taun)*w*ef[i-1]*perm[iperm]*ye1[iy]*labor+trbar-ygrowth*k1                    
                            c0 = c0/(1+tauc)
                            vw[ia,ice,iy,iperm,i-1] = bellman[iamax,ilabormax]
                            awopt[ia,ice,iy,iperm,i-1] = k1
                            cwopt[ia,ice,iy,iperm,i-1] = c0
                            lopt[ia,ice,iy,iperm,i-1] = labor 
                        
                        else:   # i<nw
                            vtemp = np.zeros((na,nlabor))
                            # contruction of the next period value function ated next-period value function for e1[ilabor]
                            for iy1 in range(ny):
                                # the na x nce variables of the productivity 
                                # type prod are stored in y     
                                y = vw[:,:,iy1,iperm,i]  
                                #y.shape = (nce,na)  
                                #y = np.transpose(y)
                                for ilabor in range(nlabor):
                                    e10=e1[ilabor]
                                    # interpolation of value function at e10
                                    # values of e below e10
                                    i1 = np.count_nonzero([e<e10])
                                    if i1<1:
                                        i1 = 1
                                    elif i1 > nce-1:
                                        i1 = nce-1
                            
                                    z1 = (e10-e[i1-1])/(e[i1]-e[i1-1])
                            
                                    # stores the linear interpolated next-period value function for e1[ilabor]
                                    v0 = (1-z1)*y[:,i1-1] + z1*y[:,i1]
                            
                                    # np.c_[a,b] column stack, np.r_[a,b] row stack
                                    if ilabor==0:
                                        vtemp0 = v0
                                    else:
                                        vtemp0 = np.c_[vtemp0,v0]
                            
                                vtemp = vtemp + py[iy,iy1]*vtemp0;
                            
                       
                            bellman = u(c0,labor1)
                            bellman = bellman+ygrowth**(gamma*(1-eta))*sp1[i-1]*beta1*vtemp
                            # index where bellman eq attains maximum
                            zmax = np.where(bellman == np.max(bellman))
                            iamax = zmax[0]
                            ilabormax = zmax[1]
                            k1 = a[iamax]
                            e1 = e1[ilabormax]
                            labor = lgrid[ilabormax]
                            c0 = (1+rb)*k0+(1-taup-taun)*w*ef[i-1]*perm[iperm]*ye1[iy]*labor+trbar-ygrowth*k1
                            c0 = c0/(1+tauc)
                            vw[ia,ice,iy,iperm,i-1] = bellman[iamax,ilabormax]
                            awopt[ia,ice,iy,iperm,i-1] = k1
                            cwopt[ia,ice,iy,iperm,i-1] = c0
                            lopt[ia,ice,iy,iperm,i-1] = labor 
                        

              
    # save results
    np.save('vw',vw)
    np.save('vr',vr)
    np.save('cwopt',cwopt)
    np.save('cropt',cropt)
    np.save('awopt',awopt)
    np.save('aropt',aropt)
    np.save('lopt',lopt)     

    print("runtime value function: --- %s seconds ---" % (time.time() - start_time))    
    sec = (time.time() - start_time)
    ty_res = time.gmtime(sec)
    res = time.strftime("%H : %M : %S", ty_res)
    print(res)     


    print("runtime value function: --- %s seconds ---" % (time.time() - policy_time))    
    sec = (time.time() - policy_time)
    ty_res = time.gmtime(sec)
    res = time.strftime("%H : %M : %S", ty_res)
    print(res)     


    # --------------------------------------------------------------
    #
    #
    # computation of the distribution of capital 
    #
    #    
    # ---------------------------------------------------------------


    distribution_time = time.time()     

    ga = np.zeros((na,nce,ny,nperm,nage)) # distribution of wealth 
                # assumption: at age nw+1, the household's productivity remains constant
    agen = np.zeros(nage)             # distribution of wealth over age
    gwealth = np.zeros(na)       # distribution of wealth 
    fy = np.zeros(nin)          # distribution of income (all households)
    gk = np.zeros((na,nw))
    fa = np.zeros(na)       # distribution of assets
    fe = np.zeros(nce)      # distribution of accumulated earnings
    
    # age-profiles
    bigl=0				    # aggregate labor 
    biga=0					# aggregate assets
    bigc=0                  # aggregate consumption
    bigcontrib=0			# aggregate contributions to pension system
    bigbequests=0			# aggregate accidental bequests
    bigpensions=0			# total pensions
    mean_labor=0
    
    
    #
    # initialization of distribution function at age 1
    #
    # all households start with zero wealth
    #
    measure_check=0        # consistency check: if all measures add up to one
    
    for iperm in range(nperm):
        for iy in range(ny):
            ga[0,0,iy,iperm,0] = muy[iy]*mass[0]/nperm   # all permanent types have same measure
        
        

    #
    # distribution at age 2,...,nage-1
    #
    for iage in range(nage-1):
        print("distribution")
        print("iage: " +str(iage))       
        for iperm in range(nperm):
            print("iperm: " +str(iperm))
            for iy in range(ny):    #  present-period idiosyncratic productivity
                for ia in range(na):
                    asset0 = a[ia]
                    for ice in range(nce):
                        e0 = e[ice]
                        measure = ga[ia,ice,iy,iperm,iage]
                        measure_check = measure_check + measure
                        agen[iage] = agen[iage] + measure*asset0
                        fa[ia] = fa[ia] + measure
                        if iage==nw:  # distribution of e at beginning of retirement
                            fe[ice] = fe[ice] + measure/mass[iage] 
                            
                        if iage<=nw-1:  # worker?
                            labor = lopt[ia,ice,iy,iperm,iage]
                            mean_labor = mean_labor + labor*measure/sum(mass[0:nw])
                        
                        if iage<=nw-1:
                            c = cwopt[ia,ice,iy,iperm,iage]
                            a1 = awopt[ia,ice,iy,iperm,iage]
                        else:
                            c = cropt[ia,ice,iage-nw]
                            a1 = aropt[ia,ice,iage-nw]
                        
                        bigc = bigc + c*measure
                        
                        if iage<=nw-1:
                            e1 = iage/(iage+1)*e0+w*perm[iperm]*ef[iage]*ye1[iy]*labor/(iage+1)
                            bigl = bigl + perm[iperm]*ef[iage]*ye1[iy]*labor*measure
                            bigcontrib = bigcontrib + taup*w*perm[iperm]*ef[iage]*ye1[iy]*labor*measure
                            # gross income gy
                            gy = w*perm[iperm]*ef[iage]*ye1[iy]*labor  + (r-delta)*asset0
                        else:
                            e1 = e0
                            bigpensions = bigpensions + pension(e0,meanincome,penmin)*measure
                            gy = pension(e0,meanincome,penmin) + (r-delta)*asset0
        
                    
                        if gy<=0:
                            fy[0] = fy[0] +  measure
                        elif gy>=incomegrid[nin-1]:
                            fy[nin-1] = fy[nin-1] + measure
                        else:   # linear interpolation between grid points
                            j0=sum(incomegrid<gy)
                            j0=min(j0,nin-1)
                            n0=(gy-incomegrid[j0-1])/(incomegrid[j0]-incomegrid[j0-1])
                            fy[j0-1] = fy[j0-1]+ (1-n0)*measure
                            fy[j0] = fy[j0]+ n0*measure
                
                        bigbequests = bigbequests + (1+rb)*a1*(1-sp1[iage])/(1+popgrowth0) * measure
                        biga = biga+asset0*measure
                        
                        
                        ia1 = sum(a<=a1)
                        ia1 = max(1,ia1)
                        ia1 = min(na,ia1)
                        
                        				
                        if e1<=0: # linear interpolation between the two adjacent grid points on a1 on agrid
                            for iy1 in range(ny):
                                ga[ia1-1,0,iy1,iperm,iage+1] = ga[ia1-1,0,iy1,iperm,iage+1] + py[iy,iy1]*sp1[iage]/(1+popgrowth0)*measure  
                                                             
                        elif e1>=emax:
                            for iy1 in range(ny):
                                ga[ia1-1,nce-1,iy1,iperm,iage+1] = ga[ia1-1,nce-1,iy1,iperm,iage+1] + py[iy,iy1]*sp1[iage]/(1+popgrowth0)*measure  
                                
                        else: # linear interpolation                  
                            ice1 = sum(e<e1)
                            lambda1 = (e1-e[ice1-1]) / (e[ice1]-e[ice1-1] )
                            for iy1 in range(ny):
                                ga[ia1-1,ice1-1,iy1,iperm,iage+1] = ga[ia1-1,ice1-1,iy1,iperm,iage+1] + (1-lambda1)*py[iy,iy1]*sp1[iage]/(1+popgrowth0)*measure  
                                ga[ia1-1,ice1,iy1,iperm,iage+1] = ga[ia1-1,ice1,iy1,iperm,iage+1] + lambda1*py[iy,iy1]*sp1[iage]/(1+popgrowth0)*measure  
                                

    #
    # distribution at age nage: last period of life
    #
    iage = nage-1
    print("distribution")
    print("iage: " +str(iage))       
    for iperm in range(nperm):
        print("iperm: " +str(iperm))
        for iy in range(ny):    #  present-period idiosyncratic productivity
            for ia in range(na):
                asset0 = a[ia]
                for ice in range(nce):
                     e0 = e[ice]
                     measure = ga[ia,ice,iy,iperm,iage]
                     measure_check = measure_check + measure
                     agen[iage] = agen[iage] + measure*asset0
                     fa[ia] = fa[ia] + measure
                            
                     c = cropt[ia,ice,iage-nw]
                     a1 = aropt[ia,ice,iage-nw]
                        
                     bigc = bigc + c*measure
                        
                     e1 = e0
                     bigpensions = bigpensions + pension(e0,meanincome,penmin)*measure
                     gy = pension(e0,meanincome,penmin) + (r-delta)*asset0
                     if gy<=0:
                         fy[0] = fy[0] +  measure
                     elif gy>=incomegrid[nin-1]:
                         fy[nin-1] = fy[nin-1] + measure
                     else:   # linear interpolation between grid points
                        j0=sum(incomegrid<gy)
                        j0=min(j0,nin-1)
                        n0=(gy-incomegrid[j0-1])/(incomegrid[j0]-incomegrid[j0-1])
                        fy[j0-1] = fy[j0-1]+ (1-n0)*measure
                        fy[j0] = fy[j0]+ n0*measure
                
                        biga = biga+asset0*measure
                        
   
    #
    # Compute new initial values
    #    


    # update of trbar from the fiscal budget constraint
    taxes=taun*w*nbar+tauk*(r-delta)*kbar+tauc*bigc    
    transfernew = taxes + bigbequests + debt*((1+popgrowth0)*ygrowth-(1+rb)) - gbar	
    trbar = phi*trbar + (1-phi)*transfernew
        
    # total savings
    omeganew = sum(agen) /sum(mass)
    ybar = production(kbar,nbar)
    debt = bybar*ybar
    gbar = gybar*ybar
    omega = phi*omega + (1-phi)*omeganew
    knew = omeganew - debt
    kbar=phi*kbar+(1-phi)*knew
    nbar = phi*nbar + (1-phi)*bigl

    # update of taup 
    taupnew = bigpensions/(w*nbar)
    taup = phi*taup + (1-phi)*taupnew
    taun = taulbar-taup    # calibration so that tau^n +tau^p = 28% as in the US

    # average pension / average wage income
    replacement_rate_new = bigpensions / sum(mass[nw:nage])*sum(mass[0:nw])/(w*nbar)
    crit=abs((kbar-knew)/kbar)
         
    print("runtime distribution function: --- %s seconds ---" % (time.time() - distribution_time))    
    sec = (time.time() - distribution_time)
    ty_res = time.gmtime(sec)
    res = time.strftime("%H : %M : %S", ty_res)
    print(res)     
    
    totalmeasure = np.sum(np.sum(ga))
    print("totalmeasure = " + str(totalmeasure))
    print("measurecheck = " + str(measure_check))
    
    
    gini_wealth = ginid(a,fa,na)
    print("gini wealth: " +str(gini_wealth))
    gini_income = ginid(incomegrid,fy,nin)
    print("gini income: " +str(gini_income))
        

    if q==0 or crit<tol: 
        # Residual: Euler equation         
        Euler_res = np.zeros((na,nce,ny,nperm,nw))
        Euler_res_old = np.zeros((na,nce,nr-1))
  
        for iage in range(nage-1):
            for iperm in range(nperm):
                for iy in range(ny):
                    for ia in range(na):
                        for ice in range(nce):
                            e0 = e[ice]
                            
                            if iage<=nw-1:
                                k1=awopt[ia,ice,iy,iperm,iage]
                                labor = lopt[ia,ice,iy,iperm,iage]
                                c=cwopt[ia,ice,iy,iperm,iage]                                     
                                e1=(iage)/(iage+1)*e0+w*perm[iperm]*ef[iage]*ye1[iy]*labor/(iage+1)  
                            else:
                                k1=aropt[ia,ice,iage-nw]
                                labor = 0
                                c=cropt[ia,ice,iage-nw]
                                e1=e0
                                
                            # computation of the Euler residual
                            #
                            # (1+g_A)^eta u_{c,t} = beta phi^i E_t{ u_{c,t+1} [1+(1-tauk) (r-delta)] }
                            #
                    
                            x=0
                            if iage<=nw-1:
                                for iy1 in range(ny):
                                    x = x + beta1 * sp1[iage]*py[iy,iy1]* (1+rb)* uc1(k1,iperm,iy,iy1,ice,e1,iage+1,labor)
                                    
                            elif iage>nw-1 and iy==0:
                                x = beta1 * sp1[iage]*(1+rb)* uc1(k1,iperm,iy,iy1,ice,e1,iage+1,labor)
                                
                            
                            if iage<=nw-1: 
                                Euler_res[ia,ice,iy,iperm,iage] = 1-x / ( ygrowth**(1-gamma*(1-eta)) * uc(c,labor))
                            elif iage>nw-1 and iy==0:
                                Euler_res_old[ia,ice,iage-nw] = 1-x / ( ygrowth**(1-gamma*(1-eta)) * uc(c,labor))
                        
                      
            
        error_Euler = np.mean(np.mean(abs(Euler_res)))
        error_Euler_old = np.mean(np.mean(abs(Euler_res_old)))
        print("Euler error, workers: " + str(error_Euler))
        print("Euler error, retirees: " + str(error_Euler_old))
   
    if abs(crit)<tol and crit_counter==0:
        crit = abs(crit)+tol
        crit_counter=1
        tol = tol_final
        na=na1        # grid numbers close to convergence of aggregate capital stock
        a =  np.linspace(kmin, kmax, na)   # asset grid policy function
        nlabor=nlabor1  # grid numbers close to convergence of aggregate capital stock
        lgrid =  np.linspace(labormin, labormax, nlabor)   # grid over labor

# plotting the wealth averages at different ages 
plt.xlabel('age')
plt.ylabel('mean wealth')
plt.plot(range(21,91),agen/mass)
plt.show()

# plotting some arbitrary policy functions
iage=20
iy=3
iperm=1
ice=0
ia=0
plt.xlabel('individual wealth')
plt.ylabel('labor supply')
plt.plot(a,lopt[:,ice,iy,iperm,iage])
plt.show()

plt.xlabel('individual wealth')
plt.ylabel('consumption')
plt.plot(a,cwopt[:,ice,iy,iperm,iage])
plt.show()


plt.xlabel('individual wealth')
plt.ylabel('savings')
plt.plot(a,awopt[:,ice,iy,iperm,iage]-a)
plt.show()


plt.xlabel('accumulated earnings')
plt.ylabel('labor')
plt.plot(e,lopt[ia,:,iy,iperm,iage])
plt.show()

# plotting the Lorenz curve WEALTH
fig, ax = plt.subplots()
label1 = 'equal distribution'
label2 = 'wealth'
ag1 = fa*a
totalwealth = np.dot(fa,a)
ag1 = ag1/totalwealth
ag1 = np.cumsum(ag1)
ag1 = np.r_[0,ag1]
# add the point (0,0) to the figure
ag0 = np.cumsum(fa)
ag0 = np.r_[0,ag0]
ax.plot(ag0,ag0, linewidth=2, label=label1)
ax.plot(ag0,ag1, linewidth=2, label=label2)
ax.legend()
plt.show()



