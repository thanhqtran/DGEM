# AK70_prog_pensions.jl
#
# computes the steady in the OLG model (Auerbach-Kotlikoff)
# with stochastic idiosyncratic income and earnings-dependent pensions
#
# Method: DISCRETE VALUE FUNCTION ITERATION
#
# in Heer/Maussner, "Dynamic General Equilibrium Modeling", 2009,
# third edition, Section 10.1.4
#
# Date: August 4, 2021
#
# I would like to thank Chris Scharrer for
# helpful assistance. All errors are mine
#
# For questions send an email to: Burkhard.Heer@wiwi.uni-augsburg.de
#


#*************************
# load packages and procs:
#*************************
# IN CASE YOU HAVE NOT ADDED the JLD-package
# or XSLX Package
#Pkg.add("JLD")
# first type "using Pkg"
# next type "Pkg.add("XLSX")" on Julia command
# Pkg.add("NLopt")
using JLD
using Printf
using Plots
using QuantEcon
using Optim
using Roots
using NLsolve
using ForwardDiff
using Revise
using Distributions
using Random
using QuantEcon
using Random
using Interpolations
import Base.Threads.@spawn
using Pkg
import XLSX
#using Dierckx
include("AK70_prog_pensions_procs.jl");
#import Distributions: pdf, Normal, quantile

#**************
# output title:
#**************
println("")
println("")
println("")
println("")
println("*******************************************")
println("** Discrete Value Function Iteration: *****")
println("*******************************************")
println("******* Heer/Maussner: DGE Modeling *******")
println("******* Section 10.1.4 (3rd edition) ******")
println("******* Auerbach-Kotlikoff Model **********")
println("******* STOCHASTIC INDIVIDUAL INCOME ******")
println("******* + EARNINGS-DEPENDENT PENSIONS *****")
println("*******************************************")
println("")


#************
# Parameters:
#************

# numerical parameters
const psi1 = 0.000001     # parameter of utility function for zero or negative consumption
const phi=0.8       # parameter for the update of initial values
const tol=0.001         # percentage deviation of final solution with na=200
const tol_final=0.001   #  .. with na =400
const nq = 50        # number of outer iterations
crit_counter = 0    # is set to =1 if accuracy increases


# asset and income grids
const na0=200    # inital number of grid points on asset grid for value + policy function
const na1=400   # final number .....
const kmin=0.0  # minimum asset level
const kmax=20.0 # maximum asset level
const emin = 0.0 # minimum accumulated earnings
const emax = 3.0 # maximum accumulated earnings
const nce = 10  # number of grid points on accumulated earnings
const e = range(emin,stop=emax,length=nce)
const labormin = 0.0    # minimum labor
const labormax = 0.60   # maximum labor
const nlabor0 = 30   # initial grid points for labor
const nlabor1 = 60  # final grid point number
const nin=100   # number of grid points on income grid -> gini coefficent
const ymin=0.0  # minimum of earnings/income
const ymax=4.0  # maximum of earnings/income
const incomegrid = range(ymin,stop=ymax,length=nin)
na = na0    # initial grid sizes
nlabor = nlabor0
#a = range(kmin,stop=kmax,length=na)
dist=(kmax-kmin)/(na-1)
a=collect(kmin:dist:kmax)
dist=(labormax-labormin)/(nlabor-1)
laborgrid=collect(labormin:dist:labormax)

# model:
const beta=1.011     # discount factor
const eta=2.00      # 1/IES (intertemporal elasticity of substitution)
const alp=0.35      # production elasticity of capital
const delta=0.083     #  rate of depreciation
dist=(kmax-kmin)/(na-1)
a=collect(kmin:dist:kmax)
dist=(labormax-labormin)/(nlabor-1)
laborgrid=collect(labormin:dist:labormax)
const ygrowth=1.02	# annual growth factor
const gam=0.33      # relative weight of consumption in utility
const nage=70       # maximum age = 70
const nw=45         # number of working years
const Rage=46       # first period of retirement
const nr=25         # number of retirement years = 25
const popgrowth0 = 0.0075400000 # population growth rate

# fiscal policy and social security
const taulbar=0.28	 # both taun+taup=0.28!!, see Mendoza, Razin (1994), p. 311
const tauk=0.36      # capital income tax rate
const tauc=0.05              # consumption tax rate
taup=0.124           # initial guess for social security contribution rate
const bybar=0.63	# debt-output ratio
const gybar=0.18	# government consumption-output ratio
const peny=0.1242   # minimum pension (relative to mean income)

# read survival probabilities and age-efficiency profile from EXCEL
# file "survival_probs.xlsx"
sp1 = zeros(eltype(gam),nage+1)
ef1 = zeros(eltype(gam),nw)
(sp1,ef1) = read_data()
sp1 = Float64.(sp1)
ef = Float64.(ef1)
age = 1:1:nage

# permanent and stochastic productivity

nperm=2        # number of permanent productivities
perm = zeros(eltype(gam),2)
perm[1]=0.57
perm[2]=1.43


ef_rho=0.96         # autoregressive parameter
sigmay1=0.38        # variance for 20-year old, log earnings */
sigmae=0.045        # earnings disturbance term variance */
ny = 5              # number of productivity types
ef_mean = 0.0       # mean of (log) productivity
ef_m=1              # width of the productivity grid
                    # -> calibrated to replicate Gini wages=0.37

#
# function "tauchen" is taken from QuantEcon
#
ef_data=tauchen(ny,ef_rho,sqrt(sigmae),ef_mean,ef_m)
py = ef_data.p
ye = ef_data.state_values
ye1 = exp.(collect(ef_data.state_values))
ye1_org=collect(ef_data.state_values)



# ---------------------------------
#
#  initialization
#
# guess:
# labor - 0.3 (working hours times labor force)
# capital - so that r=3%
# wealth - approximately 80% in capital, 20% in government bonds
# trbar - small number
#
# --------------------------------- @

rbar=0.03   # initial interest rate
nbar=0.3    # aggregate efficienct labor L
mean_labor=0.3   # average working hours
kbar=(alp/(rbar+delta))^(1/(1-alp))*nbar # foc of firms with respect to capital
omega = kbar*1.2  # approx. 80% of savings omega are spent on capital (20% on bonds)
trbar=0.01  # transfers, initial guess
w=wage(kbar,nbar)
r=interest(kbar,nbar)
rb = (1-tauk)*(r-delta)


# cohort measures
mass=ones(nage)
for i in 2:nage
    mass[i] = mass[i-1]*sp1[i-1]/(1+popgrowth0)
end
mass = mass/sum(mass)

# initial guess for dependency ratio 1/3, therefore share of workers in
# population equal to 2/3
replacement_ratio=0.352
pen=replacement_ratio*(1-taulbar)*w*mean_labor*sum(mass)/sum(mass[1:nw])
taup=pen*sum(mass[nw+1:nage])/(w*nbar)  # balanced budet social security
taun = taulbar-taup
bequests=0

# endogenous aggregate variables (updated in outer loop)
glo = zeros(7)
glo[1] = kbar
glo[2] = nbar
glo[3] = pen
glo[4] = taup
glo[5] = taun
glo[6] = trbar
glo[7] = omega


# retired agents' value functions
vr1=zeros(na,nr)    # value function with lump-sum pensions: only depends on assets a
aropt1=zeros(na,nr) # optimal asset
cropt1=zeros(na,nr) # optimal consumption

function main(xagg,ef,sp1)
    kbar = xagg[1]
    nbar = xagg[2]
    pen = xagg[3]
    taup = xagg[4]
    taun = xagg[5]
    trbar = xagg[6]
    omega = xagg[7]
    glo1 = xagg

    na = na0
    nlabor = nlabor0
    dist=(kmax-kmin)/(na-1)
    a=collect(kmin:dist:kmax)
    dist=(labormax-labormin)/(nlabor-1)
    laborgrid=collect(labormin:dist:labormax)
    crit_counter=0
    tol1 = tol

    # plot survival probabilities
    periods=1:1:nage
    p1=plot(periods,sp1[1:nage],title="Survival Probabilities",lw=3,label="UN (2015)",)
    display(p1)

    # initial distribution of idiosyncratic productivities ye1 at age 0
    muy=zeros(ny)
    w=ye[2]-ye[1]
    muy[1]=cdf.(Normal(),((ye[1]+w/2)/sqrt(sigmay1)))
    muy[ny]=1-cdf.(Normal(),((ye[ny]-w/2)/sqrt(sigmay1)))
    for i in 2:ny-1
        muy[i] = cdf.(Normal(),((ye[i]+w/2)/sqrt(sigmay1)))
        muy[i] = muy[i]-cdf.(Normal(),((ye[i]-w/2)/sqrt(sigmay1)))
    end


    # cohort measures
    mass=ones(nage)
    for i in 2:nage
        mass[i] = mass[i-1]*sp1[i-1]/(1+popgrowth0)
    end
    mass = mass/sum(mass)
    age = 1:1:nage
    p=plot(age,mass,title="Measure of Households",lw=3,label="Year 2015",)
    display(p)
    #
    # compute Gini-coefficient of wage distribution
    #


    wages= zeros(nperm,ny,nw)   # productivity
    masswages = zeros(nperm,ny,nw) # measure
    # wages and measures at age 1
    for iperm in 1:nperm
        for iy in 1:ny
            wages[iperm,iy,1]=exp(ye[iy])*perm[iperm]*ef[1]
            masswages[iperm,iy,1]=1/2*muy[iy]*mass[1]/sum(mass[1:nw])
        end
    end

    for iage in 2:nw
        for iperm in 1:nperm
            for iy in 1:ny
                wages[iperm,iy,iage]=exp(ye[iy])*perm[iperm]*ef[iage]
                for iy1 in 1:ny
                    y0 = py[iy,iy1]*sp1[iage-1]/(1+popgrowth0)*masswages[iperm,iy,iage-1]
                    masswages[iperm,iy1,iage]= masswages[iperm,iy1,iage]+y0
                end
            end
        end
    end

    # reshape matrix as a vector
    wages1 = reshape(wages,:,1)
    masswages1 = reshape(masswages,:,1)
    wdist = hcat(wages1,masswages1)     # horizontal concatenation
    wdist1 = sortslices(wdist,dims=1)    # sort matrix according to first column
    gini_wage= get_gini(wdist1[:,1],wdist1[:,2])
    println("")
    println("Gini coefficient of wage disribution: $gini_wage")

    #
    # start outer iteration over aggregate variables
    #


    # saving the aggregate variables in outer iteration q
    # kbar and nbar

    aggregatesq = zeros((nq,7))
    q=0
    crit = 1.0 + tol # initial value for divergence of final solution
    while (q<nq) & (crit>tol)
    	q=q+1
    	w = wage(kbar,nbar)
    	r = interest(kbar,nbar)
        rb = (1-tauk)*(r-delta)
        ybar = production(kbar,nbar)
        debt = bybar*ybar
        gbar = gybar*ybar
        meanincome = w*nbar/sum(mass[1:nw]) # average earnings of workers
        penmin = peny*meanincome

        # save aggregate variables in outer iteration -> check convergence
        aggregatesq[q,:]=xagg


        # retired agents' value functions
        vr=zeros(na,nce,nr)    # value function with lump-sum pensions: only depends on assets a
        aropt=zeros(na,nce,nr) # optimal asset
        cropt=zeros(na,nce,nr) # optimal consumption

        for ia in 1:na
            for ice in 1:nce
                c=a[ia]*(1+(1-tauk)*(r-delta))+pension(e[ice],meanincome,penmin)+trbar
                c = c/(1+tauc)
                vr[ia,ice,nr]=u(c,0)
                cropt[ia,ice,nr]=c
            end
        end


        # computation of the decision rules for the retired

        for iage in nr-1:-1:1

            println("outer iteration q=$q and age $iage")

            for ia in 1:na
                a0 = a[ia]
                for ice in 1:nce
                    e0 = e[ice]
                    c=(1+rb)*a0+pension(e0,meanincome,penmin)+trbar.-ygrowth*a
                    c = c/(1+tauc)
                    c[c.<=0.0].=psi1
                    labor=zeros(size(c)[1],1);
                    y=u.(c,labor)
                    y=y+ygrowth^(gam*(1-eta))*sp1[nw+iage]*beta*vr[:,ice,iage+1]

                    indexopt_all=argmax(y);
                    indexopt=indexopt_all[1]
                    # display(indexopt)
                    k1 = a[indexopt];
                    aropt[ia,ice,iage] = k1
                    vr[ia,ice,iage] = y[indexopt] # maximum in bellman equations y
                    cropt[ia,ice,iage]=( (1+rb)*a0+pension(e0,meanincome,penmin)+trbar-ygrowth*k1 )/(1+tauc)
                end
            end

            if iage==1
                println("ia=10,ice=3, iage=1: v=")
                display(vr[10,3,1])
            end

        end # iage



        # workers' value function
        vw=zeros(nperm,ny,na,nce,nw)
        awopt=zeros(nperm,ny,na,nce,nw)
        lopt=zeros(nperm,ny,na,nce,nw)
        cwopt=zeros(nperm,ny,na,nce,nw)

        for iage in nw:-1:1
            println("outer iteration q=$q and age $iage")
            println("K=$kbar and crit=$crit")


            for iperm in 1:nperm
                for iy in 1:ny
                    for ia in 1:na
                        a0 = a[ia]
                        for ice in 1:nce
                            e0=e[ice]
                            labor1 = laborgrid'.*ones(na,nlabor);
                            # next-period cumulated earnings
                            eff = ef[iage]*perm[iperm]*exp(ye[iy]) # idiosyncratic efficiency
                            e1=e0*(iage-1)/iage.+laborgrid.*w.*eff./iage
                            c=(1+rb)*a0.+(1-taun-taup).*w.*eff.*labor1.+trbar.-ygrowth*a.*ones(na,nlabor)
                            c=c/(1+tauc)
                            c[c.<=0.0].=psi1
                            bellman = u.(c,labor1)

                            vtemp = zeros(na,nlabor)

                            if iage==nw
                                # contruction of the next period value function
                                v0 = zeros(na)

                                for ilabor in 1:nlabor
                                    # interpolation of value function at e1[ilabor]
                                    e10=e1[ilabor]
                                    i1=sum(e.<e10)
                                    i1=max(1,i1)
                                    i1=min(nce-1,i1)
                                    z1=(e10-e[i1])/(e[i1+1]-e[i1])
                                    v0=(1-z1).*vr[:,i1,1].+z1.*vr[:,i1+1,1]
                                    # add the value function for the labor supply lgrid[ilabor]
                                    vtemp[1:na,ilabor] = v0
                                end # ilabor
                            else   # iage<nw
                                # contruction of the next period value function
                                v0 = zeros(na)
                                for iy1 in 1:ny
                                    vtemp0 = zeros(na,nlabor)
                                    for ilabor in 1:nlabor
                                        e10=e1[ilabor]
                                        i1=sum(e.<e10)
                                        i1=max(1,i1)
                                        i1=min(nce-1,i1)
                                        z1=(e10-e[i1])/(e[i1+1]-e[i1])
                                        v0=(1-z1).*vw[iperm,iy1,:,i1,iage+1].+z1.*vw[iperm,iy1,:,i1+1,iage+1]
                                        # add the value function for the labor supply lgrid[ilabor]
                                        vtemp0[1:na,ilabor] = v0
                                    end # ilabor
                                    vtemp = vtemp + py[iy,iy1]*vtemp0
                                end # iy1 =1,..,ny

                            end # if iage==nw

                            discount_factor = ygrowth^(gam*(1-eta))*sp1[iage]*beta
                            bellman = bellman.+discount_factor.*vtemp
                            maxbellman=[maximum(bellman[:,i]) for i in 1:1:size(bellman)[2]];
                            # display(maxbellman)
                            ilabormax=argmax(maxbellman)[1];
                            maxbellman1=bellman[:,ilabormax];
                            iamax=argmax(maxbellman1)[1];
                            k1=a[iamax];
                            labor=laborgrid[ilabormax];
                            c = (1+rb)*a0+(1-taun-taup)*w*eff*labor+trbar-ygrowth*k1
                            c = c/(1+tauc)
                            vw[iperm,iy,ia,ice,iage]=bellman[iamax,ilabormax]
                            awopt[iperm,iy,ia,ice,iage]=k1
                            cwopt[iperm,iy,ia,ice,iage]=c
                            lopt[iperm,iy,ia,ice,iage]=labor

                        end # ice
                    end # ia
                end # iy
           end # iperm

           # I compared the result with those from the PYTHON and GAUSS code
           if iage==45 || iage==44 || iage==1
               display("iage = $iage:")
               display(vw[1,5,10,5,iage])
               display(awopt[1,5,10,5,iage])
               display(cwopt[1,5,10,5,iage])
               display(lopt[1,5,10,5,iage])
           end

        end  # iage


        save("vw.jld","vw",vw)
        save("awopt.jld","awopt",awopt)
        save("cwopt.jld","cwopt",cwopt)
        save("lopt.jld","lopt",lopt)

        save("vr.jld","vr",vr)
        save("aropt.jld","aropt",aropt)
        save("cropt.jld","cropt",cropt)
        save("aggregatesq.jld","aggregatesq",aggregatesq)

        # -----------------------------------------------------------------------
        #
        # computation of the stationary distribution
        #
        #-------------------------------------------------------------------------


        gkw=zeros(nperm,ny,na,nce,nw)      # distribution of wealth among workers
        gkr=zeros(na,nce,nr)               # distribution of wealth among retirees
        agen=zeros(nage)                # distribution of wealth over age
        fa=zeros(na)           # distribution of wealth
        fe = zeros(nce)             # distribution of accumulated earnings
        gincome=zeros(nin)          # distribution of income (all households)
        gearn = zeros(nin)
        gk=zeros(na,nce,nw)    # auxiliary variable

        bigl=0.0
        bigcontrib=0.0
        mean_labor = 0.0
        bigc = 0.0
        biga = 0.0
        bequest = 0.0
        bigpension = 0.0

        # mass at age 1
        # all agents have zero wealth and zero accumulated earnings
        for iy in 1:ny
            gkw[1,iy,1,1,1] = 1/2*muy[iy]*mass[1]      # measure of perm[1] at age 1
            gkw[2,iy,1,1,1] = 1/2*muy[iy]*mass[1]     # perm[2] at age 1
        end

        #
        # distribution at age 2,...,nw
        #
        for iage in 1:nw-1
            for iperm in 1:nperm
                for iy in 1:ny
                    for ia in 1:na
                        a0 = a[ia]
                        for ice in 1:nce
                            e0 = e[ice]
                            measure = gkw[iperm,iy,ia,ice,iage]
                            a1=awopt[iperm,iy,ia,ice,iage]   # optimal a' for a=0
                            labor = lopt[iperm,iy,ia,ice,iage]
                            c = cwopt[iperm,iy,ia,ice,iage]
                            agen[iage]=agen[iage]+measure*a0
                            fa[ia]=fa[ia]+measure		# wealth distribution
                            # average working hours
                            mean_labor = mean_labor + labor*measure/sum(mass[1:nw])
                            # aggregate consumption
                            bigc = bigc + c*measure
                            biga=biga+a0*measure
                            bequest = bequest+(1+rb)*a1*measure*(1-sp1[iage])/(1+popgrowth0)

                            #
                            # distribution of earnings and income at age iage
                            #
                            x = perm[iperm]*exp(ye[iy])*ef[iage]*w*labor   # earnings
                            y = x + (r-delta)*a0                # income

                            # adding the measure of (iperm,iy,ia,iage) to distribution of earnings/income
                            if x<=0
                                gearn[1] = gearn[1] +  measure/sum(mass[1:nw])
                            elseif x>=incomegrid[nin]
                                gearn[nin] = gearn[nin] + measure/sum(mass[1:nw])
                            else   # linear interpolation between grid points
                                j0=sum(incomegrid.<x)+1
                                j0=min(j0,nin)
                                lambda=(x-incomegrid[j0-1])/(incomegrid[j0]-incomegrid[j0-1])
                                gearn[j0-1] = gearn[j0-1]+ (1-lambda)*measure/sum(mass[1:nw])
                                gearn[j0] = gearn[j0]+ lambda*measure/sum(mass[1:nw])
                            end

                            if y<0
                                gincome[1] = gincome[1] + measure
                            elseif y>=incomegrid[nin]
                                gincome[nin] = gincome[nin] + measure
                            else
                                j0=sum(incomegrid.<y)+1
                                j0=min(j0,nin)
                                n0=(y-incomegrid[j0-1])/(incomegrid[j0]-incomegrid[j0-1])
                                gincome[j0-1] = gincome[j0-1]+ (1-n0)*measure
                                gincome[j0] = gincome[j0]+ n0*measure
                            end

                            #
                            # dynamics of the distribution function gkw
                            #
                            # iage -> iage+1
                            #


                            # next-period cumulated earnings
                            e1=(iage-1)/iage*e0+perm[iperm]*exp(ye[iy])*ef[iage]*w*labor/iage
                            # Summing up aggregate labor in efficiency unit (=L)
                            bigl=bigl+labor*perm[iperm]*ef[iage]*exp(ye[iy])*measure
                            # pension contributions
                            bigcontrib=bigcontrib+taup*labor*perm[iperm]*ef[iage]*exp(ye[iy])*w*measure


                            ia1=sum(a.<=a1)
                            ice1 = sum(e.<=e1)

                            if e1<=0
                                for iy1 in 1:ny
                                    gkw[iperm,iy1,ia1,1,iage+1]=(gkw[iperm,iy1,ia1,1,iage+1]+
                                        py[iy,iy1]*sp1[iage]/(1+popgrowth0)*measure)
                                end

                            elseif ice1>=nce
                                for iy1 in 1:ny
                                    gkw[iperm,iy1,ia1,nce,iage+1]=(gkw[iperm,iy1,ia1,nce,iage+1]+
                                        py[iy,iy1]*sp1[iage]/(1+popgrowth0)*measure)
                                end

                            else	# linear interpolation between the two adjacent grid points of e1 on grid e

                                lambda=(e[ice1+1]-e1) / (e[ice1+1]-e[ice1] )
                                for iy1 in 1:ny
                                    gkw[iperm,iy1,ia1,ice1,iage+1]=(gkw[iperm,iy1,ia1,ice1,iage+1]+
                                        lambda*py[iy,iy1]*sp1[iage]/(1+popgrowth0)*measure)
                                    gkw[iperm,iy1,ia1,ice1+1,iage+1]=(gkw[iperm,iy1,ia1,ice1+1,iage+1]+
                                        (1-lambda)*py[iy,iy1]*sp1[iage]/(1+popgrowth0)*measure)
                                end
                            end # endif: e' boundary solution
                        end # ice
                    end # a[ia], ia=1,..nag
                end   # ye[iy], j=1,..,ny
            end       #  iperm =1,2
        end       # for iage=1,..,nw-1 -> measure in iage+1

        #
        # distribution in the first year of retirement at age nw+1
        #

        iage=nw
        for iperm in 1:nperm
            for iy in 1:ny
                for ia in 1:na
                    a0 = a[ia]
                    for ice in 1:nce
                        e0 = e[ice]
                        measure = gkw[iperm,iy,ia,ice,iage]
                        a1=awopt[iperm,iy,ia,ice,iage]   # optimal a' for a=0
                        labor = lopt[iperm,iy,ia,ice,iage]
                        c = cwopt[iperm,iy,ia,ice,iage]
                        agen[iage]=agen[iage]+measure*a0
                        fa[ia]=fa[ia]+measure		# wealth distribution
                        # average working hours
                        mean_labor = mean_labor + labor*measure/sum(mass[1:nw])
                        # aggregate consumption
                        bigc = bigc + c*measure
                        biga=biga+a0*measure
                        bequest = bequest+(1+rb)*a1*measure*(1-sp1[iage])/(1+popgrowth0)

                        #
                        # distribution of earnings and income at age iage
                        #
                        x = perm[iperm]*exp(ye[iy])*ef[iage]*w*labor   # earnings
                        y = x + (r-delta)*a0                # income

                        # adding the measure of (iperm,iy,ia,iage) to distribution of earnings/income
                        if x<=0
                            gearn[1] = gearn[1] +  measure/sum(mass[1:nw])
                        elseif x>=incomegrid[nin]
                            gearn[nin] = gearn[nin] + measure/sum(mass[1:nw])
                        else   # linear interpolation between grid points
                            j0=sum(incomegrid.<x)+1
                            j0=min(j0,nin)
                            n0=(x-incomegrid[j0-1])/(incomegrid[j0]-incomegrid[j0-1])
                            gearn[j0-1] = gearn[j0-1]+ (1-n0)*measure/sum(mass[1:nw])
                            gearn[j0] = gearn[j0]+ n0*measure/sum(mass[1:nw])
                        end

                        if y<0
                            gincome[1] = gincome[1] + measure
                        elseif y>=incomegrid[nin]
                            gincome[nin] = gincome[nin] + measure
                        else
                            j0=sum(incomegrid.<y)+1
                            j0=min(j0,nin)
                            n0=(y-incomegrid[j0-1])/(incomegrid[j0]-incomegrid[j0-1])
                            gincome[j0-1] = gincome[j0-1]+ (1-n0)*measure
                            gincome[j0] = gincome[j0]+ n0*measure
                        end

                        #
                        # dynamics of the distribution function gkw
                        #
                        # iage -> iage+1
                        #


                        # next-period cumulated earnings
                        e1=(iage-1)/iage*e0+perm[iperm]*exp(ye[iy])*ef[iage]*w*labor/iage
                        # Summing up aggregate labor in efficiency unit (=L)
                        bigl=bigl+labor*perm[iperm]*ef[iage]*exp(ye[iy])*measure
                        # pension contributions
                        bigcontrib=bigcontrib+taup*labor*perm[iperm]*ef[iage]*exp(ye[iy])*w*measure


                        ia1=sum(a.<=a1)
                        ice1 = sum(e.<=e1)
                        if e1<=0
                            gkr[ia1,1,1]= gkr[ia1,1,1]+sp1[iage]/(1+popgrowth0)*measure
                        elseif ice>=nce
                            gkr[ia1,nce,1]= gkr[ia1,nce,1]+sp1[iage]/(1+popgrowth0)*measure
                        else	# linear interpolation between the two adjacent grid points of e1 on grid e
                            lambda=(e[ice1+1]-e1) / (e[ice1+1]-e[ice1] )
                            gkr[ia1,ice1,1]= gkr[ia1,ice1,1]+lambda*sp1[iage]/(1+popgrowth0)*measure
                            gkr[ia1,ice1+1,1]= gkr[ia1,ice1+1,1]+(1-lambda)*sp1[iage]/(1+popgrowth0)*measure
                        end # endif: e' boundary solution
                    end # ice
                end # a[ia], ia=1,..nag
            end   # ye[iy], j=1,..,ny
        end       #  iperm =1,2

        #
        # distribution of wealth among the retiree at age nw+2, nw+3,...
        #



        for iage in 1:nr
            for ia in 1:na
                a0 = a[ia]
                for ice in 1:nce
                    e0 = e[ice]
                    measure = gkr[ia,ice,iage]
                    if iage==1 # first year of retirement
                        fe[ice] = fe[ice]+measure/mass[nw+1]
                    end

                    x = pension(e0,meanincome,penmin)   # earnings
                    a1=aropt[ia,ice,iage]   # optimal a' for a=0
                    c = cropt[ia,ice,iage]
                    labor = 0.0
                    agen[iage+nw]=agen[iage+nw]+measure*a0
                    fa[ia]=fa[ia]+measure		# wealth distribution

                    # aggregate consumption
                    bigc = bigc + c*measure
                    biga=biga+a0*measure
                    bequest = bequest+(1+rb)*a1*measure*(1-sp1[iage+nw])/(1+popgrowth0)
                    bigpension = bigpension + x*measure

                    #
                    # distribution of earnings and income at age iage
                    #
                    y = x + (r-delta)*a0                # income

                    if y<0
                        gincome[1] = gincome[1] + measure
                    elseif y>=incomegrid[nin]
                        gincome[nin] = gincome[nin] + measure
                    else
                        j0=sum(incomegrid.<y)+1
                        j0=min(j0,nin)
                        n0=(y-incomegrid[j0-1])/(incomegrid[j0]-incomegrid[j0-1])
                        gincome[j0-1] = gincome[j0-1]+ (1-n0)*measure
                        gincome[j0] = gincome[j0]+ n0*measure
                    end

                    #
                    # dynamics of the distribution function gkr
                    #
                    # iage -> iage+1
                    #
                    # dynamics of the distribution during retirement
                    #
                    if iage<nr
                        ia1 = sum(a.<=a1)
                        gkr[ia1,ice,iage+1] = gkr[ia1,ice,iage+1] + sp1[iage+nw]/(1+popgrowth0)*measure
                    end     # iage<nr
                end # ice
            end     # ia
        end     # iage =1,..,nr

        totalmeasure = sum(sum(gkw))+sum(sum(gkr))
        println("total measure = $totalmeasure")
        #
        # computation of the Gini coefficients
        #    - wealth, earnings, income
        gini_wealth = get_gini(a,fa)
        println("gini wealth:  ")
        display(gini_wealth)
        gini_earnings = get_gini(incomegrid,gearn)
        println("gini earnings:")
        display(gini_earnings)
        println("gini income:")
        gini_income = get_gini(incomegrid,gincome)
        display(gini_income)

        taxes=taun*w*nbar+tauk*(r-delta)*kbar+tauc*bigc
        transfernew=taxes+bequest+debt*((1+popgrowth0)*ygrowth-(1+(1-tauk)*(r-delta))) - gbar
        trbar=phi*trbar+(1-phi)*transfernew

        # total savings
        omeganew = biga
        ybar=production(kbar,nbar)
        debt = bybar*ybar
        gbar = gybar*ybar
        omega = phi*omega + (1-phi)*omeganew
        knew = omeganew - debt
        kbar=phi*kbar+(1-phi)*knew
        nbar = phi*nbar + (1-phi)*bigl


        taupnew = bigpension/(w*nbar)
        taup = phi*taup + (1-phi)*taupnew
        taun = taulbar-taup    # calibration so that tau^n +tau^p = 28% as in the US

        # average pension / average wage income
        replacement_rate_new=bigpension/sum(mass[nw+1:nage])*sum(mass[1:nw])/(w*nbar)





        # check accuracy of policy functions:
        # Residual: Euler equation
        #
        Euler_res = zeros(nperm,ny,na,nce,nw)
        Euler_res_old = zeros(na,nce,nr-1)

        for iage in 1:nw
            for iperm in 1:nperm
                for iy in 1:ny
                    for ia in 1:na
                        a0 = a[ia]
                        for ice in 1:nce
                            e0 = e[ice]
                            a1=awopt[iperm,iy,ia,ice,iage]
                            labor=lopt[iperm,iy,ia,ice,iage]
                            c = cwopt[iperm,iy,ia,ice,iage]
                            e1=(iage-1)/iage*e0+perm[iperm]*exp(ye[iy])*ef[iage]*w*labor/iage

                            # computation of the Euler residual
                            #
                            # (1+g_A)^eta u_{c,t} = beta phi^i E_t{ u_{c,t+1} [1+(1-tauk) (r-delta)] }
                            #

                            x=0
                            if iage<nw
                                for iy1 in 1:ny
                                    x = x + beta*sp1[iage]*py[iy,iy1]* (1+rb)*uc1(a1,iperm,iy1,iage+1,cwopt,lopt,cropt,e1)
                                end
                            else
                                x = beta*sp1[iage]*(1+rb)*uc1(a1,iperm,iy,iage+1,cwopt,lopt,cropt,e1)
                            end
                            Euler_res[iperm,iy,ia,ice,iage] = 1-x / ( ygrowth^(1-gam*(1-eta)) * uc(c,labor))

                        end
                    end
                end
            end
        end

        for iage in 1:nr-1
            for ia in 1:na
                a0 = a[ia]
                for ice in 1:nce
                    e0 =e[ice]
                    a1=aropt[ia,ice,iage]
                    c = cropt[ia,ice,iage]
                    labor=0
                    Euler_res_old[ia,ice,iage] = 1- beta*sp1[nw+iage]*(1+rb)*uc1(a1,1,1,iage+nw+1,cwopt,lopt,cropt,e0)/ (ygrowth^(1-gam*(1-eta))*uc(c,labor))
                end
            end
        end


        println("mean Euler residual young: ")
        display(mean(abs.(Euler_res)))
        println("mean Euler residual old: ")
        display(mean(abs.(Euler_res_old)))


        agen1=zeros(nage)
        agen1[1:nage]=agen[1:nage]./mass[1:nage]
        age = 1:1:nage
        p=plot(age,agen1,title="Wealth-age Profile",lw=3,label="Year 2015",)
        display(p)

        crit = abs((kbar-knew)/knew)

        println("iteration q = $q")
        println("mean working hours:  $mean_labor")
        println("new capital stock: $kbar")
        println("new labor: $nbar")


        if (crit<tol1) & (crit_counter==0)
            crit = abs(crit)+tol
            crit_counter=1
            tol1 = tol_final
            na=na1        # grid numbers close to convergence of aggregate capital stock
            nlabor = nlabor1
            dist=(kmax-kmin)/(na-1)
            a=collect(kmin:dist:kmax)
            dist=(labormax-labormin)/(nlabor-1)
            laborgrid=collect(labormin:dist:labormax)
        end


    end     # outer while loop over aggregates, q = 1, .., nq



    return
end  # function main

@time main(glo,ef,sp1)
        #vr1 = main(glo,ef,sp1)
            #display(vr1[10,1])
