# AK70_stochastic_income_main.jl
#
# computes the steady in the OLG model (Auerbach-Kotlikoff)
# with stochastic idiosyncratic income
#
# Mthod: VALUE FUNCTION ITERATION
#
# in Heer/Maussner, "Dynamic General Equilibrium Modeling", 2009,
# third edition, Section 10.1.1
#
# Date: June 12, 2021
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
using Dierckx
include("AK70_stochastic_income_procs.jl");
#import Distributions: pdf, Normal, quantile

#**************
# output title:
#**************
println("")
println("")
println("")
println("")
println("*******************************************")
println("******* Value Function Iteration: *********")
println("*******************************************")
println("******* Heer/Maussner: DGE Modeling *******")
println("******* Section 10.1.1 (3rd edition) ******")
println("******* Auerbach-Kotlikoff Model **********")
println("******* STOCHASTIC INDIVIDUAL INCOME ******")
println("*******************************************")
println("")


#************
# Parameters:
#************

# numerical parameters
const ipolt = 2     # interpolation of policy functions:
                    # 1=linear or 2=Cubic

const scalf=-1.0     # muliplier of rhs of Bellman eq to find MINIMUM instead of MAXIMUM
                     # in golden section search
const eps0=0.05       # small parameter to check if we have boundary solution for a'(a) at a=0 or a=kmax
const phi=0.80       # update aggregate variables in outer iteration over K, L, tr, taup, taun
const tol=0.0001     # percentage deviation of final solution K and L
const tol1=1e-5      # tolerance for golden section search
const neg=-1e10      # initial value for value function
const nq = 50        # number of outer iterations

# constants for interpolations
if ipolt==1
    const kx_pf=1
elseif ipolt==2
    const kx_pf=3
end
const s_pf=0.0

# asset and income grids
const na=501    # number of grid points on asset grid for value + policy function
const nag=1002  # number of grid points on asset grid for distribution
const kmin=0.0  # minimum asset level
const kmax=20.0 # maximum asset level
const a = range(kmin,stop=kmax,length=na)
const ag = range(kmin,stop=kmax,length=nag)
const nearn=501 # number of grid points on earnings/income grid -> gini coefficent
const earnmin=0.0   # minimum of earnings/income
const earnmax=5.0   # maximum of earnings/income
const earn = range(earnmin,stop=earnmax,length=nearn)
const income = range(earnmin,stop=earnmax,length=nearn)
const labormax=0.60 # upper bound on labor supply



# model:
const beta=1.011     # discount factor
const eta=2.00      # 1/IES (intertemporal elasticity of substitution)
const alp=0.35      # production elasticity of capital
const delta=0.083     #  rate of depreciation
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
const replacement_ratio=0.352	# gross replacement ratio US
const bybar=0.63	# debt-output ratio
const gybar=0.18	# government consumption-output ratio


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
# initial guess for dependency ratio 1/3, therefore share of workers in
# population equal to 2/3
pen=replacement_ratio*(1-taulbar)*w*mean_labor*3/2
taup=1/3*pen/(w*nbar)  # balanced budet social security
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
    krit_agg = 1.0 + tol # initial value for divergence of final solution
    while (q<nq) & (krit_agg>tol)
    	q=q+1
    	w = wage(kbar,nbar)
    	r = interest(kbar,nbar)
        # save aggregate variables in outer iteration -> check convergence
        aggregatesq[q,:]=xagg
        kold=kbar
        nold=nbar
        rb=(1-tauk)*(r-delta)

        # retired agents' value functions
        vr=zeros(na,nr)    # value function with lump-sum pensions: only depends on assets a
        aropt=zeros(na,nr) # optimal asset
        cropt=zeros(na,nr) # optimal consumption

        for ia in 1:na
            c=a[ia]*(1+(1-tauk)*(r-delta))+pen+trbar
            c = c/(1+tauc)
            vr[ia,nr]=u(c,0)
            cropt[ia,nr]=c
        end


        # computation of the decision rules for the retired

        for age in nr-1:-1:1
            # preparation of interplation for  next-period value function
            vfa_grid = vr[:,age+1]

            println("outer iteration q=$q and age $age")

            # prepare interpolation
            itpu=Spline1D(a,vfa_grid,k=kx_pf,s=s_pf)
            itp=(k1)->evaluate(itpu,k1)

            glo2 = zeros(4)
            glo2[2] = r
            glo2[3] = pen
            glo2[4] = trbar

            for ia in 1:na
                a0 = a[ia]
                glo2[1] = a0

                # boundary points GOLDEN SECTION Search ax,bx
                # search for [a[ia],a[ia+1]] = [ax,bx] that bracket maximum
                ax=0.0
                bx=-1.0
                v0=neg  # initial value value function
                # search for ax,bx such that value function declines again
                ia1=0
                while (ax>bx)
                    ia1 = ia1+1
                    # evaluate value function at a[ia1]
                    v1 = value1(a[ia1],glo2,itp,sp1,age)
                    if v1>v0
                        if ia1==1
                            ax=a[ia1]
                        else
                            ax=a[ia1-1]
                        end
                        v0=v1
                    else
                        bx=a[ia1]
                    end     # search over boundaries golden section search
                            # that bracket the maximum of the rhs of Bellman eq.
                end # do while loop to find [ax,bx] that brackets maximum

                if ax==a[1]  # check: a[1] is maximum point on grid?
                    cx=ax+(a[2]-a[1])*eps0
                    if value1(ax,glo2,itp,sp1,age)>value1(cx,glo2,itp,sp1,age) # boundary solution
                        a1max=a[1]
                        vmax = value1(a1max,glo2,itp,sp1,age)
                    else
                        bx=a[2]
                        res=optimize(x->scalf*value1(x,glo2,itp,sp1,age),ax,bx, GoldenSection(),abs_tol=tol1,rel_tol=tol1)
                        # optimal next-period capital stock a1 and value of the value function
                        # from Golden Section Search
                        a1max=Optim.minimizer(res)
                        vmax=Optim.minimum(res)*scalf
                    end
                elseif bx==a[na] # check: a[na] is maximum point on grid?
                    cx=a[na]-eps0*(a[na]-a[na-1])
                    if value1(bx,glo2,itp,sp1)>value1(cx,glo2,itp,sp1,age)
                        a1max=a[na]
                        vmax = value1(bx,glo2,itp,sp1,age)
                    else
                        ax=a[na-1]
                        res=optimize(x->scalf*value1(x,glo2,itp,sp1,age),ax,bx, GoldenSection(),abs_tol=tol1,rel_tol=tol1)
                        # optimal next-period capital stock a1 and value of the value function
                        # from Golden Section Search
                        a1max=Optim.minimizer(res)
                        vmax=Optim.minimum(res)*scalf
                    end
                else # interior solution ax<bx<cx
                    res=optimize(x->scalf*value1(x,glo2,itp,sp1,age),ax,bx, GoldenSection(),abs_tol=tol1,rel_tol=tol1)
                    # optimal next-period capital stock a1 and value of the value function
                    # from Golden Section Search
                    a1max=Optim.minimizer(res)
                    vmax =Optim.minimum(res)*scalf
                end

                aropt[ia,age] = a1max
                vr[ia,age] = vmax
                c = (1+(1-tauk)*(r-delta))*a0+pen+trbar-ygrowth*a1max
                cropt[ia,age] = c/(1+tauc)
            end # a[ia], l=1,..na
        end # age = nr, nr-1,...,1


        # In case you want to analyse some local variables from the function
        # main() (or any other) at the JULIA prompt
        # save("data.jld", "data", pfka)
        # -> analyse as follows: julia>load("data.jld")["data"]
        # -> julia>data
        save("vr.jld","vr",vr)
        save("aropt.jld","aropt",aropt)
        save("cropt.jld","cropt",cropt)
        save("agrid.jld","agrid",a)


        # -----------------------------------------------------------------------
        #
        # compuation of the decsion rules for the worker
        #
        #-------------------------------------------------------------------------


        # workers' value function
        vw=zeros(nperm,ny,na,nw)
        awopt=zeros(nperm,ny,na,nw)
        lopt=zeros(nperm,ny,na,nw)
        cwopt=zeros(nperm,ny,na,nw)

        for age in nw:-1:1  # all ages age=nw,nw-1,..1
             for iperm in 1:2 # all permanent productivity types
                #
                # prepare interpolation of value function at age "iage+1"
                # for use in the Bellman eq. at age "iage"
                # which needs to be passed to value2(.)
                #
                # Notice that we only need to prepare interpolation once at the
                # beginning of the iteration over age "iage" and "iperm" to save
                # computational time
                #
                # to pass on the interpolation nodes is a bit tricky
                # using the Dierckx-package
                #

                if age==nw
                    vfo=vr[:,1]
                    itpu=Spline1D(a,vfo,k=kx_pf,s=s_pf)
                    # fv[iy1]=itpu
                    itpo=(k1)->evaluate(itpu,k1)
                else
                    vfo=vw[iperm,:,:,age+1]
                    fv=Array{Spline1D}(undef,size(ye)[1],length(ye))
                    for iy1 in 1:ny
                        itpu=Spline1D(a,vfo[iy1,:],k=kx_pf,s=s_pf)
                        fv[iy1]=itpu
                    end
                    # define function that evaluates
                    # the value function spline at
                    # idiosyncratic efficiency efi
                    # and wealth level k1
                    itpo=(k1,efi)->evaluate(fv[efi],k1)

                    # test interpolation at a specific value
                    # println("Test evaluation of spline")
                    # println("interpolation value function at age $age")
                    # println("and idiosyncratic productivity 3")
                    # println("and wealth a=5.13")
                    # y = evaluate(fv[3],5.18)
                    # display(y)
                    # println("true value")
                    # y = (vw[iperm,3,130,age+1]+vw[iperm,3,131,age+1])/2.0
                    # display(y)
                end

                Threads.@threads for iy in 1:ny # all idiosyncratic productivity types at age "age"
                    ia10 = 0    # start search over optimal a' at a[ia10+1]
                    vglo = zeros(6)
                    vglo[2] = r
                    vglo[3] = w
                    vglo[4] = trbar
                    vglo[5] = taun
                    vglo[6] = taup

                    for ia in 1:na # all wealth holdings at age "age"
                        # boundary points GOLDEN SECTION Search ax,bx
                        # search for [a[ia],a[ia+1]] = [ax,bx] that bracket maximum
                        a0 = a[ia]
                        vglo[1] = a0
                        ax=0.0
                        bx=-1.0
                        v0=neg  # initial value value function
                        # search for ax,bx such that value function declines again
                        ia1=ia10
                        while (ax>bx)
                            ia1 = ia1+1

                            # evaluate value function at a[ia1]
                            v1 = value2(a[ia1],vglo,sp1,ye,py,ef,age,iperm,iy,itpo)
                            if v1>v0
                                ia0 = max(0,ia1-2)  # monotonocity of the value function
                                                    # with respect to a'(a)a[m];
                                if ia1==1
                                    ax=a[ia1]
                                else
                                    ax=a[ia1-1]
                                end
                                v0=v1
                            else
                                bx=a[ia1]
                            end     # search over boundaries golden section search
                                    # that bracket the maximum of the rhs of Bellman eq.
                            if ia1==na
                                ax = a[na-1]
                                bx = a[na]
                            end
                        end # do while loop to find [ax,bx] that brackets maximum

                        if ax==a[1]  # check: a[1] is maximum point on grid?
                            cx=ax+(a[2]-a[1])*eps0
                            if value2(ax,vglo,sp1,ye,py,ef,age,iperm,iy,itpo)>value2(cx,vglo,sp1,ye,py,ef,age,iperm,iy,itpo) # boundary solution
                                a1max=a[1]
                                vmax = value2(a1max,vglo,sp1,ye,py,ef,age,iperm,iy,itpo)
                            else
                                bx=a[2]
                                res=optimize(x->scalf*value2(x,vglo,sp1,ye,py,ef,age,iperm,iy,itpo),ax,bx, GoldenSection(),abs_tol=tol1,rel_tol=tol1)
                                # optimal next-period capital stock a1 and value of the value function
                                # from Golden Section Search
                                a1max=Optim.minimizer(res)
                                vmax=Optim.minimum(res)*scalf
                            end
                        elseif bx==a[na] # check: a[na] is maximum point on grid?
                            cx=a[na]-eps0*(a[na]-a[na-1])
                            if value2(bx,vglo,sp1,ye,py,ef,age,iperm,iy,itpo)>value2(cx,vglo,sp1,ye,py,ef,age,iperm,iy,itpo)
                                a1max=a[na]
                                vmax = value2(bx,vglo,sp1,ye,py,ef,age,iperm,iy,itpo)
                            else
                                ax=a[na-1]
                                res=optimize(x->scalf*value2(x,vglo,sp1,ye,py,ef,age,iperm,iy,itpo),ax,bx, GoldenSection(),abs_tol=tol1,rel_tol=tol1)
                                # optimal next-period capital stock a1 and value of the value function
                                # from Golden Section Search
                                a1max=Optim.minimizer(res)
                                vmax=Optim.minimum(res)*scalf
                            end
                        else # interior solution ax<bx<cx
                            res=optimize(x->scalf*value2(x,vglo,sp1,ye,py,ef,age,iperm,iy,itpo),ax,bx, GoldenSection(),abs_tol=tol1,rel_tol=tol1)
                            # optimal next-period capital stock a1 and value of the value function
                            # from Golden Section Search
                            a1max=Optim.minimizer(res)
                            vmax =Optim.minimum(res)*scalf
                        end

                        awopt[iperm,iy,ia,age] = a1max
                        vw[iperm,iy,ia,age] = vmax

                        lglo = vcat(a1max,vglo)
                        labor = optimal_labor(lglo,ye,ef,age,iperm,iy)
                        lopt[iperm,iy,ia,age] = labor

                        eff = ef[age]*perm[iperm]*exp(ye[iy]) # idiosyncratic efficiency
                        c=(1+(1-tauk)*(r-delta))*a0+(1-taun-taup)*w*eff*labor+trbar-ygrowth*a1max
                        cwopt[iperm,iy,ia,age] = c/(1+tauc)



                    end   # a: l=1,..,na
                end  # idiosyncratic productivity j=1,..,ny
            end       # iperm = 1,2
            println("computation of value function")
            println("iteration q=$q ,age=$age")
            println("K= $kbar , L= $nbar")
        end           # age i=1,..,nw


        save("vw.jld","vw",vw)
        save("awopt.jld","awopt",awopt)
        save("cwopt.jld","cwopt",cwopt)
        save("lopt.jld","lopt",lopt)


        # -----------------------------------------------------------------------
        #
        # computation of the stationary distribution
        #
        #-------------------------------------------------------------------------


        gkw=zeros(nperm,ny,nag,nw)      # distribution of wealth among workers
        gkr=zeros(nag,nr)               # distribution of wealth among retirees
        gk=zeros(nag,nw)                # distribution of wealth at ages 1,..,nw
        kgen=zeros(nage)                # distribution of wealth over age
        gwealth=zeros(nag)              # distribution of wealth
        gearn=zeros(nearn)              # distribution of earnings (workers)
        gincome=zeros(nearn)            # distribution of income (all households)
        gk=zeros(nag,nw)    # auxiliary variable
        gk[1,1]=mass[1]

        # mass at age 1
        # all agents have zero wealth
        for iy in 1:ny
            gkw[1,iy,1,1] = 1/2*muy[iy]*mass[1]      # measure of perm[1] at age 1
            gkw[2,iy,1,1] = 1/2*muy[iy]*mass[1]     # perm[2] at age 1
        end

        #
        # distribution at age 2,...,nw
        #
        for iage in 1:nw-1
            for iperm in 1:nperm
                for iy in 1:ny
                    for iag in 1:nag
                        a0 = ag[iag]
                        if a0<=kmin
                            a1=awopt[iperm,iy,1,iage]   # optimal a' for a=0
                            labor = lopt[iperm,iy,1,iage]
                        elseif a0>=kmax
                            a1=awopt[iperm,iy,na,iage] # optimal a' for a=kmax
                            labor = lopt[iperm,iy,na,iage]
                        else   # linear interpolation between grid points
                            j0=sum(a.<a0)+1
                            j0=min(j0,na)
                            n0=(a0-a[j0-1])/(a[j0]-a[j0-1])
                            a1=(1-n0)*awopt[iperm,iy,j0-1,iage] + n0*awopt[iperm,iy,j0,iage]
                            labor=(1-n0)*lopt[iperm,iy,j0-1,iage] + n0*lopt[iperm,iy,j0,iage];
                        end

                        #
                        # distribution of earnings and income at age i
                        #
                        x = perm[iperm]*exp(ye[iy])*ef[iage]*w*labor   # earnings
                        y = x + (r-delta)*a0                # income

                        # adding the measure of (iperm,iy,ia,iage) to distribution of earnings/income
                        if x<=0
                            gearn[1] = gearn[1] +  gkw[iperm,iy,iag,iage]/sum(mass[1:nw])
                        elseif x>=earn[nearn]
                            gearn[nearn] = gearn[nearn] + gkw[iperm,iy,iag,iage]/sum(mass[1:nw])
                        else   # linear interpolation between grid points
                            j0=sum(earn.<x)+1
                            j0=min(j0,nearn)
                            n0=(x-earn[j0-1])/(earn[j0]-earn[j0-1])
                            gearn[j0-1] = gearn[j0-1]+ (1-n0)*gkw[iperm,iy,iag,iage]/sum(mass[1:nw])
                            gearn[j0] = gearn[j0]+ n0*gkw[iperm,iy,iag,iage]/sum(mass[1:nw])
                        end

                        if y<0
                            gincome[1] = gincome[1] + gkw[iperm,iy,iag,iage]
                        elseif y>=income[nearn]
                            gincome[nearn] = gincome[nearn] + gkw[iperm,iy,iag,iage]
                        else
                            j0=sum(income.<y)+1
                            j0=min(j0,nearn)
                            n0=(y-earn[j0-1])/(earn[j0]-earn[j0-1])
                            gincome[j0-1] = gincome[j0-1]+ (1-n0)*gkw[iperm,iy,iag,iage]
                            gincome[j0] = gincome[j0]+ n0*gkw[iperm,iy,iag,iage]
                        end

                        #
                        # dynamics of the distribution function gkw
                        #
                        # iage -> iage+1
                        #
                        if a1<=kmin
                            for iy1 in 1:ny
                                gkw[iperm,iy1,1,iage+1]=(gkw[iperm,iy1,1,iage+1]
                                    +py[iy,iy1]*sp1[iage]/(1+popgrowth0)*gkw[iperm,iy,iag,iage])
                            end

                        elseif a1>=kmax
                            for iy1 in 1:ny
                                gkw[iperm,iy1,nag,iage+1] = (gkw[iperm,iy1,nag,iage+1]
                                    +py[iy,iy1]*sp1[iage]/(1+popgrowth0)*gkw[iperm,iy,iag,iage] )
                            end
                        else
                            j0=sum(ag.<a1)+1
                            j0=min(j0,nag)
                            n0=(a1-ag[j0-1])/(ag[j0]-ag[j0-1]);
                            for iy1 in 1:ny
                                gkw[iperm,iy1,j0-1,iage+1]= (gkw[iperm,iy1,j0-1,iage+1]
                                    +(1-n0)*py[iy,iy1]*sp1[iage]/(1+popgrowth0)*gkw[iperm,iy,iag,iage])
                                gkw[iperm,iy1,j0,iage+1] = ( gkw[iperm,iy1,j0,iage+1]
                                    +n0*py[iy,iy1]*sp1[iage]/(1+popgrowth0)*gkw[iperm,iy,iag,iage] )
                            end
                        end # endif: a' boundary solution
                    end # a[ia], ia=1,..nag
                end   # ye[iy], j=1,..,ny
            end       #  iperm =1,2

            # summing up the entries for ag[ia], ia=1,...,nag
            for iperm1 in 1:nperm
                for iy1 in 1:ny
                    meas0 = zeros(eltype(0.100),nag)
                    meas0[1:nag] = gkw[iperm1,iy1,1:nag,iage+1]
                    gk[1:nag,iage+1] = gk[1:nag,iage+1] + meas0 # distribution of wealth at iage+1
                    kgen[iage+1] = kgen[iage+1] + ag'*meas0 # wealth at age iage+1
                end
            end

        end       # for iage=1,..,nw-1 -> measure in iage+1

        #
        # distribution in the first year of retirement at age nw+1
        #

        iage=nw
        for iperm in 1:nperm
            for iy in 1:ny
                for iag in 1:nag
                    a0 = ag[iag]
                    if a0<=kmin
                        a1=awopt[iperm,iy,1,iage]
                        labor=lopt[iperm,iy,1,iage]
                    elseif a0>=kmax
                        a1=awopt[iperm,iy,na,iage]
                        labor=lopt[iperm,iy,na,iage]
                    else
                        j0=sum(a.<a0)+1
                        j0=min(j0,na)
                        n0=(a0-a[j0-1])/(a[j0]-a[j0-1])
                        a1 = (1-n0)*awopt[iperm,iy,j0-1,iage]+n0*awopt[iperm,iy,j0,iage]
                        labor = (1-n0)*lopt[iperm,iy,j0-1,iage]+n0*lopt[iperm,iy,j0,iage]
                    end


                    #
                    # distribution of earnings and income at age iage
                    #
                    x = perm[iperm]*exp(ye[iy])*ef[iage]*w*labor   # earnings
                    y = x + (r-delta)*a0                # income

                    if x<=0
                        gearn[1] = gearn[1] +  gkw[iperm,iy,iag,iage]/sum(mass[1:nw])
                    elseif x>=earn[nearn]
                        gearn[nearn] = gearn[nearn] + gkw[iperm,iy,iag,iage]/sum(mass[1:nw])
                    else   # linear interpolation between grid points
                        j0=sum(earn.<x)+1
                        j0=min(j0,nearn)
                        n0=(x-earn[j0-1])/(earn[j0]-earn[j0-1])
                        gearn[j0-1] = gearn[j0-1]+ (1-n0)*gkw[iperm,iy,iag,iage]/sum(mass[1:nw])
                        gearn[j0] = gearn[j0]+ n0*gkw[iperm,iy,iag,iage]/sum(mass[1:nw])
                    end

                    if y<0
                        gincome[1] = gincome[1] + gkw[iperm,iy,iag,iage]
                    elseif y>=income[nearn]
                        gincome[nearn] = gincome[nearn] + gkw[iperm,iy,iag,iage]
                    else
                        j0=sum(income.<y)+1
                        j0=min(j0,nearn)
                        n0=(y-earn[j0-1])/(earn[j0]-earn[j0-1])
                        gincome[j0-1] = gincome[j0-1]+ (1-n0)*gkw[iperm,iy,iag,iage]
                        gincome[j0] = gincome[j0]+ n0*gkw[iperm,iy,iag,iage]
                    end


                    if a1<=kmin
                        gkr[1,1]=gkr[1,1]+sp1[iage]/(1+popgrowth0)*gkw[iperm,iy,iag,iage]
                    elseif a1>=kmax
                        gkr[nag,1]=gkr[nag,1]+sp1[iage]/(1+popgrowth0)*gkw[iperm,iy,iag,iage]
                    else
                        j0=sum(ag.<a1)+1
                        n0=(a1-ag[j0-1])/(ag[j0]-ag[j0-1])
                        gkr[j0-1,1]=gkr[j0-1,1]+(1-n0)*sp1[iage]/(1+popgrowth0)*gkw[iperm,iy,iag,iage]
                        gkr[j0,1]=gkr[j0,1]+n0*sp1[iage]/(1+popgrowth0)*gkw[iperm,iy,iag,iage]
                    end
                end  # iag = 1,..,nag
            end     # iy = 1, ..ny
        end   # iperm=1,2

        meas0 = gkr[:,1]
        kgen[nw+1] = meas0'*ag

        #
        # distribution of wealth among the retiree at age nw+2, nw+3,...
        #

        for iage in 1:nr
            for iag in 1:nag
                a0 = ag[iag]
                if a0<=kmin
                    a1=aropt[1,iage]
                elseif a0>=kmax
                    a1=aropt[na,iage]
                else
                    j0=sum(a.<a0)+1
                    j0=min(j0,na)
                    n0=(a0-a[j0-1])/(a[j0]-a[j0-1])
                    a1=(1-n0)*aropt[j0-1,iage]+n0*aropt[j0,iage]
                end

                #
                # distribution of income
                #
                y = pen + (r-delta)*a0       # income

                if y<0
                    gincome[1] = gincome[1] + gkr[iag,iage]
                elseif y>=income[nearn]
                    gincome[nearn] = gincome[nearn] + gkr[iag,iage]
                else
                    j0=sum(income.<y)+1
                    j0=min(j0,nearn)
                    n0=(y-earn[j0-1])/(earn[j0]-earn[j0-1])
                    gincome[j0-1] = gincome[j0-1]+ (1-n0)*gkr[iag,iage]
                    gincome[j0] = gincome[j0] + n0*gkr[iag,iage]
                end

                if iage<nr
                    #
                    # dynamics of the distribution during retirement
                    #
                    if a1<=kmin
                        gkr[1,iage+1]=gkr[1,iage+1]+sp1[iage+nw]/(1+popgrowth0)*gkr[iag,iage]
                    elseif a1>=kmax
                        gkr[nag,iage+1]=gkr[nag,iage+1]+sp1[iage+nw]/(1+popgrowth0)*gkr[iag,iage]
                    else
                        j0=sum(ag.<a1)+1
                        j0=min(j0,nag)
                        n0=(a1-ag[j0-1])/(ag[j0]-ag[j0-1])
                        gkr[j0-1,iage+1]=gkr[j0-1,iage+1]+(1-n0)*sp1[iage+nw]/(1+popgrowth0)*gkr[iag,iage]
                        gkr[j0,iage+1]=gkr[j0,iage+1]+n0*sp1[iage+nw]/(1+popgrowth0)*gkr[iag,iage]
                    end
                    meas0 = gkr[1:nag,iage+1]
                    kgen[iage+nw+1]= meas0'*ag
                end     # iage<nr
            end     # iag
        end     # iage =1,..,nr


        #
        # computation of the Gini coefficients
        #    - wealth, earnings, income

        gk=hcat(gk,gkr)
        gk=gk/sum(gk)  # make sure that the sum of the measure is equal to one
                        # avoid possible rounding errors
        # computation of the gini coefficient of capital distribution
        # first compute horizontal sums of gk to get measures of ag[iag]
        gk1 = zeros(eltype(0.10),nag)
        gk1[1:nag]=sum(gk,dims=2)
        gini_wealth = get_gini(ag,gk1)
        println("gini wealth:  ")
        display(gini_wealth)

        gini_earnings = get_gini(earn,gearn)
        println("gini earnings:")
        display(gini_earnings)
        println("gini income:")
        gini_income = get_gini(earn,gincome)
        display(gini_income)


        # total savings
        omeganew=sum(kgen)/sum(mass)
        ybar=production(kbar,nbar)
        debt = bybar*ybar
        gbar = gybar*ybar
        omega = phi*omega + (1-phi)*omeganew
        knew = omeganew - debt
        kbar=phi*kold+(1-phi)*knew

        # aggregate variablec L, Beq, C, mean working hours mean_labor
        nnew = 0   # aggregate efficient labor L
        bequests = 0 # aggregate bequests
        bigc = 0   # aggregate consumption
        mean_labor=0  # average working hours -> input into pensions formula
        cgen = zeros(nage)  # consumption-age profile
        lgen = zeros(nw)    # labor supply-age profile

        # check accuracy of policy functions:
        # Residual: Euler equation
        #
        Euler_res = zeros(nperm,ny,nag,nw)
        Euler_res_old = zeros(nag,nr-1)

        for iage in 1:nw
            for iperm in 1:nperm
                for iy in 1:ny
                    for iag in 1:nag
                        a0 = ag[iag]
                        if a0<=kmin
                            a1=awopt[iperm,iy,1,iage]
                            labor=lopt[iperm,iy,1,iage]
                            c = cwopt[iperm,iy,1,iage]
                        elseif a0>=kmax
                            a1=awopt[iperm,iy,na,iage]
                            labor=lopt[iperm,iy,na,iage]
                            c = cwopt[iperm,iy,na,iage]
                        else   # linear interpolation between grid points
                            j0=sum(a.<a0)+1;
                            j0=min(j0,na)
                            n0=(a0-a[j0-1])/(a[j0]-a[j0-1])
                            a1=(1-n0)*awopt[iperm,iy,j0-1,iage]+n0*awopt[iperm,iy,j0,iage]
                            labor=(1-n0)*lopt[iperm,iy,j0-1,iage]+n0*lopt[iperm,iy,j0,iage]
                            c = (1-n0)*cwopt[iperm,iy,j0-1,iage]+n0*cwopt[iperm,iy,j0,iage]
                        end

                        nnew = nnew + ef[iage]*exp(ye[iy])*perm[iperm]*labor*gkw[iperm,iy,iag,iage]
                        mean_labor = mean_labor + labor/(sum(mass[1:nw]))*gkw[iperm,iy,iag,iage]
                        cgen[iage] = cgen[iage] + c*gkw[iperm,iy,iag,iage]/mass[iage]
                        lgen[iage]= lgen[iage] + labor*gkw[iperm,iy,iag,iage]/mass[iage]
                        bequests = bequests + a1*(1+(1-tauk)*(r-delta))*(1-sp1[iage])/(1+popgrowth0)*gkw[iperm,iy,iag,iage]
                        bigc = bigc + c*gkw[iperm,iy,iag,iage]
                        gwealth[iag] = gwealth[iag]+gkw[iperm,iy,iag,iage]

                        # computation of the Euler residual
                        #
                        # (1+g_A)^eta u_{c,t} = beta phi^i E_t{ u_{c,t+1} [1+(1-tauk) (r-delta)] }
                        #

                        x=0
                        if iage<nw
                            for iy1 in 1:ny
                                x = x + beta*sp1[iage]*py[iy,iy1]* (1+rb)*uc1(a1,iperm,iy1,iage+1,cwopt,lopt,cropt)
                            end
                        else
                            x = beta*sp1[iage]*(1+rb)*uc1(a1,iperm,iy,iage+1,cwopt,lopt,cropt)
                        end
                        Euler_res[iperm,iy,iag,iage] = 1-x / ( ygrowth^(1-gam*(1-eta)) * uc(c,labor))

                    end
                end
            end
        end

        for iage in 1:nr
            for iag in 1:nag
                a0 = ag[iag]
                if a0<=kmin
                    a1=aropt[1,iage]
                    c = cropt[1,iage]
                elseif a0>=kmax
                    a1=aropt[na,iage]
                    c = cropt[na,iage]
                else
                    j0=sum(a.<a0)+1
                    j0=min(j0,na)
                    n0=(a0-a[j0-1])/(a[j0]-a[j0-1])
                    a1=(1-n0)*aropt[j0-1,iage]+n0*aropt[j0,iage]
                    c=(1-n0)*cropt[j0-1,iage]+n0*cropt[j0,iage]
                end

                labor=0
                if iage<nr # there is no Euler eq in the last period of life: all income is consumed
                    Euler_res_old[iag,iage] = 1- beta*sp1[nw+iage]*(1+rb)*uc1(a1,1,1,iage+nw+1,cwopt,lopt,cropt)/ (ygrowth^(1-gam*(1-eta))*uc(c,labor))
                end

                bequests = bequests + a1*(1+(1-tauk)*(r-delta))*(1-sp1[iage+nw])/(1+popgrowth0)*gkr[iag,iage]
                bigc = bigc + c*gkr[iag,iage]
                cgen[nw+iage] = cgen[nw+iage]+ c*gkr[iag,iage]/mass[nw+iage]
            end
        end


        println("mean Euler residual young: ")
        display(mean(abs.(Euler_res)))
        println("mean Euler residual old: ")
        display(mean(abs.(Euler_res_old)))


        kgen1=zeros(nage)
        kgen1[1:nage]=kgen[1:nage]./mass[1:nage]
        age = 1:1:nage
        p=plot(age,kgen1,title="Wealth-age Profile",lw=3,label="Year 2015",)
        display(p)

        p=plot(age,cgen,title="Consumption-age Profile",lw=3,label="Year 2015",)
        display(p)

        age = 1:1:nw
        p=plot(age,lgen,title="Labor-supply age Profile",lw=3,label="Year 2015",)
        display(p)


        # tatonnement process to update L: simple dampening iterative scheme as described in Judd (1998), Section 3.9
        nbar=phi*nbar+(1-phi)*nnew
        taxes=taun*w*nbar+tauk*(r-delta)*kbar+tauc*bigc
        # update of L, K, transfer, kshare, and pen
        transfernew=taxes+bequests+debt*((1+popgrowth0)*ygrowth-(1+(1-tauk)*(r-delta))) - gbar
        trbar=phi*trbar+(1-phi)*transfernew
        pennew=replacement_ratio*w*mean_labor
        # social security contributions are calculated so that the social security budget balances
        taupnew=pennew*sum(mass[nw+1:nage])/(w*nbar)
        taup=phi*taup+(1-phi)*taupnew
        taunnew=taulbar-taup
        taun = phi*taun + (1-phi)*taunnew
        pen = phi*pen + (1-phi)*pennew

        krit_agg1 = abs((kbar-kold)/kold)
        krit_agg2 = abs((nbar-nold)/nold)
        krit_agg = max(krit_agg1,krit_agg2)


        println("iteration q = $q")
        println("mean working hours:  $mean_labor")
        println("new capital stock: $kbar")
        println("new labor: $nbar")
        println("new tau: $taun")
        println("new pensions: $pen")
        println("new transfers: $trbar")


    end # outer while loop over aggregates, q = 1, .., nq

    return
end

@time main(glo,ef,sp1)
#vr1 = main(glo,ef,sp1)
#display(vr1[10,1])
