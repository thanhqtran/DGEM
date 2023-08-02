# AK60_value_main.jl
#
# computes the steady in the OLG model (Auerbach-Kotlikoff)
# with VALUE FUNCTION ITERATION
#
# in Heer/Maussner, "Dynamic General Equilibrium Modeling", 2009,
# second edition, Section 9.1.2
#
# Date: June 11, 2021
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
#Pkg.add("JLD")
using JLD
using Printf
using Plots
using QuantEcon
using Optim
using Roots
using NLsolve
using ForwardDiff
using Revise
using Interpolations
import Base.Threads.@spawn
include("AK60_value_procs.jl");


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
println("******* Section 9.2.1 (2nd edition) *******")
println("******* Auerbach-Kotlikoff Model **********")
println("*******************************************")
println("")


#************
# Parameters:
#************

# model:
const beta=0.96     # discount factor
const eta=2.00      # 1/IES (intertemporal elasticity of substitution)
const alp=0.36      # production elasticity of capital
const del=0.10    #  rate of depreciation
const tw=40         # years as a worker
const tr=20         # years as a retiree
const rep=0.30      # replacement rate of pensions
const gam=2.0       # weight: disutility from working
# numerical parameters
const valprob=-1e8  # large negative value to initialize rhs of Bellman equation
#const ftol_nl=1e-8      # tolerance final solution
#const xtol_nl=1e-12
const abs_tol_agg=1e-4 	# absolute tolerance for percentage deviation
						# of aggregate solution for K in two consecutive iterations
const abs_tol_gs=1e-8   # (absolute) tolerance golden section search
const rel_tol_gs=1e-12  # (relative) tolerance golden section search
const nq=100 		# maximum number of iterations over aggregate capital stock
const phi=0.80 		# updating parameter for aggregate variables

#
# parameters that are only used in the function main()
#

# DO NOT CHANGE iplot
ipolt=1   #1: linear  -- not implemented for cubic interpolation yet


# grid on individual wealth
kmin=0.0       # lower limit on individual capital grid
kmax=10.0       # upper limit of capital grid
gp_ngrid = 200        # number of grid points on assets
ln_grid=0       # 0 -- equispaced grid, 1 -- logarithmic grid
n_max = 0.60    # maximum labor supply

# initial guesses of global variables
nbar=0.2 		# initial value for aggregate labor
r0=0.045		# initial guess for interest rate
kbar=(alp/(r0+del))^(1/(1-alp))*nbar  # guess for aggregate capital

function main(gp_kgrid,kbar,nbar)

    # Step 1:
	# prepare grids

	# grid over individual capital stock
	kdist=(kmax-kmin)/(gp_kgrid)
	if ln_grid==0 # linear grid
		k_grid=kmin:kdist:kmax
	else	# logarithmic grid
		k_grid=log_grid(kmin,kmax,gp_kgrid)
	end

	#p=plot(k_grid,label="capital stock",xlabel="grid point number")
	#display(p)

	# grid = 1, 2, ..., tw+tr for all ages
	age_grid=1:1:(tw+tr)

	# array for the value function and policy functions
	# arrays of the dimension (tw+tr) x gp_ngrid

	# value function
	vfa=zeros(eltype(k_grid[1]),length(age_grid),length(k_grid))
	# optimal next-period capital stock
	pfka=zeros(eltype(k_grid[1]),length(age_grid),length(k_grid))
	# optimal labor supply
	pflab=zeros(eltype(k_grid[1]),tw,length(k_grid))

	# saving the aggregate variables in outer iteration q
	# kbar and nbar
	kq = zeros((nq,2))
	glo = zeros(4)
	q=0
	krit_agg = 1.0 + abs_tol_agg
	# Step 2:
	# computation of aggregate factor prices
	while (q<nq) & (krit_agg>abs_tol_agg)
		q=q+1
		w = wage(kbar,nbar)
		r = interest(kbar,nbar)
		tau =rep/(2+rep)     # income tax rate */
		pen0=rep*(1-tau)*w*nbar*(tw+tr)/tw

		# saving factor prices w,r and pensions pen as input
		# into the functions for the rhs of the Bellman eq.
		glo[1] = w
		glo[2] = r
		glo[3] = pen0
		glo[4] = tau
		kold=kbar
		nold=nbar
		kq[q,1]=kbar
		kq[q,2]=nbar

		#
		# Step 3: value function in the last period
		# of life at age tw+tr
		#
		# oldest household:
		age=tw+tr
		# create a vector of the same element type (e.g. Float64)
		# depending on the computer and of the same length,
		# all with zero entries
		grid_ret = zeros(eltype(k_grid[1]),length(k_grid))

		# value function with last period capital stock k_grid[ki]
		for ki in 1:1:length(k_grid)
			k=k_grid[ki,1]
			c=pen0+k*(1+r-del)
			if c>0
				grid_ret[ki]=util(c,0)
			else # large negative value for negative consumption
				grid_ret[ki]=valprob*(1+c^2)
			end
		end
		vfa[age,:]=grid_ret

		# value function iteration
		# age tw+tr-1,..,1

		for age in tw+tr-1:-1:1
			println("outer iteration q: $q, value function age: $age, K= $kbar")

			# to compute the value function at age "age", we need
			# to interpolate the value function at age "age+1"
			# we need to take care of the cases where we have
			# not stared a number in the respective entry
			# => in this case, the row [k,value(k)] is deleted
			int_grid=[k_grid vfa[age+1,:]]
			int_grid=int_grid[(isnan.(int_grid[:,2]).==false),:]
			k1_grid=int_grid[:,1]
			vfa_grid=int_grid[:,2]

			# depending on the interpolation (linear/cubic)
			# and the grid (equi-spaced/cubic) we need to use
			# different entries in the interpolation function
			# interpolate()
			# itp = interpolate(A, interpmode, gridstyle, λ, k)
			# Interpolate an array A in the mode determined by interpmode
			# and gridstyle with regularization following [1], of order k
  			# and constant λ


			if (ln_grid==0) & (ipolt==1) # linear interpolation with equispaced grid
				k1_grid=k1_grid[1,1]:kdist:k1_grid[length(k1_grid),1]
				itp=interpolate((k1_grid,),vfa_grid,Gridded(Linear()))
			elseif ln_grid==1	# linear interpolation with logarithmic grid
				itp=interpolate((k1_grid,),vfa_grid,Gridded(Linear()))
			end

			# choose value functions worker or retiree:
			if age>tw
				func = valf_ret
			else
				func = valf_wor
			end
			(vfa0,pfka0) = func_gs(func,itp,k_grid,k1_grid,age,glo)
			vfa[age,:] = vfa0
			pfka[age,:] = pfka0

			# computation of the optimal policy function for labor
			if age<tw+1
				# optimal labor supply initilization
				pflab0 = zeros(eltype(k1_grid[1,1]),length(k_grid),1)
				# optimal labor supply
				# that follows from the first-order condition of the household
				# with respect to labor supply

				for k0i in 1:1:length(k_grid)
					k0 = k_grid[k0i,1]
					k1 = pfka[age,k0i]
		    		n0= 1/(1+gam)*(1-gam/((1-tau)*w)*((1+r-del)*k0-k1))
		    		if n0<0.0
		        		n0=0.0
		    		elseif n0>n_max
		        		n0=n_max
		    		end
					pflab0[k0i,1] = n0
				end
				pflab[age,:] = pflab0
			end
		end

		#
		# computation of the aggregate capital stock and employment nbar
		#

		# capital stock over the life-time
		k_age=zeros(eltype(k_grid[1]),length(age_grid))
		c_age=zeros(eltype(k_grid[1]),length(age_grid))
		# labor supply id denoted by "l" in the book; since
		# this is easily confused with a "1" in the editing process
		# I suggest to use "n" or "labor" for labor instead
		labor_age=zeros(eltype(k_grid[1]),tw)

		# households starts with zeros capital (first grid point on k_grid)
		k_age[1] = 0.0
		k1 = pfka[1,1]
		labor = pflab[1,1]
		labor_age[1] = labor
		c_age[1] = (1-tau)*w*labor-k1 # budget constraint of the 1-year old worker
 		k_age[2] = k1
		k0 = k1

		# computation of the capital-age and labor-age profiles,
		# starting at age 2 (in first period, k=0)
		for age in 2:1:tw+tr
			# interpolation of policy functions at age "age" for capital stock "k0"
			if (ln_grid==0) & (ipolt==1) # linear interpolation with equispaced grid
				# interpolation optimal next-period capital stock
				itp_k=interpolate((k_grid,),pfka[age,:],Gridded(Linear()))
				# interpolation optimal labor supply at age "age"
				if age<tw+1
					itp_lab=interpolate((k_grid,),pflab[age,:],Gridded(Linear()))
				end
			elseif (ln_grid==0) & (ipolt==2) # cubic interpolation, equispaced grid
				itp2=interpolate(pfka[age,:],BSpline(Cubic(Line(OnGrid()))))
				itp_k=scale(itp2,k_grid)
				# interpolation optimal labor supply at age "age"
				if age<tw+1
					itp2=interpolate(pflab[age,:],BSpline(Cubic(Line(OnGrid()))))
					itp_lab=scale(itp2,k_grid)
				end

			elseif ln_grid==1	# linear interpolation with logarithmic grid
				itp_k=interpolate((k_grid,),pfka[age,:],Gridded(Linear()))
				if age<tw+1
					itp_lab=interpolate((k_grid,),pflab[age,:],Gridded(Linear()))
				end
			end

			k1 = itp_k(k0)
			if age<tw+tr
				k_age[age+1] = k1
			end
			if age<tw+1
				labor = itp_lab(k0)
				labor_age[age] = labor
				c0 = (1-tau)*w*labor+(1+r-del)*k0-k1
			else
				c0 = pen0+(1+r-del)*k0-k1
			end
			c_age[age] = c0
			k0 = k1
		end

		#
		# Aggregation: each cohort has equal measure 1/(tw+tr)
		#
		knew = mean(k_age)
		nnew = mean(labor_age)*tw/(tw+tr) # remember the measure of workers is tw/(tw+tr)

		# update of the aggregates
		kbar=phi*kold+(1-phi)*knew
		nbar=phi*nold+(1-phi)*nnew
		krit_agg = abs((kbar-kold)/kbar)



		if (q==nq) | (krit_agg<abs_tol_agg)
			println("")
			println("solution: K=$kbar , L=$nbar")
			println("number of outer iterations: q=$q")
			println("")
			p1=plot(age_grid,c_age,title="Consumption-Age Profile",lw=3,label="Consumption",)
			display(p1)

			p2=plot(age_grid,k_age,title="Capital-Age Profile",lw=3,label="Capital",)
			display(p2)

			age_grid1 = 1:1:tw
			p3 = plot(age_grid1,labor_age,title="Labor-Age Profile",lw=3,label="Labor",)
			display(p3)
		end	# if q==nq
	end  # while q<nq
	# In case you want to analyse some local variables from the function
	# main() (or any other) at the JULIA prompt
	# save("data.jld", "data", pfka)
	# -> analyse as follows: julia>load("data.jld")["data"]
	# -> julia>data
end

@time main(gp_ngrid,kbar,nbar)
