# AK60_value_procs.jl
#
# input into "AK60_value_main.jl"
# please make sure to save the two files in the same directory
#
# contains the functions used in "AK60_value_main.jl"
# for the computation of the solution to the Auerbach-Kotlikoff
# model in Heer/Maussner, DGE modeling, 2nd ed., 2009, Section 9.2.1
#
# Date: June 11, 2021
# Author: Burkhard Heer
#
# I would like to thank Chris Scharrer for helpful assistance.
# All errors are mine.
# For questions send an email to: Burkhard.Heer@wiwi.uni-augsburg.de
#


#******************
# utility function:
# eq. (9.2) in Heer/Maussner (2009)
#******************
function util(c,n)

	if (c>0.0) .& (eta!=1.0)
    	return 1/(1-eta).*((c*(1-n)^gam)^(1-eta)-1)
	elseif (c>0.0) .& (eta==1.0)
    	return log(c)+gam*log(1-n)
	end
end

#******************
# logarithmic grid:
#
# input:
# s0 - lower boundary
# s1 - upper boundary
# gp - number of grid points
#******************
function log_grid(s0,s1,gp)

	dist=s1-s0+1
	grid=exp.(range(log(1),stop=log(dist),length=gp))
	grid=grid.+s0.-1
	return grid
end



# wage function
function wage(capital,labor)
	y=(1-alp)*capital^alp*labor^(-alp)
	return y
end

# interest rate function
function interest(capital,labor)
	y=alp*capital^(alp-1)*labor^(1-alp)
	return y
end

# ****************************************
# value function iteration for the retired
# at age tw+tr-1,...tw+1
#
# inputs:
# func -- value function: valf_ret (retired), valf_wor(worker)
# itp -- interpolation nodes/values for value function next period
# 		at grid points k1_grid
# k_grid -- grid over capital stock at age "age"
# k1_grid -- grid over next-period capital stock at age "age+1"
# age -- age of the household
#
# ****************************************

function func_gs(func,itp,k_grid,k1_grid,age,glo)
	w = glo[1]
	r = glo[2]
	pen0 = glo[3]
	tau = glo[4]
	check_k0=0.0125 # check if at the boundaries kmin+check_k0*dist
						# (or kmax-check_k0*dist) implies lower rhs of
						# Bellman equation than the boundary point
	min_c0=0.0001	# constant: minimum consumption must be min_c0
						# For this value, we can compute the maximum possible
						# next period capital stock k1 which serves as
						# the boundary in the golden section search
	kdist=k1_grid[2]-k1_grid[1] # distance at lower grid point
	# value function and policy function
	# at age "age"  with capital k_grid[k0i]
	vfa = zeros(eltype(k1_grid[1,1]),length(k_grid),1)
	pfka = zeros(eltype(k1_grid[1,1]),length(k_grid),1)
	lower=minimum(k1_grid) # starting point for iteration over next-period wealth k1
	max_k1=maximum(k1_grid) # upper grid point of k1
	min_k1=minimum(k1_grid) # lowest grid point of k1
	glo1 = zeros(5)
	glo1[2:5] = glo
	for k0i in 1:1:length(k_grid)
		k0=k_grid[k0i,1]
		glo1[1]= k0

			# maximum crit point for k1: consumption at least min_c0
		if age>tw
			upper=k0*(1+r-del)+pen0-min_c0
		else
	       	upper=k0*(1+r-del)+(1-tau)*w*n_max-min_c0
		end

		# upper grid point on next-period wealth k1
		# 1) implied by c>c_min (upper)
		# 2) implied by maximum grid point with finite value
		#    of the value function
		upper=minimum([upper;max_k1])

		# check if the lower boundary is smaller than the
		# upper boundary for golden section search
		if lower>upper
			upper=minimum([lower+kdist;max_k1])
		end

		if lower<=upper
			# Golden Section Search searches minimum
			# we look for maximum => we multiply the function by -1 and
			# minimize it
			sc_valf=-1.00
			# check if boundary solution k1 = min_k1
			if lower==min_k1
				val0 = func(min_k1,glo1,itp)
				val1 = func(min_k1+check_k0*kdist,glo1,itp)
				if val0>val1  # boundary solution
					k1 = min_k1
					v0 = val0
				else # not a boundary solution even though the lower point is min_k1
					res=optimize(x->sc_valf*func(x,glo1,itp), lower, upper, GoldenSection(),abs_tol=abs_tol_gs,rel_tol=rel_tol_gs)
					# optimal next-period capital stock k1 and value of the value function
					# from Golden Section Search
					k1=Optim.minimizer(res)
					v0=Optim.minimum(res)*sc_valf
					if (isnan(Optim.minimum(res))==true) | (Optim.converged(res)==false) | isnan(k1)
						lower=min_k1
					else
						lower=k1
					end
				end
			elseif upper>=kmax
				val0 = func(kmax,glo1,itp)
				val1 = func(kmax-check_k0*kdist,glo1,itp)
				if val0>val1  # boundary solution
					k1 = kmax
					v0 = val0
				else # not a boundary solution even though the lower point is min_k1
					upper = kmax
					res=optimize(x->sc_valf*func(x,glo1,itp), lower, upper, GoldenSection(),abs_tol=abs_tol_gs,rel_tol=rel_tol_gs)
					# optimal next-period capital stock k1 and value of the value function
					# from Golden Section Search
					k1=Optim.minimizer(res)
					v0=Optim.minimum(res)*sc_valf
					if (isnan(Optim.minimum(res))==true) | (Optim.converged(res)==false) | isnan(k1)
						lower=min_k1
					else
						lower=k1
					end
				end

			else # not a boundary solution because the lower point is not min_k1
				res=optimize(x->sc_valf*func(x,glo1,itp), lower, upper, GoldenSection(),abs_tol=abs_tol_gs,rel_tol=rel_tol_gs)
				# optimal next-period capital stock k1 and value of the value function
				# from Golden Section Search
				k1=Optim.minimizer(res)
				v0=Optim.minimum(res)*sc_valf
				if (isnan(Optim.minimum(res))==true) | (Optim.converged(res)==false) | isnan(k1)
					lower=min_k1
				else
					lower=k1
				end
			end


			# exploiting the monotonocity of the value function
			# => the lower boundary for golden section search for k[ik+1]>k[ik]
			# is equal to the optimal value of k1 in the previous iteration


			vfa[k0i,1]=v0
			pfka[k0i,1]=k1
		else
			vfa[k0i,1]=valprob
			pfka[k0i,1]=valprob
		end

	end
	return vfa,pfka
end

#*************************
# value function retirees:
# input:
# k1 - next-period capital stock
# glo - global variables
# glo[1] - capital stock this period
# itp -- interpolation nodes value function next period
# output:
# right-hand side of the Bellman equation
#*************************
function valf_ret(k1,glo,itp)

	k0=glo[1]
	w=glo[2]
	r=glo[3]
	pen0=glo[4]
	tau = glo[5]

    c0=pen0+k0*(1+r-del)-k1
	if k1>0.0
    	v1=itp(k1)
	else # numerical inaccuracies may result in k1-5.4e-20
		v1 = itp(0.0)
	end

	if c0>0
		return util(c0,0)+beta*v1
	else
		return valprob*(1+c0^2)+beta*v1
	end

end


#************************
# value function workers:
#
# input:
# k1 - next-period capital stock
# glo - global variables
# glo[1] - capital stock this period
# itp -- interpolation nodes value function next period
# output:
# right-hand side of the Bellman equation
#
#************************

function valf_wor(k1,glo,itp)

	k0=glo[1]
	w=glo[2]
	r=glo[3]
	pen0=glo[4]
	tau = glo[5]

	# optimal labor supply
	# that follows from the first-order condition of the household
	# with respect to labor supply
    n0= 1/(1+gam)*(1-gam/((1-tau)*w)*((1+r-del)*k0-k1))
    if n0<0.0
        n0=0.0
    elseif n0>n_max
        n0=n_max
    end

	# consumption from budget constraint
	c0=k0*(1+r-del)+(1-tau)*w*n0-k1

	if (c0>0) & (n0>0)
		#display([k0 k1 n0 c0])
		if k1>0.0
    		v1=itp(k1)
		else
			v1=itp(0.0)
		end

		return util(c0,n0)+beta*v1
	else
		if k1>0
			v1=itp(k1)
		else
			v1=itp(0.0)
		end
		return valprob*(1+c0^2)+beta*v1
    end
end
