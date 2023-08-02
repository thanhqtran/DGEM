# AK70_stochastic_income_procs.jl
#
# input into "AK70_stochastic_income_main.jl"
# please make sure to save the two files in the same directory
#
# contains the functions used in "AK70_stochastic_income_main.jl"
# for the computation of the solution to the Auerbach-Kotlikoff
# model in Heer/Maussner, DGE modeling, 3rd ed., 2022(?), Section 10.1.1
#
# Date: June 12, 2021
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
function u(c,n)
	if (c>0.0) .& (eta!=1.0)
    	return 1/(1-eta).*((c^gam*(1-n)^(1-gam))^(1-eta))
	elseif (c>0.0) .& (eta==1.0)
    	return gam*log(c)+(1-gam)*log(1-n)
	end
end

# marginal utility from consumption
function uc(c,l)
	y=gam*c^(gam*(1-eta)-1)*(1-l)^((1-gam)*(1-eta))
	return y
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


# production function
function production(capital,labor)
	y=capital^alp*labor^(1-alp)
	return y
end

#****************************
# calculate gini coefficient:
#****************************
function get_gini(x,m)
	if m==0
		m=ones(Float64,size(x)[1],1)/size(x)[1]
	end

	if abs(sum(m)-1.0)>1e-12
		my_println("")
		my_println("*********************************************")
		my_println("total mass in get_gini() is not eqqal to one!")
		str="problem: sum(mass)="*string(sum(m))*"!"
		my_println(str)
		my_println("*********************************************")
	end

	dat=[collect(1:1:length(x)) m x]
	inc_s=dat[sortperm(dat[:,3],rev=false),:]
	big_h=zeros(Float64,size(inc_s)[1],1);
	small_h=inc_s[:,2];
	small_q=inc_s[:,3].*inc_s[:,2]./sum(inc_s[:,3].*inc_s[:,2]);
	big_q=zeros(Float64,size(inc_s)[1],1);
	big_h[1,1]=inc_s[1,2];
	big_q[1,1]=small_q[1,1];

	for i in 2:1:size(inc_s)[1]
		big_h[i,1]=big_h[i-1,1]+small_h[i,1];
		big_q[i,1]=big_q[i-1,1]+small_q[i,1];
	end

	dum=0.0;
	for i in 1:1:size(inc_s)[1]
		if i==1;
			dum=big_q[i,1]*small_h[i,1];
		else;
			dum=dum+(big_q[i-1,1]+big_q[i,1])*small_h[i,1];
		end
	end
	gini=1.0-dum;
	gini=(size(inc_s)[1]/(size(inc_s)[1]-1))*gini;

# mu=sum(x.*m)
# dum=0.0;
# for i in 1:1:size(inc_s)[1]
# for j in 1:1:size(inc_s)[1]
# dum=dum+m[i]*m[j]*abs(x[i]-x[j])
# end
# end
# gini2=dum/(2*mu)
# gini2=(size(inc_s)[1]/(size(inc_s)[1]-1))*gini2;
# display([gini gini2])
	return gini
end

function read_data()
    sp = XLSX.readdata("survival_probs.xlsx", "Tabelle1","A2:A77")
    ef = XLSX.readdata("survival_probs.xlsx", "Tabelle1","B2:B46")
    return sp,ef
end

function value1(a1,glo,itp0,sp1,age)
	a0 = glo[1]
	r = glo[2]
	pen = glo[3]
	trbar = glo[4]
    c=(1+(1-tauk)*(r-delta))*a0+pen+trbar-ygrowth*a1
    c=c/(1+tauc)
	# next period value function at x
	val1 = itp0(a1)
	y = u(c,0)+ygrowth^(gam*(1-eta))*sp1[nw+age]*beta*val1
    return y
end

#
# value2: value function of the worker
#
# valuef contains the next-period value function that needs to be interpolated
function value2(x,vglo,sp1,ye,py,ef,age,iperm,iy,itpo)
	a1 = x
	a0 = vglo[1]
	r = vglo[2]
	w = vglo[3]
	trbar = vglo[4]
	taun = vglo[5]
	taup = vglo[6]
    glo = vcat(a1,vglo)
    labor = optimal_labor(glo,ye,ef,age,iperm,iy)

	eff = ef[age]*perm[iperm]*exp(ye[iy]) # idiosyncratic efficiency
    c=(1+(1-tauk)*(r-delta))*a0+(1-taun-taup)*w*eff*labor+trbar-ygrowth*a1
    c=c/(1+tauc)



    if age==nw
        y=u(c,labor)+ygrowth^(gam*(1-eta))*sp1[age]*beta*itpo(a1)
    else
        y=u(c,labor)
		for iy1 in 1:ny
			y=y+py[iy,iy1]*ygrowth^(gam*(1-eta))*sp1[age]*beta*itpo(a1,iy1)
		end
    end
	return y
end



#
# computes optimal labor supply
#
function optimal_labor(yglo,ye,ef,age,iperm,iy)
	a1 = yglo[1]
    a0 = yglo[2]
	r = yglo[3]
	w = yglo[4]
	trbar = yglo[5]
	taun = yglo[6]
	taup = yglo[7]

    y= ef[age]
	y=perm[iperm]
	y=exp(ye[iy])
	eff = ef[age]*perm[iperm]*exp(ye[iy]) # idiosyncratic efficiency
    w0=(1-taun-taup)*w*eff
    labor = gam-((1+(1-tauk)*(r-delta))*a0+trbar-ygrowth*a1)*(1-gam)/w0
	# value of labor in admissable space?
    if labor<0
		labor = 0
	elseif labor>labormax
		labor = labormax
	end
    return labor
end



# marginal utility from consumption next period given a',iperm,j,i
function uc1(a1,iperm,iy,iage,cwopt,lopt,cropt)
    if iage<=nw
        if a1<=kmin
            c1 = cwopt[iperm,iy,1,iage]
            labor = lopt[iperm,iy,1,iage]
        elseif a1>=kmax
            c1 = cwopt[iperm,iy,na,iage]
            labor =  lopt[iperm,iy,na,iage]
        else
            j0=sum(a.<a1)+1
            j0=min(j0,na)
            n0=(a1-a[j0-1])/(a[j0]-a[j0-1])
			# linear interpolation
            c1 = (1-n0)*cwopt[iperm,iy,j0-1,iage] + n0*cwopt[iperm,iy,j0,iage]
            labor = (1-n0)*lopt[iperm,iy,j0-1,iage]+n0*lopt[iperm,iy,j0,iage]
        end
    else
        labor = 0
        if a1<=kmin
            c1 = cropt[1,iage-nw]
        elseif a1>=kmax
            c1 = cropt[na,iage-nw]
        else
            j0=sum(a.<a1)+1
            j0=min(j0,na)
            n0=(a1-a[j0-1])/(a[j0]-a[j0-1])
            c1=(1-n0)*cropt[j0-1,iage-nw]+n0*cropt[j0,iage-nw]
        end
    end
    return uc(c1,labor)
end
