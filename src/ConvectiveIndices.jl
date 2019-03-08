module ConvectiveIndices
export trapz
export saturation_vapor_pressure_liquid
export moistlaps
export virtualtemp
export liftparcel
export thermo_rh
export rh_to_r
export r_to_rh
export q_to_rh
export dewpoint_to_q
export calc_CAPE_thetae
export calc_BCL
export calc_entropy
export calc_entropy_inverse
export calc_dCAPE
export calc_dilute_CAPE

using Interpolations
using LinearAlgebra
using Statistics

"""
numerical integration of y (Vector) with respect to the coordinates specified by x using the trapezoidal rule
"""
function trapz(x::Vector{R}, y::Vector{R}) where R<:Real

	len = length(y)
    if (len != length(x))
        error("Vectors must be of same length")
    end
    r = R(0.1)
    @inbounds for i in 2:len
       r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
	r/R(2.0)
	
end

"""
Find the element in an array closest to a value, return the index and difference
"""
function mindist(x, val)
	# using enumerate to avoid indexing
	min_i = 0
	min_x = Inf
	for (i, xi) in enumerate(x)
		dist = abs(xi - val)
		if dist < min_x
			min_x = dist
			min_i = i
		end
	end
	return (min_i, min_x)
end

"""
saturation_vapor_pressure_liquid(T)
Calculate the saturation vapor pressure temperature [K] according to a formula presented in:
D. M. MURPHY and T. KOOP: Review of the vapour pressures of ice and supercooled water for atmospheric applications
Q. J. R. Meteorol. Soc. (2005), 131, pp. 1539-1565
"""
function saturation_vapor_pressure_liquid(T::Real)
	#if T < 123 || T > 332   
	#	@warn "Temperature out of range [123-332] K."
    #end
    temp = (54.842763f0 - (6763.22f0 / T) - 4.210f0  * log(T) + 0.000367f0  * T) + (tanh( 0.0415f0  * (T - 218.8f0 ) )*(53.878f0  - (1331.22f0 /T) - 9.44523f0 *log(T) + 0.014025f0  * T))
	
	es = exp(temp);
end

function saturation_vapor_pressure_liquid(T::Vector{R}) where R<:Real

    #if ~(all(123 .< T) && all(T .< 332)) 
	#	@warn "Temperature out of range [123-332] K."
	#end
	
	N = length(T)
	es = Vector{R}(undef,N)

	for i = 1:N
    	es[i] = exp((54.842763f0 - (6763.22f0 / T[i]) - 4.210f0 * log(T[i]) + 0.000367f0  * T[i]) + (tanh( 0.0415f0 * (T[i] - 218.8f0 ) )*(53.878f0 - (1331.22f0 /T[i]) - 9.44523f0 *log(T[i]) + 0.014025f0 * T[i])))
	end
end



"""
virtualtemp( tk, r )
Return virtual temperature given temperature in K and mixing ratio in g/g
"""
function virtualtemp(tk::Real, r::Real)
	gc_ratio =  287.04/461.5;#epsilon=Rd/Rv, where Rd=gas const. for dry air and Rv=gas const. for water vap.
	tk_v = tk.*(1.0+r/gc_ratio)./(1.0+r); #tk*(1 + 0.608*r);
end

"""
rh_to_r( tk, p, r )
Return mixing ratio [kg/kg] given temperature [K], pressure [hPa] and relative humidity [%] 
"""
function rh_to_r(tk::Real, p::Real, RH::Real)
	c = 18.0152/28.9644;
	es = saturation_vapor_pressure_liquid(tk);
	e = (RH/100)*es;
	Pd = p*100 - e;
	r = (e/Pd)*c
end

"""
r_to_rh( tk, p, r )
Return relative humidity [%] given temperature [K], pressure [hPa] and mixing ratio [kg/kg]
"""
function r_to_rh(tk::Real, p::Real, r::Real)
	es = saturation_vapor_pressure_liquid(tk);
	c = 18.0152/28.9644;
	e = (r*p*100)/(r + c);
	RH = 100 * e/es;
end

"""
q_to_rh( tk, p, r )
Return relative humidity [%] given temperature [K], pressure [hPa] and specific humidity 
"""
function q_to_rh(tk::Real, p::Real, q::Real)
	es = saturation_vapor_pressure_liquid(tk);
	c = 18.0152f0/28.9644f0;
	e = (q*p*100f0)/(c + (1 - c)*q);
	RH = 100 * e/saturation_vapor_pressure_liquid(tk);
end

function q_to_rh(tk::Vector{F},p::Vector{F},q::Vector{F}) where F<:AbstractFloat
	RH = @. 100 * ((q*p*100f0)/((18.0152f0/28.9644f0) + (1 - 18.0152f0/28.9644f0 )*q))  /saturation_vapor_pressure_liquid(tk);
end

"""
dewpoint_to_q( tk, p, td )
Return specific humidity [kg/kg] given temperature [K], pressure [hPa] and dewpoint temperature [K]
"""
function dewpoint_to_q(tk::Real, p::Real, td::Real)
	c = 18.0152/28.9644;
	e = saturation_vapor_pressure_liquid(td);
	Pd = p*100 - e;
	q = (e*c)/(e*c+Pd);
end


"""
function (theta, thetae, thetaes, p_lcl, tk_lcl, r) = thermo_rh(tk,p,rh)
thermo_rh generates thermodynamic variables from t(K) p(Pa) rh(%) 
output [theta, thetae, r, tk_lcl, p_lcl] in K, K, kg/kg, K, hPa; r = water vapour mixing ratio
"""
function thermo_rh(tk::Real,p::Real,rh::Real)

	# changed rsat denominator from p to p-es (3 may 07)
	# Realised thetae calculation is incorrect (17 Nov 2009)
	# Changed to used formulae of Bolton (1980), for pseudo-equivalent p.t.

	#convert t,p to SI units
	#p = p.*100f0; #Pa

	p0 = 100000f0;	#reference pressure in Pa
	R = 287f0;		#gas constant
	K = R/1004f0;
	epsi = 0.62198f0;

	a = ((p/p0)^(-K));
	theta = tk*a;

	es = saturation_vapor_pressure_liquid(tk); # Saturation vapor pressure
	e = ((rh/100)*es);                              # vapour pressure (Pa)

	# Calculate water vapour mixing ratio r and q specific humidity
	r = (epsi*e)/(p-e);
	rsat = (epsi*es)/(p-es);
	#ri = r - rsat; #liquid water content of air (Zhang et al 1990)
	k = 0.2854f0 * (1 - 0.28f0*r);
	# change units from g/g to g/kg
	rg = r*1f3; rsat = rsat*1f3;

	# calculate pseudo-equivalent potential temperature, from Bolton, Mon Wea Rev, 1980
	# r = is g/kg
	# Firstly calculate Temp at LCL, note e must be in mb.
	tk_lcl = ( 2840 / (3.5f0 * log(tk) - log(e/100) - 4.805f0) ) + 55;               # eqn 21, Bolton
	thetae = theta*exp(( (3.376f0/tk_lcl) - 0.00254f0)*rg*(1+0.81f0 * rg *1f-3));   # eqn 38, Bolton
	thetaes = theta*exp(( (3.376f0/tk_lcl) - 0.00254f0)*rsat*(1+0.81f0 * rg *1f-3));   # eqn 38, Bolton

	#LCL height using Poisson's equation
	p_lcl =  0.01f0*p*((tk_lcl/tk)^(1.0f0/k));

	return theta, thetae, thetaes, p_lcl, tk_lcl, r 

end

"""
function (theta, thetae, thetaes, p_lcl, tk_lcl, r) = thermo_rh(tks,ps,rhs)
thermo_rh generates thermodynamic variables from t(K) p(Pa) rh(%) N-element arrays 
output [theta, thetae, r, tk_lcl, p_lcl] in K, K, kg/kg, K, hPa; r = water vapour mixing ratio
"""
function thermo_rh(tk::Vector{F},p::Vector{F},rh::Vector{F}) where F<:AbstractFloat

	# changed rsat denominator from p to p-es (3 may 07)
	# Realised thetae calculation is incorrect (17 Nov 2009)
	# Changed to used formulae of Bolton (1980), for pseudo-equivalent p.t.
	#p = p.*100f0; #Pa
	p0 = F(100000);	#reference pressure in Pa
	R = F(287);		#gas constant
	K = R/F(1004);
	epsi = F(0.62198);

    len = length(tk)
    if (len != length(p))
		error("All vectors must be of same length")
    end
		
    theta = zero(tk)
    thetae = zero(tk)
	thetaes = zero(tk)
	p_lcl = zero(tk)
	tk_lcl = zero(tk)
	#r = zero(tk)
	#rsat = zero(tk)

	@inbounds for i in 1:len
		#p[i] = 100*p[i]  # convert t,p to SI units
		# a = ((p[i]/p0)^(-K))
		theta[i] = tk[i] * ((p[i]/p0)^(-K)) #a
		es = saturation_vapor_pressure_liquid(tk[i]) # Saturation vapor pressure
		e = (rh[i]/ 100 )*es  # vapour pressure (Pa)
		
		# Calculate water vapour mixing ratio r and q specific humidity
		r = (epsi*e)/(p[i]-e)
		rsat = (epsi*es)/(p[i]-es)

		#ri = r[i] - rsat[i]; #liquid water content of air 	

		# calculate pseudo-equivalent potential temperature, from Bolton, Mon Wea Rev, 1980
		# r = is g/kg
		# Firstly calculate Temp at LCL, note e must be in mb.
		tk_lcl[i] = ( F(2840) / (F(3.5) * log(tk[i]) - log(e/100) - F(4.805)) ) + 55;               # eqn 21, Bolton
		thetae[i] = theta[i]  * exp( ((F(3.376)/tk_lcl[i]) - F(0.00254)) * r*1000 * (1+F(0.81) *r*1000 *F(0.001)));   # eqn 38, Bolton
		thetaes[i] = theta[i] * exp( ((F(3.376)/tk_lcl[i]) - F(0.00254)) * rsat*1000 * (1+F(0.81) *r*1000 *F(0.001)));   # eqn 38, Bolton

		#LCL height using Poisson's equation
		k = F(0.2854) * (1 - F(0.28)*r);
		p_lcl[i] =  F(0.01)*p[i]*((tk_lcl[i]/tk[i])^(1/ k ));

    end

	return theta, thetae, thetaes, p_lcl, tk_lcl#, r, rsat

end


"""
	calc_CAPE_thetae (ps[hPa], tks[K], qs[kg/kg], zs[m], parcel, dp_mix[hPa], dp_intp[hPa], kiss) 

Calculate convective indices such as **CAPE** and **CIN** for a parcel of choice (surface/most unstable, with/without vertical mixing), using a theta-e formulation.

# Examples

```jldoctest 
julia> LI, CAPE, CIN = calc_CAPE_theta(ps,tks,qs,zs) # most unstable parcel, mixed over 50 hPa (default)
(-8.94582333506736, 1613.7159227760612, 327.257167221434))
julia> LI, CAPE, CIN = calc_CAPE_theta(ps,tks,qs,zs, parcel = 2, dp_mix = 0) # surface parcel, not mixed
(-12.416924139871522, 2428.182537242374, 85.85516940477973)
julia> LI, CAPE, CIN, pLCL, zBCL, CAPECIN_ALCL, CIN_LCL, MRH_ALCL, MRH1, MRH2 = calc_CAPE_thetae(ps,tks,qs,zs, FULL = 1) # full calculations
(-8.94582333506736, 1613.7159227760612, 327.257167221434, 936.6429885118564, 1230.0, -189.68905798724995, 128.5705360872618, 69.90722164805184, 56.290565968008316, 30.494525283693054)
```

**OUTPUT** for FULL=1: \n
Lifted Index [°C], CAPE [J/kg], CIN [J/kg], pLCL [hPa], Buoyant Condensation Level [m], CAPE-CIN above the LCL [J/kg], CIN below LCL [J/kg], MRH (mean RH%) above the LCL [%], MRH 600-800 hPa, MRH 300-600 hPa \n

**INPUT**:
(N-element ARRAYs) **ps**,**tks**,**qs**,**zs** = vertical profiles of pressure, temperature, specific humidity and geopotential height \n
OPTIONAL keyword arguments: \n
`parcel = 1` : the most unstable parcel in the lowest 350 hPa (default) \n
`parcel = 2` : surface parcel, or parcel from the lowest level \n
`dp_mix = 0...100` : pressure layer depth [hPa] for mixing the source parcel (default 50, use 0 for no mixing) \n
`dp_intp = 5` linearly interpolate to a uniform pressure grid with resolution dp (default 5). Use 0 to skip. This is a lot faster, but disables mixing and FULL option,
and is not recommended for low-resolution input. \n
`FULL = 1`: Full calculations to include also less known convective predictors, which were found useful in [1]. \n

This routine uses an `equivalent potential temperature formulation` for all indices (similarly to ECMWF CAPE), avoiding vertical loops altogether. 
This results in larger absolute values (e.g. for CAPE, 30% larger) than classic computations, but in no worse correlation with observed convection [1].

NOTE: Default option `parcel=1` and `dp_mix=50` corresponds to a hybrid mixed-layer most-unstable parcel similar to the one used by ECMWF. The MLMU-Lifted Index was the overall
thunderstorm index for Europe in [1].

[1] Ukkonen and Mäkelä (2019): Evaluation of machine learning classifiers for predicting deep convection
"""
function calc_CAPE_thetae(ps::Vector{F},tks::Vector{F},qs::Vector{F},zs::Vector{F}; parcel::Integer=1,dp_mix::Real=50,dp_intp::Real=5,FULL::Integer=0) where F<:AbstractFloat

	g = F(9.80665)

	#Check that profiles are descending, otherwise flip the arrays
	if ps[end] < ps[1]
		tks = reverse(tks,1); rhs = reverse(rhs,1); ps = reverse(ps,1)
	end
  
	sp = ps[end]
	
	rhs = q_to_rh(tks,ps,qs)
	theta,thetae,thetaes = thermo_rh(tks,100*ps,rhs)

	if dp_intp == 0 # No interpolation to a uniform pressure grid 
		FULL = 0   # Full calculations are skipped 
		dp_mix = 0 # No linear parcel mixing is done
		thetae_env = thetae 
		thetaes_env = thetaes
		tk_env = tks; rh_env = rhs; z_env = zs
		pres = ps
		nlevs = length(pres)
	else # Interpolation
		pres = collect(ceil(ps[1]):dp_intp:floor(sp))
		nlevs = length(pres)
	
		# use "Interpolations" package for linear interpolation on a non-uniform grid
		knots = (ps,);
	
		itp = interpolate(knots,tks,Gridded(Linear())) 
		tk_env = itp(pres)
		itp = interpolate(knots,rhs,Gridded(Linear()))
		rh_env = itp(pres)
		itp = interpolate(knots,zs,Gridded(Linear()))
		z_env = itp(pres)
	
		itp = interpolate(knots,thetae,Gridded(Linear())); thetae_env = itp(pres)
		itp = interpolate(knots,thetaes,Gridded(Linear())); thetaes_env = itp(pres)
	end

	# FIND LAUNCH PARCEL
	if parcel == 1 #	--> MOST UNSTABLE
		# Find level with the highest theta-e in the lowest 350 hPa
		thetae_env[pres.<(sp - 350)] .= 0; # Set elements above the lowest 350 hPa to 0
		thetae_max,iPL = findmax(thetae_env); #Index of parcel level
	else #				--> SURFACE PARCEL
		iPL = nlevs # 
	end

	# PARCEL MIXING WITH ENVIRONMENT
	# Mix the properties over depth specified by dp_mix (e.g. 50 hPa)
	@views if dp_mix > 0
		#how many indices correspond to the mixing depth divided by two? (50 hPa mixing depth = 25 below, 25 above)
		dind = floor(Int,(dp_mix/2)/dp_intp)
		if nlevs > (iPL+dind)
			p0 = mean(pres[iPL-dind:iPL+dind])
			tk0 = mean(tk_env[iPL-dind:iPL+dind])
			rh0 = mean(rh_env[iPL-dind:iPL+dind])
		else
			p0 = mean(pres[iPL-2*dind:end]) 
			tk0 = mean(tk_env[iPL-2*dind:end])
			rh0 = mean(rh_env[iPL-2*dind:end])
		end
	else # NO MIXING
		p0 = pres[iPL]
		tk0 = tk_env[iPL]
		rh0 = rh_env[iPL]
	end
	
	# CALCULATE THETA-E AND LCL
	_,thetae_parc,_,pLCL = thermo_rh(tk0,100*p0,rh0)

	# LIFTED INDEX
	i500 = mindist(pres,500)[1]
	LI = thetaes_env[i500] - thetae_parc # T (around 500 hPa ) - T (parcel lifted to 500 hPa)

	# CAPE AND CIN
	#the quantity being integrated in theta-e formulation for CAPE
	tdiff = (thetae_parc .- thetaes_env)./thetaes_env

	iLCL = findlast(pres.<pLCL)
	iLFC = findlast(tdiff[1:iLCL].>0)
	iEL = findfirst(tdiff.>0)

	if iLFC == nothing
		CAPE = F(0)
		CIN = NaN
		iLFC = NaN
	else
		if (iLFC == iEL)
			iEL=iEL-1;
		end
		#Only positive regions in sounding contribute to CAPE?
		tdiff_cape = copy(tdiff)
		tdiff_cape[tdiff_cape.<0] .= 0
		CAPE = -g*trapz(z_env[iEL:iLFC],tdiff_cape[iEL:iLFC])
		CIN = g*trapz(z_env[iLFC:iPL],tdiff[iLFC:iPL])
	end
	
	# RETURN EARLY?
	if FULL == 0 # Stop here and return with only LI,CAPE,CIN 
		return LI,CAPE,CIN
	end
	
	# CAPE and RH in the 250-hPa depth layer ABOVE LCL	
	idp_250 = div(250,dp_intp)
	if iLCL-idp_250 > 0
		MRH_ALCL = @views mean(rh_env[iLCL-idp_250:iLCL])
		CAPECIN_ALCL = -g*trapz(z_env[iLCL-idp_250:iLCL],tdiff[iLCL-idp_250:iLCL])
	else
		MRH_ALCL = @views mean(rh_env[1:iLCL])
		CAPECIN_ALCL = -g*trapz(z_env[1:iLCL],tdiff[1:iLCL])
	end

	if iLCL < iPL
		CIN_LCL = g*trapz(z_env[iLCL:iPL],tdiff[iLCL:iPL])
	else
		CIN_LCL = F(0)
	end

	i300 = mindist(pres,300)[1]
	i600 = mindist(pres,600)[1]
	i800 = mindist(pres,800)[1]

	# MEAN RELATIVE HUMIDITIES 
	MRH1 =  @views mean(rh_env[i600:i800])
	MRH2 =  @views mean(rh_env[i300:i600])
	
	# INHIBITON
	# Calculate the "Buoyant Condensation Level" (BCL) which is "the level at which saturation would occur
	# through buoyant mixing alone due to sensible heating from the surface"
	# This is an inhibition measure by Tawfik et al. (2014):
	#  A process-based framework for quantifying the atmospheric preconditioning of surface-triggered convection
	
	zBCL = calc_BCL(qs,rhs,zs)
	
	return LI, CAPE, CIN, pLCL, zBCL, CAPECIN_ALCL, CIN_LCL, MRH_ALCL, MRH1, MRH2

end


"""
```jldoctest 
julia> CAPE_dilute = calc_dilute_CAPE(ps,tks,qs,zs)
```
where the input variables are `descending columns` of: pressure [hPa], temperature [K], specific humidity [kg/kg] and height [m],
returns the `diluted CAPE` using a constant mass entrainment rate. This is a simplified version of subroutine buoyan_dilute (by Richard Neale) in the CAM model.

The parcel is launched from the level with the highest theta-e in the lowest 350 hPa of the atmosphere. Optionally, the index can be specified by the user: 
calc_dilute_CAPE(ps,tks,qs,zs,parcel_index=30), would use a surface parcel if length(ps) = 30.

"""

function calc_dilute_CAPE(ps::Vector{F},tks::Vector{F},qs::Vector{F},zs::Vector{F}; parcel_index::Integer=0,freezing::Integer=0) where F<:AbstractFloat
	# p = pressure [hPa], tk = temp [K], q = spec. hum [kg/kg], z = height [m], s = entropy (J/kg)
	sp = ps[end]
	rhs = q_to_rh(tks,ps,qs)

	if parcel_index > 0 # The parcel already decided
		ilaunch = parcel_index
	else
		theta,thetae = thermo_rh(tks,100*ps,rhs)
		thetae[ps.<(sp - 350)] .= 0; # Set elements above the lowest 350 hPa to 0
		thetae_max,ilaunch = findmax(thetae); #Index of parcel level
	end
	println("ilaunch: ",ilaunch)
	# Initialize values
	N = length(ps)

	tfreez  = 	F(273.15) 	#K
	g = F(9.80665)			# Gravity
	R = F(287.04)			# gas constant
	tk_mix = zero(tks)		# Tempertaure of the entraining parcel
	qt_mix = zero(tks)		# Total water of the entraining parcel
	qsat_mix = zero(tks)	# Saturation mixing ratio of the entraining parcel
	s_mix = zero(tks)		# Entropy of the entraining parcel
	tdiff = zero(tks) 	 	# (virtual) temperature difference between environment and parcel
	qtp = F(0)  			# Parcel total water	
	sp = F(0)				# Parcel entropy
	mp = F(0)				# Parcel relative mass flux

	# At parcel launch level
	qtp0 = qs[ilaunch]		# Parcel launch total water
	mp0 = F(1)				# Parcel launch relative mass flux	
	tk_guess = tks[ilaunch] # Parcel launch temperature 
	sp0 = calc_entropy(tk_guess,ps[ilaunch],qtp0)	# Parcel launch entropy
	s_mix[ilaunch] = sp0
	qt_mix[ilaunch] = qtp0
	tk_mix[ilaunch],qsat_mix[ilaunch] =  calc_entropy_inverse(sp0,ps[ilaunch],qtp0,tk_guess)


	for i in range(ilaunch-1,1,step=-1)
		# ps[i] is the current level, loop traverses back the array = upwards in the atmos. column
		dp = ps[i]-ps[i+1] # dp is negative

		p_env = F(0.5) * (ps[i]+ps[i+1])			# Environmental pressure
		tk_env = F(0.5) * (tks[i]+tks[i+1])			# Environmental temperature	
		qt_env = F(0.5) * (qs[i]+qs[i+1])			# Environmental total water
		s_env = calc_entropy(tk_env,p_env,qt_env)	# Environmental entropy

		dmpdz	= -F(1e-3) 				# Entrainment rate (/m)
		dpdz 	= -(p_env*g)/(R*tk_env) # hPa / m
		dzdp 	= 1/dpdz 				# m / hPa      
		dmpdp 	= dmpdz*dzdp 			# Entrainment rate (/hPa)
		dmpdp = 0
		
		# Sum entrainment to current level
		# entrains q,s out of intervening dp layers, in which linear variation is assumed
		# so really it entrains the mean of the 2 stored values.
         sp  = sp  - dmpdp*dp*s_env 
         qtp = qtp - dmpdp*dp*qt_env 
		 mp  = mp  - dmpdp*dp
		 
		 #println(" p_env			sp				qtp				mp")
		 #println("$p_env	$sp		$qtp	$mp")

		 # Entrain s and qt to next level.
         s_mix[i]  = (sp0  +  sp) / (mp0 + mp)
		 qt_mix[i] = (qtp0 + qtp) / (mp0 + mp)
	
		# The new parcel entropy and total water due to mixing are now calculated

		# Invert entropy from s and q to determine T and qsat of mixture 
		# t[i,k] used as a first guess so that it converges faster.

         tk_guess = tk_mix[i+1]
		 tk_mix[i],qsat_mix[i] =  calc_entropy_inverse(s_mix[i],ps[i],qt_mix[i],tk_guess) 
		 
		 #tk_mix_v = tk_mix[i]
		 #tk_env_v = tks[i]
		 tk_mix_v = virtualtemp(tk_mix[i],(qt_mix[i]/(1-qt_mix[i])))
		 tk_env_v = virtualtemp(tks[i],(qs[i]/(1-qs[i])))
		 tdiff[i] = (tk_mix_v + 0.5 - tk_env_v)/tk_env_v
		 #println("p, tmix, qmix")
		 #println([p_env,tk_mix[i],qsat_mix[i]])
	end
	#return tk_mix,qsat_mix,qt_mix
	if freezing==1
		xsh2o = F(0)
		ds_xsh2o = F(0)		# Entropy change due to loss of condensate
		ds_freeze = F(0)	# Entropy change dut to freezing of precipitation

		# ! Set parcel values at launch level assume no liquid water. 
		tp[ilaunch] = tk_mix[ilaunch]
		qstp[ilaunch] = qtp0
		tkv[ilaunch] = virtualtemp(tp[ilaunch],qstp[ilaunch])

		for i in range(ilaunch-1,1,step=-1)
			# Iterature 2 times for s,qt changes
			for j in [1,2]
				# Rain (xsh2o) is excess condensate, bar LWMAX (Accumulated loss from qtmix).
				xsh20[i] = maximum(F(0),(qt_mix[i] - qsat_mix[i] -lwmax) )

				# contribution to entropy from precip loss of condensate  (Accumulated change from smix).(-ve)
				ds_xsh2o[i] = ds_xsh2o[i] - cpliq * log (tk_mix[i]/tfreez) * maximum(F(0),(xsh2o[i] - xsh2o[i+1]) )
				
				# Entropy of freezing: latice times amount of water involved divided by T
				 
				if ( (tk_mix[i] <= tfreez+tscool) & (ds_freeze[i+1] == F(0)) ) # One off freezing of condensate. 
					ds_freeze[k]) = (latice/tk_mix[i]) * maximum(F(0),qt_mix[i]-qsat_mix[i]-xsh2o[i]) ! Gain of LH
				end if
				 
				if ( (tk_mix[i] <= tfreez+tscool) & (ds_freeze[i+1] != F(0)) ) #  Continual freezing of additional condensate.
					ds_freeze[k] = ds_freeze[i+1]+(latice/tk_mix[i]) * maximum(F(0),(qsat_mix[i+1]-qsat_mix[i]))
				end if

				# Adjust entropy and accordingly to sum of ds (be careful of signs).
				new_s = s_mix[i] + ds_xsh2o[i] + ds_freeze[i] 
				# Adjust liquid water
				new_q = qt_mix[i] - xsh2o[i]

				# Invert entropy to get updated tk_mix and qs_mix of parcel
				tk_guess = tk_mix[i]

				tk_mix[i],qsat_mix[i] =  calc_entropy_inverse(new_s,ps[i],new_q,tk_guess) 

			end
		end

	end

	iLCL = findlast(qt_mix .>= qsat_mix)
	iLFC = findlast(tdiff[1:iLCL].>0)
	iEL = findfirst(tdiff.<0)
	#return iLCL,iLFC,iEL
	if iLFC == nothing
		CAPE = F(0)
		#CIN = NaN
		iLFC = NaN
	else
		if (iLFC == iEL)
			iEL=iEL-1;
		end
		if iEL == nothing
			iEL = 1
		end
		#Only positive regions in sounding contribute to CAPE??
		tdiff_cape = copy(tdiff)
		tdiff_cape[tdiff_cape.<0] .= 0
		CAPE = -g*trapz(zs[iEL:iLFC],tdiff_cape[iEL:iLFC])
		#CIN = g*trapz(zs[iLFC:ilaunch],tdiff[iLFC:ilaunch])
	end
	#return CAPE
	return iLCL,iLFC,iEL,CAPE #,CIN,tdiff
end

function calc_dCAPE(ps::Vector{F},tks::Vector{F},qs::Vector{F},zs::Vector{F},adv_q::Vector{F},adv_tk::Vector{F}) where F<:AbstractFloat

	g = F(9.80665)
	#CAPE calculated "normally", in this case diluted CAPE:
	CAPE = calc_dilute_CAPE(ps,tks,qs,zs)
	#CAPE where environmental (but only environmental) tk,q gain an increment from advection
	#T = T0 + adv(T)*dt = T0 + u*(dT/dx)*dt + v*(dT/dy)*dt
	# adv_T is (u*dT/dx) = [ m/s * K/m] = [K/s]. Multiply by 3600 to get hourly increment
	theta,thetae = thermo_rh(tk_env,100*ps,rhs)

	# itp = interpolate(knots,thetae,Gridded(Linear())); 
	# thetae_env = itp(pres)
	# itp = interpolate(knots,thetaes,Gridded(Linear())); 
	# thetaes_env = itp(pres)

	thetae_env = thetae; 

	thetae_env[ps.<(sp - 350)] .= 0; # Set elements above the lowest 350 hPa to 0
	tmp,iparc = findmax(thetae_env);

	tk0 = tks[iparc]; q0 = qs[iparc];
 	tks = tks .+ adv_tk*3600
	qs = qs .+ adv_q*3600
	tks[ind] = tk0; qs[ind] = q0

	CAPE_adv = calc_dilute_CAPE(ps,tks,qs,zs,iparc)

	#dCAPE = CAPE - CAPE_new

end

function calc_entropy(tk::F,p::F,qtot::F) where F<:AbstractFloat
	# %T(K), p(hPa), qtot (kg/kg) 

	tfreez  = 	F(273.15) #K
	pref 	= 	1000 # hPa
	rl 		=	F(2.501e6) # J/kg latent heat of vaporization at 0C
	cpliq	=	F(4.188e3) # J/kg/K specific heat of liquid water
	cpwv	=	F(1.810e3) # J/kg/K specific heat of atmosphere water (wv)
	cpres	=	F(1.00464e3) # J/kgK
	c1		=	F(6.112)
	c2		=	F(17.67)
	c3		=	F(243.5)
	epsl	=	F(0.62197) #ratio of gas constants of air and wv
	rgas	=	F(287.04) #J/kg/K gas constant for dry air
	rh2o	=	F(461.5) #J/kg/K gas constant for water vapor
	#eref=6.106 #sat p at tfreez (mb)
	#temperature dependent latent heat of vaporization
	L= rl - (cpliq-cpwv)*(tk-tfreez); #T will be in celsius for this calculation
	
	esat = F(0.01)*saturation_vapor_pressure_liquid(tk)
	#esat= c1*exp(c2* (tk-tfreez)./(c3+tk-tfreez)) # This formula is not accurate when T < 0 degC

	qsat = epsl*esat/(p-esat) ; #saturation mixing ratio (in kg/kg)
	#sat vapor pressure equation from weather.gov used in zm_conv:
	#https://www.weather.gov/media/epz/wxcalc/rhTdFromWetBulb.pdf

	qv=min(qtot,qsat);  #partition qtot into vapr part only..not sure on why they do this?
	e=qv*p / (epsl+qv); #partial pressure exerted by wv ..is this ignoring -q*epsilon?
	#https://svn.ssec.wisc.edu/repos/bennartz_group/LIBRARY/idl/std_libs/colib/humidity.pro

	s_entropy= (cpres + qtot*cpliq)*log(tk/tfreez) - rgas*log((p-e)/pref) + L*qv/tk - (qv*rh2o)*log(qv/qsat);
	return s_entropy
	#with qt=qsat

	end

function calc_entropy_inverse(s::F,p::F,qtot::F,tk_guess::F) where F<:AbstractFloat
	# Given the entropy, pressure and total water, iteratively solve the entropy eq. for temperature using Newtons method
	# Also returns the saturated vapor mixing ratio qsat
	tfreez  = 	F(273.15) #K
	pref 	= 	1000 # hPa
	rl 		=	F(2.501e6) # J/kg latent heat of vaporization at 0C
	cpliq	=	F(4.188e3) # J/kg/K specific heat of liquid water
	cpwv	=	F(1.810e3) # J/kg/K specific heat of atmosphere water (wv)
	cpres	=	F(1.00464e3) # J/kgK
	c1		=	F(6.112)
	c2		=	F(17.67)
	c3		=	F(243.5)
	epsl	=	F(0.62197) #ratio of gas constants of air and wv
	rgas	=	F(287.04) #J/kg/K gas constant for dry air
	rh2o	=	F(461.5) #J/kg/K gas constant for water vapor

	loopmax = 100
	tk = tk_guess
	i = 0
	while i < loopmax
		L = latent_heat_vaporization(tk)
		esat = F(0.01)*saturation_vapor_pressure_liquid(tk)
		qsat = epsl*esat/(p-esat)
		qv = min(qtot,qsat)
		e = qv*p / (epsl+qv)
		
		#ds1 = calc_entropy(tk,p,qtot)
		ds1 = (cpres + qtot*cpliq)*log(tk/tfreez) - rgas*log((p-e)/pref) + L*qv/tk - (qv*rh2o)*log(qv/qsat) - s
		L = latent_heat_vaporization(tk-1)
		esat = F(0.01)*saturation_vapor_pressure_liquid(tk-1)
		qsat = epsl*esat/(p-esat)
		qv = min(qtot,qsat)
		e = qv*p / (epsl+qv)

		ds2 = (cpres + qtot*cpliq)*log((tk-1)/tfreez) - rgas*log((p-e)/pref) + L*qv/(tk-1) - (qv*rh2o)*log(qv/qsat) - s
		dTs = ds1/(ds2 -ds1)
		tk = tk + dTs

		if abs(dTs) < F(0.001)
			break
		end
		i += 1
		
	end
	esat = F(0.01)*saturation_vapor_pressure_liquid(tk)
	qsat = epsl*esat/(p-esat)
	return tk,qsat

end

function latent_heat_vaporization(tk::F) where F<:AbstractFloat
	return F(2.501e6) - (F(4.188e3)-F(1.810e3))*(tk-F(273.15))
end

"""
calc_BCL
Calculate the 'Buoyant Condensation Level' (**BCL**) which is 'the level at which saturation would occur
through buoyant mixing alone due to sensible heating from the surface'

This is an inhibition measure for quantifying how preconditioned the atmosphere is to moist convection. \n
Described in Tawfik et al. (2014):
A process-based framework for quantifying the atmospheric preconditioning of surface-triggered convection
"""
function calc_BCL(qs::Vector{F},rhs::Vector{F},zs::Vector{F}) where F<:AbstractFloat

	#qsats = F(100.0) ./rhs .* qs
	
	N = length(qs)
	qmix = Vector{F}(undef,N)
	qsats = Vector{F}(undef,N)

	@views @inbounds for i = 1:N
		qsats[i] = (100/rhs[i]) * qs[i]
		qmix[i] = mean(qs[i:N])
	end
	
	z = vcat( collect(10000:-300:6100), collect(6000:-100:3000), collect(3000:-30:zs[end]) )

	knots = (reverse(zs,1),);

	itp = interpolate(knots,reverse(qmix,1),Gridded(Linear())) 
	qmix_int = itp(z)
	itp = interpolate(knots,reverse(qsats,1),Gridded(Linear())) 
	qsat_int = itp(z)
	iBCL = findlast(qmix_int.>qsat_int)
	
	zBCL = z[iBCL]
	return zBCL

end

"""
liftparcel( tk0, p0, rh0, dp)
Calculate the vertical profile of a lifted parcel based on pseudo-adiabatic ascent
INPUT: tk0,p0,rh0 = parcel temperature [K], pressure [hPa], RH[%]
dp = vertical resolution for simulated parcel ascent (default 5 hPa, decreasing this value increases computational time)
OUTPUT: ts = temperature profile (from 100 hPa to p0 with dp increments)
tvs = virtual temp. profile [K],
ps = pressure profile [hPa],
rs = mixing ratio profile [kg/kg],
plcl = LCL pressure [hPa],
tlcl = LCL temperature [K]
"""
function liftparcel(tk0::Real, p0::Real, rh0::Real,dp::Real=5)

	eps = 0.62198; 	#Rd/Rv

	#Calculate LCL height in hPa p_lcl
	_,_,_,p_lcl,tk_lcl,r0 = thermo_rh(tk0,100*p0,rh0);

	pres = 100:dp:p0;

	n = length(pres) - 1;
	r_parc = zeros(n); #mixing ratio
	tk_parc = zeros(n);  
	tkv_parc = zeros(n); 

	i_lcl = findlast(pres.<p_lcl);

	#From p0 to LCL, use dry adiabatic lapse rate:
	r_parc[i_lcl:end] .= r0;
	tk_parc[i_lcl:end] = tk0.*(pres[i_lcl:end-1]./p0).^0.286;
	tkv_parc[i_lcl:end] = virtualtemp.(tk_parc[i_lcl:end],r_parc[i_lcl:end]);

	i = i_lcl - 1;
	# LCL reached, change to moist adiabatic lapse rate

	if (i_lcl == n+1)
	    tk_parc[end] = tk0
	    tkv_parc[end] = virtualtemp(tk0,r0)
	    i = i - 1;
	end

	while i > 0
		lapserate = moistlaps((pres[i+1]-dp/2),tk_parc[i+1]);
		tk_parc[i] = tk_parc[i+1] - 100*dp*lapserate;
		es = saturation_vapor_pressure_liquid(tk_parc[i]);
		r_parc[i] = eps*es./(100*pres[i]-es);
		# virtual temp of saturated air parcel uses saturation mixing ratio
		tkv_parc[i] = virtualtemp(tk_parc[i],r_parc[i]);
		i = i-1;
	end

	pres = pres[1:end-1];

	return tk_parc, tkv_parc,  pres, p_lcl,tk_lcl,r_parc

end


"""
moistlaps(p, tk)
Calculate moist adiabatic lapse rate (deg/Pa) based on pressure [hPa] and temperature [K]
"""
function moistlaps(p::Real, tk::Real)

	# constants
	Rd       = 287.04;      # gas constant for dry air [J/kg/K]
	Rv       = 461.5;       # gas constant for water vapor [J/kg/K]
	cpd      = 1005.7;      # specific heat dry air [J/kg/K]
	cpv      = 1870;        # specific heat water vapor [J/kg/K]
	g        = 9.80665;     # [m/s^2]
	gc_ratio = Rd/Rv;#epsilon
	p = 100*p; #to Pa


	# saturation vapor pressure [Pa]
	Tc = tk-273.15;
	#es = 611.20*exp(17.67*Tc./(Tc+243.5)); # Bolton 1980 equation 10
	es = saturation_vapor_pressure_liquid(tk);

	# latent heat of condensation [J/kg]
	L = (2.501-0.00237*Tc)*1e6; # Bolton 1980 equation 2

	# saturation mixing ratio
	rs = gc_ratio*es/(p-es);

	# density
	tk_virtual = virtualtemp(tk,rs);
	rho = p/(Rd*tk_virtual);

	# moist adiabatic lapse rateff
	malr1 = ( ( (g/cpd)*(1+rs) )/(1+(cpv/cpd*rs)) )*(1 + (L*(rs/Rd)/tk))/(1+(((L^2)*rs*(1+(rs/gc_ratio)))/(Rv*(tk^2)*(cpd+rs*cpv))));

	# convert C/m to C/Pa
	malr1 = malr1/(rho*g);
end

end