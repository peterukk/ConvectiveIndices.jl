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
returns the `diluted CAPE` using a constant mass entrainment rate. The code is adapted from the subroutine buoyan_dilute (by Richard Neale) 
in the CAM model. The effect of freezing/condensation on parcel temperature can be accounted for by toggling freezing=1. This is about 
twice as slow and the resulting dilute CAPE actually had lower correlation with observed thunderstorms in Sri Lanka (than using freezing=0).

The parcel is launched from the level with the highest theta-e in the lowest 350 hPa of the atmosphere. Optionally, the index can be specified by the user: 
calc_dilute_CAPE(ps,tks,qs,zs,parcel_index=30), would use a surface parcel if length(ps) = 30 and the last level is the surface.

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
	#println("ilaunch: ",ilaunch)

	# Initialize values
	N = length(ps)

	g = F(9.80665)			# Gravity
	R = F(287.04)			# gas constant
	cpliq =	F(4.188e3)		# J/kg/K specific heat of liquid water (fresh)
	tfreez = F(273.15) 		# K 
	tscool = F(0)			# super cooled temperature offset (in degC) (eg -35).
	latice = F(3.337e5)		# specific heat of fusion
	lwmax = F(1.e-3)		# Maximum condesate that can be held in cloud before rainout.
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
    # Calculate tempereature and humidity from entropy - this step is unnecessary, same values returned?
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
		#dmpdp = 0
		
		# Sum entrainment to current level
		# entrains q,s out of intervening dp layers, in which linear variation is assumed
		# so really it entrains the mean of the 2 stored values.
        sp  = sp  - dmpdp*dp*s_env 
        qtp = qtp - dmpdp*dp*qt_env 
		mp  = mp  - dmpdp*dp
		 
		 # Entrain s and qt to next level.
        s_mix[i]  = (sp0  +  sp) / (mp0 + mp)
		qt_mix[i] = (qtp0 + qtp) / (mp0 + mp)
	
		# The new parcel entropy and total water due to mixing are now calculated

		# Invert entropy from s and q to determine T and qsat of mixture 
		# t[i,k] used as a first guess so that it converges faster.

        tk_guess = tk_mix[i+1]
		tk_mix[i],qsat_mix[i] =  calc_entropy_inverse(s_mix[i],ps[i],qt_mix[i],tk_guess) 
		 
		 #tk_mix_v = tk_mix[i]
		 #tk_env_v = tks[i])
		 
		 #println("p, tmix, qmix")
		 #println([p_env,tk_mix[i],qsat_mix[i]])
	end
	#return tk_mix,qsat_mix,qt_mix

	if freezing==1	# ACCOUNT FOR FREEZING/CONDENSATION EFFECT ON TEMPERATURE

	# So we now have a profile of entropy and total water of the entraining parcel Varying with height 
	# from the launch level klaunch parcel=environment. To the top allowed level for the existence of convection.
	# Now we have to adjust these values such that the water held in vaopor is < or = to qsmix. Therefore, we assume
	# that the cloud holds a certain amount of condensate (lwmax) and the rest is rained out (xsh2o). This, obviously 
	# provides latent heating to the mixed parcel and so this has to be added back to it. 
	# But does this also increase qsmix as well? Also freezing processes

		xsh2o = zero(tks)		# Condensate lost from parcel
		ds_xsh2o = zero(tks)	# Entropy change due to loss of condensate
		ds_freeze = zero(tks)	# Entropy change dut to freezing of precipitation

		# ! Set parcel values at launch level assume no liquid water. 
		qp = zero(tks) 			# Parcel water vapour (saturated value above lcl)
		qp[ilaunch] = qtp0		# Just the wv at launch level

		# tk_v[ilaunch] = virtualtemp(tk_mix[ilaunch],qp[ilaunch])
		new_s = F(0)
		new_q = F(0)

		for i in range(ilaunch-1,1,step=-1)
			# Iterate 2 times for s,qt changes
			for j in [1,2]
				# Rain (xsh2o) is excess condensate, bar LWMAX (Accumulated loss from qtmix).
				
				xsh2o[i] = maximum( [F(0),(qt_mix[i] - qsat_mix[i] - lwmax)] )

				# contribution to entropy from precip loss of condensate  (Accumulated change from smix).(-ve)
				ds_xsh2o[i] = ds_xsh2o[i] - cpliq * log(tk_mix[i]/tfreez) * maximum( [F(0),(xsh2o[i] - xsh2o[i+1])] )
				
				# Entropy of freezing: latice times amount of water involved divided by T
				 
				if ( (tk_mix[i] <= tfreez+tscool) & (ds_freeze[i+1] == F(0)) ) # One off freezing of condensate. 
					ds_freeze[i] = (latice/tk_mix[i]) * maximum( [F(0),qt_mix[i]-qsat_mix[i]-xsh2o[i]] ) # Gain of LH
				end
				 
				if ( (tk_mix[i] <= tfreez+tscool) & (ds_freeze[i+1] != F(0)) ) #  Continual freezing of additional condensate.
					ds_freeze[i] = ds_freeze[i+1]+(latice/tk_mix[i]) * maximum( [F(0),(qsat_mix[i+1]-qsat_mix[i])] )
				end

				# Adjust entropy and accordingly to sum of ds (be careful of signs).
				new_s = s_mix[i] + ds_xsh2o[i] + ds_freeze[i] 
				# Adjust liquid water
				new_q = qt_mix[i] - xsh2o[i]

				# Invert entropy to get updated tk_mix and qs_mix of parcel
				tk_guess = tk_mix[i]

				tk_mix[i],qsat_mix[i] =  calc_entropy_inverse(new_s,ps[i],new_q,tk_guess) 

			end # Iteration loop for freezing processes
	
			if (new_q > qsat_mix[i]) # Super-saturated so condensate present - reduces buoyancy.
				qp[i] = qsat_mix[i]
			 else                    # Just saturated/sub-saturated - no condensate virtual effects
				qp[i] = new_q 		 # in this case new_q is simply the old qt_mix (xsh2o is 0)
			end
			
		end		# vertical loop for freezing processes

	else # No freezing/condensation processes
		qp = qt_mix # parcel water vapour = total water of the entraining parcel
	end 

	# Virtual temperature
	tiedke_add = F(0)
	tk_parcel_v = virtualtemp.(tk_mix, (qp ./ (1 .- qp)))
	tk_env_v = virtualtemp.(tks, (qs ./ (1 .- qs)))
	# Buoyancy
	tdiff = (tk_parcel_v .+ tiedke_add .- tk_env_v)./tk_env_v

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
	return CAPE #,ilaunch,iLCL,iLFC,iEL #,CIN,tdiff
end


function calc_dilute_dCAPE(ps::Vector{F},tks::Vector{F},qs::Vector{F},zs::Vector{F},adv_tk::Vector{F},adv_q::Vector{F}) where F<:AbstractFloat

	g = F(9.80665)
	sp = ps[end]
	rhs = q_to_rh(tks,ps,qs)

	#CAPE calculated "normally", in this case diluted CAPE:

	theta,thetae = thermo_rh(tks,100*ps,rhs)
	thetae[ps.<(sp - 350)] .= 0; # Set elements above the lowest 350 hPa to 0
	thetae_max,ilaunch = findmax(thetae); #Index of parcel level

	# Calculate dilute CAPE using the previously determined parcel (most unstable in the lowest 350 hPa)
	CAPE = calc_dilute_CAPE(ps,tks,qs,zs,parcel_index=ilaunch)

	# Store the parcel values
	tk0 = tks[ilaunch]; q0 = qs[ilaunch];

	# in dCAPE, the environmental (but only environmental) temp and humidity gain an increment from 
	tks = tks .+ adv_tk 
	qs = qs .+ adv_q 
	# adv_T is (u*dT/dx) = [ m/s * K/m] = [K/s]. Multiply by 3600 to get hourly increment

	# The advective increment is in units of temperature or moisture per second.
	# if multiplied my e.g. 3600 to get the hourly change, the updated q values are sometimes far too small, resulting in negative relative hum
 	tks = tks .+ adv_tk 
	qs = qs .+ adv_q 
	tks[ilaunch] = tk0; qs[ilaunch] = q0

	CAPE_new = calc_dilute_CAPE(ps,tks,qs,zs,parcel_index=ilaunch)

	dCAPE = (CAPE_new - CAPE)*3600 # multiply by 3600 here instead to get the hourly change in CAPE due to large-scale advection
	return dCAPE,CAPE
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