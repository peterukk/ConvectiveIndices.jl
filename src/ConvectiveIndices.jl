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
saturation_vapor_pressure_liquid(T)
Calculate the saturation vapor pressure temperature [K] according to a formula presented in:
D. M. MURPHY and T. KOOP: Review of the vapour pressures of ice and supercooled water for atmospheric applications
Q. J. R. Meteorol. Soc. (2005), 131, pp. 1539-1565
"""
function saturation_vapor_pressure_liquid(T::Real)

	if T < 123 || T > 332   
		@warn "Temperature out of range [123-332] K."
    end

    temp = (54.842763f0 - (6763.22f0 / T) - 4.210f0  * log(T) + 0.000367f0  * T) + (tanh( 0.0415f0  * (T - 218.8f0 ) )*(53.878f0  - (1331.22f0 /T) - 9.44523f0 *log(T) + 0.014025f0  * T))
    es = exp(temp);

end

function saturation_vapor_pressure_liquid(T::Vector{R}) where R<:Real

    if ~(all(123 .< T) && all(T .< 332)) 
		@warn "Temperature out of range [123-332] K."
	end
	
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
	RH = 100 * e/es;
end

"""
q_to_rh( tk, p, r )
Return 1D array of relative humidity [%] given 1D arrays of temperature [K], pressure [hPa] and specific humidity 
"""
function q_to_rh(tk::Vector{F}, p::Vector{F}, q::Vector{F}) where F<:AbstractFloat
	N = length(tk)
	RH = Vector{F}(undef,N)
	c = F(18.0152)/F(28.9644);
	for i = 1:N
		es = saturation_vapor_pressure_liquid(tk[i])
		e = (q[i]*p[i]*F(100)) /(c + (1 - c)*q[i]);
		RH[i] = 100 * e/es;
	end
	return RH
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
		
    theta = Vector{F}(undef,len)
    thetae = Vector{F}(undef,len)
	thetaes = Vector{F}(undef,len)
	p_lcl = Vector{F}(undef,len)
	tk_lcl = Vector{F}(undef,len)
	r = Vector{F}(undef,len)
	rsat = Vector{F}(undef,len)

	@inbounds for i in 1:len
		#p[i] = 100*p[i]  # convert t,p to SI units
		a = ((p[i]/p0)^(-K))
		theta[i] = tk[i]*a
		es = saturation_vapor_pressure_liquid(tk[i]) # Saturation vapor pressure
		e = (rh[i]/ 100 )*es  # vapour pressure (Pa)
		
		# Calculate water vapour mixing ratio r and q specific humidity
		r[i] = (epsi*e)/(p[i]-e)
		rsat[i] = (epsi*es)/(p[i]-es)

		#ri = r[i] - rsat[i]; #liquid water content of air 
		
		# change units from g/g to g/kg
		rg = r[i]*1000; rsatg = rsat[i]*1000;

		# calculate pseudo-equivalent potential temperature, from Bolton, Mon Wea Rev, 1980
		# r = is g/kg
		# Firstly calculate Temp at LCL, note e must be in mb.
		tk_lcl[i] = ( F(2840) / (F(3.5) * log(tk[i]) - log(e/100) - F(4.805)) ) + 55;               # eqn 21, Bolton
		thetae[i] = theta[i] *  exp(( (F(3.376)/tk_lcl[i]) - F(0.00254)) * rg * (1+F(0.81) *rg *F(0.001)));   # eqn 38, Bolton
		thetaes[i] = theta[i] * exp(( (F(3.376)/tk_lcl[i]) - F(0.00254)) * rsatg * (1+F(0.81) *rg *F(0.001)));   # eqn 38, Bolton

		#LCL height using Poisson's equation
		k = F(0.2854) * (1 - F(0.28)*r[i]);
		p_lcl[i] =  F(0.01)*p[i]*((tk_lcl[i]/tk[i])^(1/k));

    end

	return theta, thetae, thetaes, p_lcl, tk_lcl, r, rsat

end


"""
	calc_CAPE_thetae (ps[hPa], tks[K], qs[kg/kg], zs[m], parcel, dp_mix[hPa], dp[hPa], kiss) 

Calculate convective indices such as **CAPE** and **CIN** for a parcel of choice (surface/most unstable, with/without vertical mixing).

# Examples

```jldoctest 
julia> LI,CAPE,CIN = calc_CAPE_theta(ps,tks,qs,zs,sp, parcel = 2, dp_mix = 100, kiss= 1)
(-27.416924139871526, 4428.182537242374, 137.85516940477973)
julia> LI, CAPE, CIN, pLCL, zBCL, CAPECIN_ALCL, CIN_LCL, MRH_ALCL, MRH1, MRH2 = calc_CAPE_thetae(ps,tks,qs,zs)
(-1.6502346944216129, 120.80558885439602, 23.64198824254466, 787.8515322945883, 351.837890625, -23.998722156796717, 0, 63.845851443325564, 76.3582759152618, 56.28591549989976)
```

**OUTPUT** by default, following Float32 values are returned (kiss=0): \n
Lifted Index [°C], CAPE [J/kg], CAPE-CIN above the LCL [J/kg], MRH (mean RH%) above the LCL [%], CIN below LCL [J/kg], MRH 600-800 hPa, MRH 300-600 HP, LCL [hPa], CIN [J/kg] \n
Toggle kiss=1 to only return CAPE, Lifted Index and CIN. \n

**INPUT**:
(N-element ARRAYs) **ps**,**tks**,**qs**,**zs** = vertical profiles of pressure, temperature, specific humidity and geopotential height \n
OPTIONAL keyword arguments: \n
`parcel = 1` : the most unstable parcel in the lowest 350 hPa (default) \n
`parcel = 2` : surface parcel \n
`dp_mix = 0...100` : pressure layer depth [hPa] for mixing the source parcel (default 50, use 0 for no mixing) \n
`kiss = 1`: keep it simple, stupid - output only CAPE, LI, and CIN (default 0). \n

This routine uses a THETA-E formulation for all indices (similarly to ECMWF CAPE), thereby skipping explicit parcel computations. 
This results in different values (e.g. for CAPE, 30% larger) than classic computations, but in no worse correlation with observed convection [1].

TIP: Use `parcel=1` and `dp_mix=50` for a hybrid mixed-layer most-unstable parcel similar to the one used by ECMWF. The MLMU-Lifted Index was the overall
thunderstorm predictor in Europe in [1].

[1] Ukkonen et al. (2018)
"""
function calc_CAPE_thetae(ps::Vector{F},tks::Vector{F},qs::Vector{F},zs::Vector{F}; parcel::Integer=1,dp_mix::Real=50,dp::Real=5,kiss::Integer=0) where F<:AbstractFloat

	g = F(9.80665)

	#Check that profiles are descending, otherwise flip the arrays
	if ps[end] < ps[1]
		tks = reverse(tks,1); rhs = reverse(rhs,1); ps = reverse(ps,1)
	end

	sp = ps[end]
	
	rhs = q_to_rh(tks,ps,qs)
	
	pres = collect(ceil(ps[1]):dp:floor(sp))
	nlevs = size(pres,1)

	# Interpolation: use "Interpolations" package for linear interpolation on a non-uniform grid
	knots = (ps,);

	itp = interpolate(knots,tks,Gridded(Linear())) 
	tk_env = itp(pres)
	itp = interpolate(knots,rhs,Gridded(Linear()))
	rh_env = itp(pres)
	itp = interpolate(knots,zs,Gridded(Linear()))
	z_env = itp(pres)

	theta,thetae,thetaes = thermo_rh(tks,100*ps,rhs)	
	itp = interpolate(knots,thetae,Gridded(Linear())); thetae_env = itp(pres)
	itp = interpolate(knots,thetaes,Gridded(Linear())); thetaes_env = itp(pres)

	# FIND PARCEL
	if parcel == 1 #	--> MOST UNSTABLE
		# Find level with the highest theta-e in the lowest 350 hPa
		thetae_env[pres.<(sp - 350)] .= 0.0; # Set elements above the lowest 350 hPa to 0
		thetae_max,iPL = findmax(thetae_env); #Index of parcel level
	else #				--> SURFACE PARCEL
		iPL = nlevs # 
	end

	# PARCEL MIXING WITH ENVIRONMENT
	# Mix the properties over depth specified by dp_mix (e.g. 50 hPa)
	@views if dp_mix > 0
		#how many indices correspond to the mixing depth divided by two? (50 hPa mixing depth = 25 below, 25 above)
		dind = floor(Int,(dp_mix/2)/dp)
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
	i500 = findmin(abs.(pres .- 500))[2]
	LI = thetaes_env[i500] - thetae_parc

	# CAPE AND CIN
	#the quantity being integrated in theta-e formulation for CAPE
	tdiff = (thetae_parc .- thetaes_env)./thetaes_env
	tdiff_cape = copy(tdiff)
	#Only positive regions in sounding contribute to CAPE
	tdiff_cape[tdiff_cape.<0] .= 0

	iLCL = findlast(pres.<pLCL)
	iLFC = @views findlast(tdiff[1:iLCL].>0)
	iEL = @views findfirst(tdiff.>0)

	if iLFC == nothing
		CAPE = F(0)
		CIN = NaN
		iLFC = NaN
	else
		if (iLFC == iEL)
			iEL=iEL-1;
		end
		CAPE = -g*trapz(z_env[iEL:iLFC],tdiff_cape[iEL:iLFC])
		CIN = g*trapz(z_env[iLFC:iPL],tdiff[iLFC:iPL])
	end
	
	# RETURN EARLY?
	if kiss == 1 #Stop here and return with only LI,CAPE,CIN 
		return LI,CAPE,CIN
	end
	
	# CAPE and RH in the 250-hPa depth layer ABOVE LCL	
	idp_250 = floor(Int,250/dp)
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

	i300 = findmin(abs.(pres .- 300))[2]
	i600 = findmin(abs.(pres .- 600))[2]
	i800 = findmin(abs.(pres .- 800))[2]

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

	@inbounds for i = 1:N
		qsats[i] = (100/rhs[i]) * qs[i]
		qmix[i] = @views mean(qs[i:N])
	end
	
	z = vcat( collect(10000:-300:6100), collect(6000:-100:3000), collect(3000:-20:zs[end]) )

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
	r_parc = zeros((n,)); #mixing ratio
	tk_parc = zeros((n,));  
	tkv_parc = zeros((n,)); 

	i_lcl = findlast(pres.<p_lcl);

	#From p0 to LCL, use dry adiabatic lapse rate:
	r_parc[i_lcl:end] = r0;
	tk_parc[i_lcl:end] = tk0 *(pres[i_lcl:end-1]./p0).^0.286;
	tkv_parc[i_lcl:end] = virtualtemp(tk_parc[i_lcl:end],r_parc[i_lcl:end]);

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