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

using Interpolations
using LinearAlgebra
using Statistics

"""
numerical integration of y (Vector) with respect to the coordinates specified by x using the trapezoidal rule
"""

function trapz(x::Vector{R}, y::Vector{R}) where R<:Real
    local len = length(y)
    if (len != length(x))
        error("Vectors must be of same length")
    end
    r = 0.0
    for i in 2:len
       r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    r/2.0
end


"""
saturation_vapor_pressure_liquid(T)
Calculate the saturation vapor pressure temperature [K] according to a formula presented in:
D. M. MURPHY and T. KOOP: Review of the vapour pressures of ice and supercooled water for atmospheric applications
Q. J. R. Meteorol. Soc. (2005), 131, pp. 1539-1565
"""
function saturation_vapor_pressure_liquid(T::Number)

        if ~(all(123 .< T) && all(T .< 332))
            @warn "Temperature out of range [123-332] K."
        end

        temp = (54.842763 - (6763.22./ T) - 4.210 * log(T) + 0.000367 * T) + (tanh( 0.0415 * (T - 218.8) )*(53.878 - (1331.22/T) - 9.44523*log(T) + 0.014025 * T))
        es = exp(temp);

end

"""
moistlaps(p, tk)
Calculate moist adiabatic lapse rate (deg/Pa) based on pressure [hPa] and temperature [K]
"""
function moistlaps(p::Number, tk::Number)

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

"""
virtualtemp( tk, r )
Return virtual temperature given temperature in K and mixing ratio in g/g
"""
function virtualtemp(tk, r)

	gc_ratio =  287.04/461.5;#epsilon=Rd/Rv, where Rd=gas const. for dry air and Rv=gas const. for water vap.
	tk_v = tk.*(1.0+r/gc_ratio)./(1.0+r); #tk*(1 + 0.608*r);

end

"""
rh_to_r( tk, p, r )
Return mixing ratio [kg/kg] given temperature [K], pressure [hPa] and relative humidity [%] 
"""
function rh_to_r(tk::Number, p::Number, RH::Number)
	local c = 18.0152/28.9644;
	es = saturation_vapor_pressure_liquid(tk);
	e = (RH/100)*es;
	Pd = p*100 - e;
	r = (e/Pd)*c
end

"""
r_to_rh( tk, p, r )
Return relative humidity [%] given temperature [K], pressure [hPa] and mixing ratio [kg/kg]
"""
function r_to_rh(tk::Number, p::Number, r::Number)
	es = saturation_vapor_pressure_liquid(tk);
	local c = 18.0152/28.9644;
	e = (r*p*100)/(r + c);
	RH = 100 * e/es;
end

"""
q_to_rh( tk, p, r )
Return relative humidity [%] given temperature [K], pressure [hPa] and specific humidity 
"""
function q_to_rh(tk::Number, p::Number, q::Number)
	es = saturation_vapor_pressure_liquid(tk);
	local c = 18.0152/28.9644;
	e = (q*p*100)/(c + (1 - c)*q);
	RH = 100 * e/es;
end

"""
dewpoint_to_q( tk, p, td )
Return specific humidity [kg/kg] given temperature [K], pressure [hPa] and dewpoint temperature [K]
"""
function dewpoint_to_q(tk::Number, p::Number, td::Number)
	local c = 18.0152/28.9644;
	e = saturation_vapor_pressure_liquid(td);
	Pd = p*100 - e;
	q = (e*c)/(e*c+Pd);
end



"""
function (theta, thetae, thetaes, p_lcl, tk_lcl, r) = thermo_rh(tk,p,rh)
thermo_rh generates thermodynamic variables from t(K) p(hPa) rh(%) 
output [theta, thetae, r, tk_lcl, p_lcl] in K, K, kg/kg, K, hPa; r = water vapour mixing ratio
"""
function thermo_rh(tk::Number,p::Number,rh::Number)

	# changed rsat denominator from p to p-es (3 may 07)
	# Realised thetae calculation is incorrect (17 Nov 2009)
	# Changed to used formulae of Bolton (1980), for pseudo-equivalent p.t.

	#convert t,p to SI units
	p = p*100; #Pa

	#calculate the potential temperature from temperature array
	p0=1000*100;	#reference pressure in Pa
	R=287;		#gas constant
	cp=1004;	#specific heat wrt pressure
	K= R/cp;
	a = ((p/p0)^(-K));
	theta = tk*a;

	#calculate equivalent potential temperature
	#Lv = 2.5e6; 	#latent heat of vapourisation
	eps = 0.62198; 	#Rd/Rv

	es = saturation_vapor_pressure_liquid(tk); # Saturation vapor pressure
	e = ((rh/100)*es);                              # vapour pressure (Pa)

	# Calculate water vapour mixing ratio r and q specific humidity
	r = (eps*e)/(p-e);
	rsat = (eps*es)/(p-es);
	#ri = r - rsat; #liquid water content of air (Zhang et al 1990)
	k = 0.2854 * (1 - 0.28*r);
	# change units from g/g to g/kg
	rg = r*1e3; rsat = rsat*1e3;

	# calculate pseudo-equivalent potential temperature, from Bolton, Mon Wea Rev, 1980
	# r = is g/kg
	# Firstly calculate Temp at LCL, note e must be in mb.
	tk_lcl = ( 2840 / (3.5 * log(tk) - log(e/100) - 4.805) ) + 55;               # eqn 21, Bolton
	thetae = theta*exp(( (3.376/tk_lcl) - 0.00254)*rg*(1+0.81 * rg *1e-3));   # eqn 38, Bolton
	thetaes = theta*exp(( (3.376/tk_lcl) - 0.00254)*rsat*(1+0.81 * rg *1e-3));   # eqn 38, Bolton

	#LCL height using Poisson's equation
	p_lcl =  0.01*p.*((tk_lcl./tk).^(1.0/k));

	return theta, thetae, thetaes, p_lcl, tk_lcl, r 

end

"""
function (theta, thetae, thetaes, p_lcl, tk_lcl, r) = thermo_rh(tks,ps,rhs)
thermo_rh generates thermodynamic variables from t(K) p(hPa) rh(%) N-element arrays 
output [theta, thetae, r, tk_lcl, p_lcl] in K, K, kg/kg, K, hPa; r = water vapour mixing ratio
"""
function thermo_rh(tk::Vector{F},p::Vector{F},rh::Vector{F}) where F<:AbstractFloat

	# changed rsat denominator from p to p-es (3 may 07)
	# Realised thetae calculation is incorrect (17 Nov 2009)
	# Changed to used formulae of Bolton (1980), for pseudo-equivalent p.t.

	#convert t,p to SI units
	p = p.*100; #Pa

	local p0 = 1000*100;	#reference pressure in Pa
	local R = 287;		#gas constant
	local cp = 1004;	#specific heat wrt pressure
	local K = R./cp;
	local eps = 0.62198;

    local len = length(tk)
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

    for i in 1:len
		a = ((p[i]/p0).^(-K));
		theta[i] = tk[i]*a;
		es = saturation_vapor_pressure_liquid(tk[i]); # Saturation vapor pressure
		e = (rh[i]/100)*es;   # vapour pressure (Pa)
		
		# Calculate water vapour mixing ratio r and q specific humidity
		r[i] = (eps*e)./(p[i]-e);
		rsat[i] = (eps*es)./(p[i]-es);

		#ri = r[i] - rsat[i]; #liquid water content of air 
		k = 0.2854 * (1 - 0.28*r[i]);
		# change units from g/g to g/kg
		rg = r[i]*1e3; rsatg = rsat[i]*1e3;

		# calculate pseudo-equivalent potential temperature, from Bolton, Mon Wea Rev, 1980
		# r = is g/kg
		# Firstly calculate Temp at LCL, note e must be in mb.
		tk_lcl[i] = ( 2840 / (3.5 * log(tk[i]) - log(e/100) - 4.805) ) + 55;               # eqn 21, Bolton
		thetae[i] = theta[i] *  exp(( (3.376/tk_lcl[i]) - 0.00254) * rg * (1+0.81 *rg *1e-3));   # eqn 38, Bolton
		thetaes[i] = theta[i] * exp(( (3.376/tk_lcl[i]) - 0.00254) * rsatg * (1+0.81 *rg *1e-3));   # eqn 38, Bolton

		#LCL height using Poisson's equation
		p_lcl[i] =  0.01*p[i]*((tk_lcl[i]/tk[i])^(1/k));

    end

	return theta, thetae, thetaes, p_lcl, tk_lcl, r, rsat

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
function liftparcel(tk0::Number, p0::Number, rh0::Number, dp::Number=5)

	local eps = 0.62198; 	#Rd/Rv

	#Calculate LCL height in hPa p_lcl
	_,_,_,p_lcl,tk_lcl,r0 = thermo_rh(tk0,p0,rh0);

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
	calc_CAPE_thetae (ps[hPa], tks[K], qs[kg/kg], zs[m], parcel, dp_mix[hPa], dp[hPa], kiss) 

Calculate convective predictors such as **CAPE** and **CIN** for a parcel of choice (surface/most unstable, with/without vertical mixing).

# Examples

```jldoctest 
julia> LI,CAPE,CIN = calc_CAPE_theta(ps,tks,qs,zs,sp, parcel = 2, dp_mix = 100, kiss= 1)
(-27.416924139871526, 4428.182537242374, 137.85516940477973)
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
function calc_CAPE_thetae(ps::Vector{F},tks::Vector{F},qs::Vector{F},zs::Vector{F}; parcel::Integer=1,dp_mix::Number=50,dp::Number=5,kiss::Integer=0) where F<:AbstractFloat

	local g = 9.80665

	#Check that profiles are descending, otherwise flip the arrays
	if ps[end] < ps[1]
		tks = reverse(tks,1); rhs = reverse(rhs,1); ps = reverse(ps,1)
	end

	sp = ps[end]
	
	rhs = q_to_rh.(tks,ps,qs)
	
	pres = collect(100:dp:sp)
	nlevs = size(pres,1)

	# Interpolation: use "Interpolations" package for linear interpolation on a non-uniform grid
	knots = (ps,);


	itp = interpolate(knots,tks,Gridded(Linear())); tk_env = itp(pres)
	itp = interpolate(knots,rhs,Gridded(Linear())); rh_env = itp(pres)
	itp = interpolate(knots,zs,Gridded(Linear())); z_env = itp(pres)

	theta,thetae,thetaes = thermo_rh(tks,ps,rhs)
	
	# !!!!! Still to add theta_def?
	
	itp = interpolate(knots,thetae,Gridded(Linear())); thetae_env = itp(pres)
	itp = interpolate(knots,thetaes,Gridded(Linear())); thetaes_env = itp(pres)

	# FIND PARCEL
	if parcel == 1 #	--> MOST UNSTABLE
		# Find level with the highest theta-e in the lowest 350 hPa
		thetae_env[pres.<(sp - 350)] = 0; # Set elements above the lowest 350 hPa to 0
		thetae_max,iPL = findmax(thetae_env); #Index of parcel level
	else #				--> SURFACE PARCEL
		iPL = nlevs # 
	end

	# PARCEL MIXING WITH ENVIRONMENT
	# Mix the properties over depth specified by dp_mix (e.g. 50 hPa)
	if dp_mix > 0
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
	_,thetae_parc,_,pLCL = thermo_rh(tk0,p0,rh0)

	# LIFTED INDEX
	i500 = findmin(abs.(pres - 500))[2]
	LI = thetaes_env[i500] - thetae_parc

	# CAPE AND CIN
	#the quantity being integrated in theta-e formulation for CAPE
	tdiff = (thetae_parc - thetaes_env)./thetaes_env
	tdiff_cape = copy(tdiff)
	#Only positive regions in sounding contribute to CAPE
	tdiff_cape[tdiff_cape.<0] = 0


	iLCL = findlast(pres.<pLCL)
	iLFC = findlast(tdiff[1:iLCL].>0)
	iEL = findfirst(tdiff.>0)

	CAPE = 0
	CIN = NaN
	if iLFC>0
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
		MRH_ALCL = mean(rh_env[iLCL-idp_250:iLCL])
		CAPECIN_ALCL = -g*trapz(z_env[iLCL-idp_250:iLCL],tdiff[iLCL-idp_250:iLCL])
	else
		MRH_ALCL = mean(rh_env[1:iLCL])
		CAPECIN_ALCL = -g*trapz(z_env[1:iLCL],tdiff[1:iLCL])
	end

	CIN_LCL = NaN
	if iLCL < iPL
		CIN_LCL = g*trapz(z_env[iLCL:iPL],tdiff[iLCL:iPL])
	end

	i300 = findmin(abs.(pres - 300))[2]
	i600 = findmin(abs.(pres - 600))[2]
	i800 = findmin(abs.(pres - 800))[2]

	# MEAN RELATIVE HUMIDITIES 
	MRH1 =  mean(rh_env[i600:i800])
	MRH2 =  mean(rh_env[i300:i600])
	
	# INHIBITON
	# Calculate the "Buoyant Condensation Level" (BCL) which is "the level at which saturation would occur
	# through buoyant mixing alone due to sensible heating from the surface"
	# This is an inhibition measure by Tawfik et al. (2014):
	#  A process-based framework for quantifying the atmospheric preconditioning of surface-triggered convection
	
	zBCL = calc_BCL(qs,rhs,zs)
	
	return LI,CAPE,CAPECIN_ALCL,MRH_ALCL,CIN_LCL,MRH1,MRH2,pLCL,zBCL,CIN

end

"""
calc_BCL
Calculate the 'Buoyant Condensation Level' (**BCL**) which is 'the level at which saturation would occur
through buoyant mixing alone due to sensible heating from the surface'

This is an inhibition measure for quantifying how preconditioned the atmosphere is to moist convection. \n
Described in Tawfik et al. (2014):
A process-based framework for quantifying the atmospheric preconditioning of surface-triggered convection
"""
function calc_BCL(qs,rhs,zs)

	qsats = 100.0/rhs .* qs;
	qm = repeat(qs,1,length(qs))
	tril!(qm);  # zeros above diagonal (upper right corner)
	qm[qm.==0] = NaN; #replace with NaNs

	#This will be unnecessary in Julia 1.0
	nanmean(x) = mean(filter(!isnan,x))
	nanmean(x,y) = mapslices(nanmean,x,y)
	#Calculate mixed-layer means for different depths: qm(i) = mean(qm_sfc,..qmi-1,qmi)
	qmix = nanmean(qm,1)[:]
	
	z = collect(zs[1]:-10:zs[end]) 

	knots = (reverse(zs,1),);
	itp = interpolate(knots,reverse(qmix,1),Gridded(Linear())); qmix_int = itp[z]
	itp = interpolate(knots,reverse(qsats,1),Gridded(Linear())); qsat_int = itp[z]
	iBCL = findlast(qmix_int.>qsat_int)
	
	zBCL = z[iBCL]
	return zBCL
	
	#itp = interpolate(knots,thetas,Gridded(Linear())); thetas_int = itp[z]
	#theta_BM = thetas_int(iBCL);
	#theta_def = theta_BM - thetas(end);

end


end

