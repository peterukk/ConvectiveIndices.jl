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

function latent_heat_vaporization(tk::F) where F<:AbstractFloat
	return F(2.501e6) - (F(4.188e3)-F(1.810e3))*(tk-F(273.15))
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