module ConvectiveIndices
	using Interpolations,
	LinearAlgebra,
	Statistics

	export saturation_vapor_pressure_liquid,
	moistlaps,
	virtualtemp,
	liftparcel,
	thermo_rh,
	rh_to_r,
	r_to_rh,
	q_to_rh,
	dewpoint_to_q,
	calc_CAPE_thetae,
	calc_BCL,
	calc_entropy,
	calc_entropy_inverse,
	calc_dilute_dCAPE,
	calc_dilute_CAPE

	include("thermodyn.jl")
	include("utils.jl")
	include("main_functions.jl")
end