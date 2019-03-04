# ConvectiveIndices.jl
A Julia package for calculating convective indices (e.g. CAPE) from atmospheric sounding data.

The core function is ```calc_CAPE_thetae```, which outputs parameters such as CAPE, Lifted Index and CIN from input columns of pressure, temperature, specific humidity and geometric height, for a user-defined parcel (surface/most unstable, with/without vertical mixing).

```
help?> calc_CAPE_thetae
search: calc_CAPE_thetae

calc_CAPE_thetae (ps[hPa], tks[K], qs[kg/kg], zs[m], parcel, dp_mix[hPa], dp_intp[hPa], kiss)

Calculate convective indices such as CAPE and CIN for a parcel of choice (surface/most unstable, with/without vertical mixing). Uses a theta-e formulation which increases absolute values.

  Examples
  ≡≡≡≡≡≡≡≡≡≡

julia> LI, CAPE, CIN = calc_CAPE_theta(ps,tks,qs,zs) # most unstable parcel, mixed over 50 hPa (default)
  (-8.94582333506736, 1613.7159227760612, 327.257167221434))
julia> LI, CAPE, CIN = calc_CAPE_theta(ps,tks,qs,zs, parcel = 2, dp_mix = 0) # surface parcel, not mixed
  (-12.416924139871522, 2428.182537242374, 85.85516940477973) 
julia> LI, CAPE, CIN, pLCL, zBCL, CAPECIN_ALCL, CIN_LCL, MRH_ALCL, MRH1, MRH2 = calc_CAPE_thetae(ps,tks,qs,zs, FULL = 1) # full calculations
  (-8.94582333506736, 1613.7159227760612, 327.257167221434, 936.6429885118564, 1230.0, -189.68905798724995, 128.5705360872618, 69.90722164805184, 56.290565968008316, 30.494525283693054)


OUTPUT for FULL=1:

Lifted Index [°C], CAPE [J/kg], CIN [J/kg], pLCL [hPa], Buoyant Condensation Level [m], CAPE-CIN above the LCL [J/kg], CIN below LCL [J/kg], MRH (mean RH%) above the LCL [%], MRH 600-800 hPa, MRH 300-600 hPa

INPUT: (N-element ARRAYs) ps,tks,qs,zs = vertical profiles of pressure, temperature, specific humidity and geopotential height

  OPTIONAL keyword arguments:

  parcel = 1 : the most unstable parcel in the lowest 350 hPa (default)

  parcel = 2 : surface parcel, or parcel from the lowest level

  dp_mix = 0...100 : pressure layer depth [hPa] for mixing the source parcel (default 50, use 0 for no mixing)

  dp_intp = 5 linearly interpolate to a uniform pressure grid with resolution dp (default 5). Use 0 to skip. This is a lot faster, but disables mixing and FULL option, and is not recommended for low-resolution input.

  FULL = 1: Full calculations to include also less known convective predictors, which were found useful in [1].

This routine uses a equivalent potential temperature formulation for all indices (similarly to ECMWF CAPE), avoiding vertical loops altogether. This results in larger absolute values (e.g. for CAPE, 30% larger) than classic computations, but in no worse correlation with observed convection [1].

NOTE: Default option parcel=1 and dp_mix=50 corresponds to a hybrid mixed-layer most-unstable parcel similar to the one used by ECMWF. The MLMU-Lifted Index was the overall thunderstorm index for Europe in [1].

  [1] Ukkonen and Mäkelä (2019): Evaluation of machine learning classifiers for predicting deep convection, JAMES.
```

The code is fast yet easy to read and modify thanks to Julias language design. 
Processing 6.2 million reanalysis pseudosoundings on an 8-core CPU:

```
@time @sync @distributed for i = 1:nlon
               for j = 1:nlat
                       for k = 1:ntime
                               tks = tk[i,j,:,k]; qs = q[i,j,:,k]; ps = p[i,j,:,k]; zs = z[i,j,:,k];
                                LI[i,j,k],CAPE[i,j,k],CIN[i,j,k],PLCL[i,j,k],ZBCL[i,j,k],CAPECIN_ALCL[i,j,k], CIN_LCL[i,j,k],        
                                MRH_ALCL[i,j,k],MRH1[i,j,k],MRH2[i,j,k]  = calc_CAPE_thetae(ps,tks,qs,zs,dp=2);
                       end
               end
       end
208.778669 seconds (3.02 M allocations: 148.865 MiB, 0.03% gc time)
```
