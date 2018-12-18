# ConvectiveIndices.jl
Julia package for calculating convective indices (e.g. CAPE) from atmospheric sounding data.

The core function is ```calc_CAPE_thetae```, which outputs parameters such as CAPE, Lifted Index and CIN from input columns of pressure, temperature, specific humidity and geometric height.

```
help?> calc_CAPE_thetae
search: calc_CAPE_thetae


  Calculate convective indices such as CAPE and CIN for a parcel of choice (surface/most unstable, with/without vertical mixing).

  Examples
  ≡≡≡≡≡≡≡≡≡≡

  julia> LI,CAPE,CIN = calc_CAPE_theta(ps,tks,qs,zs,sp, parcel = 2, dp_mix = 100, kiss = 1)
  (-27.416924139871526, 4428.182537242374, 137.85516940477973)
  julia> LI, CAPE, CIN, pLCL, zBCL, CAPECIN_ALCL, CIN_LCL, MRH_ALCL, MRH1, MRH2 = calc_CAPE_thetae(ps,tks,qs,zs)
  (-1.6502346944216129, 120.80558885439602, 23.64198824254466, 787.8515322945883, 351.837890625, -23.998722156796717, 0, 63.845851443325564, 76.3582759152618, 56.28591549989976)


  OUTPUT: by default, following Float32 values are returned (kiss=0):

  Lifted Index [°C], CAPE [J/kg], CAPE-CIN above the LCL [J/kg], MRH (mean RH%) above the LCL [%], CIN below LCL [J/kg], MRH 600-800 hPa, MRH 300-600 HP, LCL [hPa], CIN [J/kg]

  Toggle kiss=1 (Keep It Simple, Stupid) to only return CAPE, Lifted Index and CIN.

  INPUT: (N-element ARRAYs) ps,tks,qs,zs = vertical profiles of pressure [hPa], temperature [K], specific humidity [kg/kg] and geopotential height [m].

  OPTIONAL keyword arguments:

  parcel = 1 : the most unstable parcel in the lowest 350 hPa (default)

         = 2 : surface parcel

  dp_mix = 0...100 : pressure layer depth [hPa] for mixing the source parcel (default 50, use 0 for no mixing)

  kiss = 0 : full outputs
       = 1 : only CAPE, LI and CIN.

  This routine uses a THETA-E formulation for all indices (similarly to how CAPE is calculated in ECMWFs IFS model), thereby skipping explicit parcel computations. This results in different values (for CAPE, roughly 30% larger) than
  classic computations, but in no worse correlation with observed convection [1].

  TIP: Use parcel=1 and dp_mix=50 for a hybrid mixed-layer most-unstable parcel similar to the one used by ECMWF. The MLMU-Lifted Index was the best overall thunderstorm predictor in Europe in [1].

  [1] Ukkonen and Mäkelä (2018): Evaluation of machine learning classifiers for predicting deep convection, JAMES.
```

The code is computationally efficient yet simple and easy to modify thanks to Julias language design. 
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
