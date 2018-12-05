# ConvectiveIndices.jl
Julia package for calculating convective indices (e.g. CAPE) from atmospheric sounding data.


Processing 6.2 million pseudosoundings on an 8-core CPU:

@time @sync @distributed for i = 1:nlon
               for j = 1:nlat
                       for k = 1:ntime
                               tks = tk[i,j,:,k]; qs = q[i,j,:,k]; ps = p[i,j,:,k]; zs = z[i,j,:,k];
                               LI[i,j,k],CAPE[i,j,k],CIN[i,j,k],PLCL[i,j,k],ZBCL[i,j,k],CAPECIN_ALCL[i,j,k],CIN_LCL[i,j,k],MRH_ALCL[i,j,k],MRH1[i,j,k                                  ],MRH2[i,j,k] = calc_CAPE_thetae(ps,tks,qs,zs,dp=2);
                       end
               end
       end
208.778669 seconds (3.02 M allocations: 148.865 MiB, 0.03% gc time)
