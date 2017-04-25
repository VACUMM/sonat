# Creates a subsample of ensemble file created by test_ens
ncks -O -v temp,temp_surf,u_surf,v_surf -d member,,4 -d lon_b,,,2 -d lat_b,,,2  ../test/test_ens_generate_pseudo_ensemble.nc ../data/ens.nc
ncatted -O -a coordinates,temp,o,c,"lat_b lon_b" ../data/ens.nc
