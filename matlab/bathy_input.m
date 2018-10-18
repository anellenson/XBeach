function bathy_input(matfilename) 


ncfilename = ['/home/server/pi/homes/aellenso/Research/XBeach/ArmyCorps/data/FRF_20170327_1132_FRF_NAVD88_LARC_GPS_UTC_v20170408.nc']
lat = ncread(ncfilename, 'lat');
lon = ncread(ncfilename, 'lon');
bathy = ncread(ncfilename,'elevation');
profile_nums = ncread(ncfilename,'profileNumber');

p_inds  = find(profile_nums == 960);
[x,y] = latLon2FRF(lat,lon);

z = bathy(p_inds);
x_transect = x(p_inds);
x_transect = x_transect-x_transect(1);
z = flipud(z);


[xgrd,zgrd] = xb_grid_xgrid(x_transect,z);

ygrd = zeros(size(xgrd));

xgrd_fname = ([basedir,'FRF.20150327.p960.x.grd'])
fid = fopen(xgrd_fname,'w')
fprintf(fid,'%6.6f  ',xgrd)
fclose(fid)

ygrd_fname = ([basedir,'FRF.20150327.p960.y.grd'])
fid = fopen(ygrd_fname,'w')
fprintf(fid,'%6.6f  ',ygrd)
fclose(fid)

bed_fname = ([basedir,'FRF.20150327.p960.bed.dep'])
fid = fopen(bed_fname,'w')
fprintf(fid,'%6.6f  ',zgrd)
fclose(fid)
