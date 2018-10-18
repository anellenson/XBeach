function bathy_input(ncfilename) 
cd('xbeach_release_10May2013')
basedir = '../../model/grids/'

load([basedir, '3.2017.simulated_transect.mat'])
z = flipud(z);

addpathfast('.')

[xgrd,zgrd] = xb_grid_xgrid(x,z);

ygrd = zeros(size(xgrd));

xgrd_fname = ([basedir,'p960.201703.x.grd'])
fid = fopen(xgrd_fname,'w')
fprintf(fid,'%6.6f  ',xgrd)
fclose(fid)

ygrd_fname = ([basedir,'p960.201703.y.grd'])
fid = fopen(ygrd_fname,'w')
fprintf(fid,'%6.6f  ',ygrd)
fclose(fid)

bed_fname = ([basedir,'p960.201703.bed.dep'])
fid = fopen(bed_fname,'w')
fprintf(fid,'%6.6f  ',zgrd)
fclose(fid)
