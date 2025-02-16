%OPENDAP_SUBSETTING_WITH_MATLAB_TUTORIAL how to benefit from OPeNDPA subsetting in Matlab
%
%See also: OPeNDAP_access_with_Matlab_tutorial 

% $Id: OPeNDAP_subsetting_with_Matlab_tutorial.m 5480 2011-11-14 17:46:41Z boer_g $
% $Date: 2011-11-14 18:46:41 +0100 (ma, 14 nov 2011) $
% $Author: boer_g $
% $Revision: 5480 $
% $HeadURL: https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/io/opendap/OPeNDAP_subsetting_with_Matlab_tutorial.m $
% $Keywords: $

% This document is also posted on a wiki: http://public.deltares.nl/display/OET/OPeNDAP+subsetting+with+matlab

%% Add snctools http://mexcdf.sourceforge.net/, shipped with OpenEarthTols
% run('...\openearthtools\matlab\oetsettings.m')

%% Define data on an opendap server
% for converting non-open access GECBO grids to same netCDF structure: see https://repos.deltares.nl/repos/OpenEarthRawData/trunk/gebco/
url_grid{1} = 'http://geoport.whoi.edu/thredds/dodsC/bathy/etopo1_bed_g2';
url_grid{2} = 'http://geoport.whoi.edu/thredds/dodsC/bathy/etopo2_v2c.nc';
url_grid{3} = 'http://geoport.whoi.edu/thredds/dodsC/bathy/srtm30plus_v1.nc';
url_grid{4} = 'http://geoport.whoi.edu/thredds/dodsC/bathy/srtm30plus_v6';
url_grid{5} = 'http://geoport.whoi.edu/thredds/dodsC/bathy/smith_sandwell_v9.1.nc';
url_grid{6} = 'http://geoport.whoi.edu/thredds/dodsC/bathy/smith_sandwell_v11';
url_grid{7} = 'F:\checkouts\OpenEarthRawData\gebco\raw\gebco_30sec.nc';
url_grid{8} = 'F:\checkouts\OpenEarthRawData\gebco\raw\gebco_1min.nc';

url_line    = 'http://opendap.deltares.nl/thredds/dodsC/opendap/noaa/gshhs/gshhs_i.nc';

%% Get line data: 1D vectors are small, so we can get all data
L.lon    = nc_varget(url_line,'lon');
L.lat    = nc_varget(url_line,'lat');

%% Define bounding box
boundingbox.lon = [ 0 10];
boundingbox.lat = [50 55];

for i=1:length(url_grid)

   ncfile = url_grid{i}

   %% Get full lat,lon vectors: 1D vectors are small, so we can get all data
   nc_dump(ncfile)
   G.lon    = nc_varget(ncfile,'lon' ); % 1D
   G.lat    = nc_varget(ncfile,'lat' ); % 1D
   
   %% Find the subset-indices within the bounding box
   ilon     = find(G.lon > boundingbox.lon(1) & G.lon < boundingbox.lon(2));
   ilat     = find(G.lat > boundingbox.lat(1) & G.lat < boundingbox.lat(2));
   
   %% Translate subset-indices to netCDF argument: [start,count,stride]
   stride   = [1 1]; % additionally specify a stride when the subset is still too big
   start    = [min(ilat)-1 min(ilon)-1]; % subtract one as netCDF is 0-based, whereas matlab is 1-bases
   count    = ceil([length(ilat) length(ilon)]./stride); % use ceil to cover at least bounding box area
   
   %% Request data subset
   G.lat    = nc_varget(ncfile,'lat' ,start(1),count(1),stride(1)); % 1D
   G.lon    = nc_varget(ncfile,'lon' ,start(2),count(2),stride(2)); % 1D
   G.topo   = nc_varget(ncfile,'topo',start(:),count(:),stride(:)); % 2D
   G.title  = nc_attget(ncfile,nc_global,'title');
   
   %% Plot data subset
   figure(i)
   h = pcolorcorcen(G.lon,G.lat,double(G.topo))
   hold on
   plot(L.lon,L.lat,'k')
   axis([boundingbox.lon boundingbox.lat])
   tickmap('ll','texttype','text','dellast',1)
   axislat % sets aspect ratio
   grid on
   clim ([-50 150])
   title(mktex(G.title))
   print2screensize(mkvar(G.title))

   %% Plot data source subset (GEBCO only, from version created by OpenEarthRawData conversion script)
   if nc_isvar(ncfile,'SID')
   G.topo   = nc_varget(ncfile,'SID' ,start(:),count(:),stride(:)); % 2D
   delete(h)
   h = pcolorcorcen(G.lon,G.lat,double(G.topo));
   ctick = nc_attget(ncfile,'SID','flag_values');
   colormap(colormap_cpt('Paired 12',length(ctick)))
   clim ([min(ctick) max(ctick)]+[-.5 .5])
   [ax,c]=colorbarwithvtext('source',ctick);
   set(ax,'yticklabel',strtokens2cell(nc_attget(ncfile,'SID','flag_short_labels'),';'))
   set(ax,'FontSize',5)
   print2screensize([mkvar(G.title),'_SID'])
   end
   
end   