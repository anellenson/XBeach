function xb = xb_read_bcflist(filename, varargin)
%XB_READ_BCFLIST  Reads bcflist files generated by XBeach
%
%   Reads ebcflist or qbcflist files generated by XBeach. The files contain
%   references to other files containing realized wave and flux fields. The
%   referred files are read as well. The result is retruned in the form of
%   a XBeach structure.
%
%   Syntax:
%   xb = xb_read_bcflist(filename, varargin)
%
%   Input:
%   filename    = filename of the ebcflist or qbcflist file to be read
%   varargin    = range:    unity-based numerical range of files to be read
%                           (e.g. [4 5], 4 or [1 10])
%
%   Output:
%   xb          = XBeach structure array
%
%   Example
%   xb = xb_read_bcflist(filename)
%
%   See also xb_read_waves

%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2010 Deltares
%       Bas Hoonhout
%
%       bas.hoonhout@deltares.nl	
%
%       Rotterdamseweg 185
%       2629HD Delft
%
%   This library is free software: you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation, either
%   version 2.1 of the License, or (at your option) any later version.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library. If not, see <http://www.gnu.org/licenses/>.
%   --------------------------------------------------------------------

% This tool is part of <a href="http://OpenEarth.nl">OpenEarthTools</a>.
% OpenEarthTools is an online collaboration to share and manage data and 
% programming tools in an open source, version controlled environment.
% Sign up to recieve regular updates of this function, and to contribute 
% your own tools.

%% Version <http://svnbook.red-bean.com/en/1.5/svn.advanced.props.special.keywords.html>
% Created: 24 Nov 2010
% Created with Matlab version: 7.9.0.529 (R2009b)

% $Id: xb_read_bcflist.m 6208 2012-05-15 15:30:24Z hoonhout $
% $Date: 2012-05-15 17:30:24 +0200 (Tue, 15 May 2012) $
% $Author: hoonhout $
% $Revision: 6208 $
% $HeadURL: https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/applications/xbeach/xb_io/xb_read_bcflist.m $
% $Keywords: $

%% read options

OPT = struct( ...
    'range', [] ...
);

OPT = setproperty(OPT, varargin{:});

if length(OPT.range) == 1
    OPT.range = OPT.range*[1 1];
end

%% read bcf list

if ~exist(filename, 'file')
    error(['File does not exist [' filename ']'])
end

fdir = fileparts(filename);

[bcendtime rt dt trep mainang data] = deal([]);
filenames = {};

data = [];
tlength = 1;
fid = fopen(filename);
while ~feof(fid)
    fline = fgetl(fid);
    if isempty(fline); continue; end;

    [bcendtime(tlength) rt(tlength) dt(tlength) trep(tlength) mainang(tlength) fname] = ...
        strread(fline, '%f%f%f%f%f%s\n', 'delimiter', ' ');

    fname = fullfile(fdir, [fname{:}]);
    filenames = [filenames {fname}];
    
    if isempty(OPT.range) || (tlength >= OPT.range(1) && tlength <= OPT.range(2))
        datai = xb_read_bcffile(fname);
        data = cat(ndims(datai), data, datai);
    end
    
    tlength = tlength+1;
end
fclose(fid);

%% create xbeach struct

xb = xs_empty();
xb = xs_set(xb, '-units', 'time', {bcendtime 's'}, 'duration', {rt 's'}, ...
    'timestep', {dt 's'}, 'Trep', {trep 's'}, 'main_angle', {mainang 'degrees'}, ...
    'data', {data 'J/m^2'});
xb = xs_consolidate(xb);
xb = xs_meta(xb, mfilename, 'boundaryconditions', [{filename} filenames]);