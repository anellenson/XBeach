function varargout=wlgrid(cmd,varargin)
%WLGRID Read/write a Delft3D grid file.
%   GRID = WLGRID('read',FILENAME) reads a Delft3D or SWAN grid file and
%   returns a structure containing all grid information.
%
%   [X,Y,ENC,CS,nodatavalue] = WLGRID('read',FILENAME) reads a Delft3D
%   or SWAN grid file and returns X coordinates, Y coordinates, grid
%   enclosure and coordinate system string and nodatavalue separately.
%
%   OK = WLGRID('write','PropName1',PropVal1,'PropName2',PropVal2, ...)
%   writes a grid file. The following property names are accepted when
%   following the write command
%     'FileName'        : Name of file to write
%     'X'               : X-coordinates (can be a 1D vector 'stick' of an
%                         orthogonal grid)
%     'Y'               : Y-coordinates (can be a 1D vector 'stick' of an
%                         orthogonal grid)
%     'Enclosure'       : Enclosure array
%     'CoordinateSystem': Coordinate system 'Cartesian' (default) or
%                         'Spherical'
%     'Format'          : 'NewRGF'   - keyword based, double precision file
%                                      (default)
%                         'OldRGF'   - single precision file
%                         'SWANgrid' - SWAN formatted single precision file
%     'MissingValue'    : Missing value to be used for NaN coordinates
%
%   Accepted without property name: x-coordinates, y-coordinates and
%   enclosure array (in this order), file name, file format, coordinate
%   system strings 'Cartesian' and 'Spherical' (non-abbreviated).
%
%   See also ENCLOSURE, WLDEP.

%----- LGPL --------------------------------------------------------------------
%
%   Copyright (C) 2011-2012 Stichting Deltares.
%
%   This library is free software; you can redistribute it and/or
%   modify it under the terms of the GNU Lesser General Public
%   License as published by the Free Software Foundation version 2.1.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%   Lesser General Public License for more details.
%
%   You should have received a copy of the GNU Lesser General Public
%   License along with this library; if not, see <http://www.gnu.org/licenses/>.
%
%   contact: delft3d.support@deltares.nl
%   Stichting Deltares
%   P.O. Box 177
%   2600 MH Delft, The Netherlands
%
%   All indications and logos of, and references to, "Delft3D" and "Deltares"
%   are registered trademarks of Stichting Deltares, and remain the property of
%   Stichting Deltares. All rights reserved.
%
%-------------------------------------------------------------------------------
%   http://www.deltaressystems.com
%   $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/trunk/src/tools_lgpl/matlab/quickplot/progsrc/wlgrid.m $
%   $Id: wlgrid.m 2105 2013-01-13 22:34:24Z jagers $

if nargin==0
    if nargout>0
        varargout=cell(1,nargout);
    end
    return
end

switch lower(cmd)
    case {'read','open'}
        Grid=Local_read_grid(varargin{:});
        if nargout<=1
            varargout{1}=Grid;
        else
            varargout={Grid.X Grid.Y Grid.Enclosure Grid.CoordinateSystem Grid.MissingValue};
        end
    case {'write'}
        Out=Local_write_grid('newrgf',varargin{:});
        if nargout>0
            varargout{1}=Out;
        end
    case {'writeold'}
        Out=Local_write_grid('oldrgf',varargin{:});
        if nargout>0
            varargout{1}=Out;
        end
    case {'writeswan'}
        Out=Local_write_grid('swangrid',varargin{:});
        if nargout>0
            varargout{1}=Out;
        end
    otherwise
        error('Unknown command')
end


function GRID=Local_read_grid(filename)
GRID.X                = [];
GRID.Y                = [];
GRID.Enclosure        = [];
GRID.FileName         = '';
GRID.CoordinateSystem = 'Unknown';
GRID.MissingValue     = 0;

if (nargin==0) || strcmp(filename,'?')
    [fname,fpath]=uigetfile('*.*','Select grid file');
    if ~ischar(fname)
        return
    end
    filename=fullfile(fpath,fname);
end

% detect extension
[path,name,ext]=fileparts(filename);
if isempty(ext)
    if ~exist(filename)
        % add default .grd extension
        ext='.grd';
    end
end
filename=fullfile(path,[name ext]);
basename=fullfile(path,name);
GRID.FileName=filename;

% Grid file
gridtype='RGF';
fid=fopen(filename);
if fid<0
    error('Couldn''t open file: %s.',filename)
end
%
try
    prevlineoffset = 0;
    line=fgetl(fid);
    if ~isempty(strfind(lower(line),'spherical'))
        GRID.CoordinateSystem = 'Spherical';
    end
    while 1
        values = sscanf(line,'%f');
        if ~ischar(line)
            fclose(fid);
            error('The file is empty: %s.',filename)
        elseif ~isempty(line) && line(1)=='*'
            prevlineoffset = ftell(fid);
            line=fgetl(fid);
        elseif strcmp('x coordinates',deblank(line)) || strcmp('x-coordinates',deblank(line))
            gridtype='SWAN';
            GRID.CoordinateSystem='Cartesian';
            break
        elseif prevlineoffset == 0 && length(values)>3
            gridtype='ECOMSED-corners';
            break
        else
            EqualSign = strfind(line,'=');
            if ~isempty(EqualSign)
                keyword = deblank2(line(1:EqualSign(1)-1));
                switch keyword
                    case 'Coordinate System'
                        GRID.CoordinateSystem=deblank2(line(EqualSign(1)+1:end));
                    case 'Missing Value'
                        GRID.MissingValue=str2double(line(EqualSign(1)+1:end));
                end
                prevlineoffset = ftell(fid);
                line=fgetl(fid);
            else
                % Not a line containing a keyword. This can happen if it is an old
                % grid file with the first line starting with a *, in which case
                % this line should not yet have been read. Or it can be an old
                % file,where the first line does not start with a *, in which case
                % this line should be skipped anyway. To distinguish between these
                % two cases, check prevlineoffset: if it is zero, then this is the
                % first line and it should be skipped, otherwise unskip this line.
                if prevlineoffset>0
                    fseek(fid,prevlineoffset,-1);
                end
                break
            end
        end
    end
    %
    grdsize=[];
    switch gridtype
        case 'RGF'
            while 1
                line=fgetl(fid);
                if ~ischar(line)
                    break
                end
                if isempty(deblank(line)) && isempty(grdsize)
                elseif line(1)=='*'
                    % skip comment
                else
                    grdsize_f=transpose(sscanf(line,'%f'));
                    grdsize=transpose(sscanf(line,'%i'));
                    %
                    % don't continue if there are floating point numbers on the
                    % grid size line.
                    %
                    if length(grdsize_f) > length(grdsize)
                        error('Floating point number in line that should contain grid dimensions\n%s',line)
                    end
                    %
                    line = fgetl(fid); % read xori,yori,alfori
                    pars = sscanf(line,'%f');
                    if length(pars)>3
                        error('Too many parameters on origin line: %s',line)
                    end
                    if length(grdsize)>2 % the possible third element contains the number of subgrids
                        if grdsize(3)>1
                            for i=1:grdsize(3) % read subgrid definitions
                                fgetl(fid); % Npart =       1 
                                fgetl(fid); %        1      54     415     578       1       0       0     414 
                            end
                        end
                    end
                    if isempty(grdsize)
                        error('Cannot determine grid size.')
                    end
                    grdsize=grdsize(1:2);
                    floc=ftell(fid);
                    str=fscanf(fid,'%11c',1);
                    fseek(fid,floc,-1);
                    cc = sscanf(str,' %*[Ee]%*[Tt]%*[Aa] = %i',1);
                    readETA=0;
                    if ~isempty(cc)
                        readETA=1;
                    end
                    GRID.X=-999*ones(grdsize);
                    for c=1:grdsize(2)
                        if readETA
                            fscanf(fid,' %*[Ee]%*[Tt]%*[Aa] = %i',1); % skip line header ETA= and read c
                        else
                            fscanf(fid,'%11c',1); % this does not include the preceding EOL
                        end
                        GRID.X(:,c)=fscanf(fid,'%f',[grdsize(1) 1]);
                        if ~readETA
                            fgetl(fid); % read optional spaces and EOL
                        end
                    end
                    GRID.Y=-999*ones(grdsize);
                    for c=1:grdsize(2)
                        if readETA
                            fscanf(fid,' %*[Ee]%*[Tt]%*[Aa] = %i',1); % skip line header ETA= and read c
                        else
                            fscanf(fid,'%11c',1); % this does not include the preceding EOL
                        end
                        GRID.Y(:,c)=fscanf(fid,'%f',[grdsize(1) 1]);
                        if ~readETA
                            fgetl(fid); % read optional spaces and EOL
                        end
                    end
                    break
                end
            end
        case 'SWAN'
            %SWANgrid
            xCoords={};
            firstline=1;
            NCoordPerLine=0;
            while 1
                line=fgetl(fid);
                if firstline
                    [Crd,cnt]=sscanf(line,'%f',[1 inf]);
                    NCoordPerLine=cnt;
                    firstline=0;
                else
                    [Crd,cnt]=sscanf(line,'%f',[1 NCoordPerLine]);
                end
                if cnt>0
                    xCoords{end+1}=Crd;
                end
                if cnt<NCoordPerLine
                    break
                end
            end
            NCnt=-1;
            if cnt>0
                NCnt=(length(xCoords)-1)*NCoordPerLine+length(xCoords{end});
                while 1
                    line=fgetl(fid);
                    [Crd,cnt]=sscanf(line,'%f',[1 NCoordPerLine]);
                    if cnt>0
                        xCoords{end+1}=Crd;
                    else
                        break
                    end
                end
            end
            if ~strcmp('y coordinates',deblank(line)) && ~strcmp('y-coordinates',deblank(line))
                error('y coordinates string expected.')
            end
            xCoords=cat(2,xCoords{:});
            yCoords=zeros(size(xCoords));
            offset=0;
            while 1
                line=fgetl(fid);
                if ~ischar(line)
                    break;
                else
                    [Crd,cnt]=sscanf(line,'%f',[1 NCoordPerLine]);
                    if cnt>0
                        yCoords(offset+(1:cnt))=Crd;
                        offset=offset+cnt;
                    else
                        break
                    end
                end
            end
            if NCnt<0
                %
                % Determine grid dimensions by finding the first large change in the
                % location of successive points. Does not work if the grid contains
                % "missing points" such as the FLOW grid.
                %
                dist=sqrt(diff(xCoords).^2+diff(yCoords).^2);
                NCnt=min(find(dist>max(dist)*0.9));
            end
            xCoords=reshape(xCoords,[NCnt length(xCoords)/NCnt]);
            yCoords=reshape(yCoords,[NCnt length(yCoords)/NCnt]);
            GRID.X=xCoords;
            GRID.Y=yCoords;
        case 'ECOMSED-corners'
            fseek(fid,0,-1);
            GRID = read_ecom_corners(fid);
            GRID.CoordinateSystem='Unknown';
            GRID.MissingValue=0;
        otherwise
            error('Unknown grid type: %s',gridtype)
    end
    fclose(fid);
catch
    fclose(fid);
    rethrow(lasterror)
end
if isempty(GRID.X)
    error('File does not match Delft3D grid file format.')
end
%
% XX=complex(GRID.X,GRID.Y);
% [YY,i1,i2]=unique(XX(:));
% ZZ=sparse(i2,1,1);
% [YYm,n]=max(ZZ);
% if YYm>1
%    Missing = XX==YY(n);
%    GRID.X(Missing)=NaN;
%    GRID.Y(Missing)=NaN;
% end
%
notdef=(GRID.X==GRID.MissingValue) & (GRID.Y==GRID.MissingValue);
GRID.X(notdef)=NaN;
GRID.Y(notdef)=NaN;
GRID.Type = gridtype;

% Grid enclosure file
fid=fopen([basename '.enc']);
if fid>0
    Enc = [];
    while 1
        line=fgetl(fid);
        if ~ischar(line)
            break
        end
        Enc=[Enc; sscanf(line,'%i',[1 2])];
    end
    GRID.Enclosure=Enc;
    fclose(fid);
else
    %warning('Grid enclosure not found.');
    [M,N]=size(GRID.X);
    GRID.Enclosure=[1 1; M+1 1; M+1 N+1; 1 N+1; 1 1];
end


function OK = Local_write_grid(varargin)
% GRDWRITE writes a grid file
%       GRDWRITE(FILENAME,GRID) writes the GRID to files that
%       can be used by Delft3D and TRISULA.

fileformat = 'newrgf';
%
i          = 1;
filename   = '';
nparset    = 0;
Grd.CoordinateSystem='Cartesian';
while i<=nargin
    if ischar(varargin{i})
        switch lower(varargin{i})
            case 'oldrgf'
                fileformat      ='oldrgf';
            case 'newrgf'
                fileformat      ='newrgf';
            case 'swangrid'
                fileformat      ='swangrid';
            case 'cartesian'
                Grd.CoordinateSystem='Cartesian';
            case 'spherical'
                Grd.CoordinateSystem='Spherical';
            otherwise
                Cmds = {'CoordinateSystem','MissingValue','Format', ...
                    'X','Y','Enclosure','FileName'};
                j = ustrcmpi(varargin{i},Cmds);
                if j>0
                    i=i+1;
                    switch Cmds{j}
                        case 'CoordinateSystem'
                            Cmds = {'Cartesian','Spherical'};
                            j = ustrcmpi(varargin{i},Cmds);
                            if j<0
                                error('Unknown coordinate system: %s',varargin{i})
                            else
                                Grd.CoordinateSystem=Cmds{j};
                            end
                        case 'MissingValue'
                            Grd.MissingValue = varargin{i};
                        case 'Format'
                            fileformat    = lower(varargin{i});
                        case 'X'
                            Grd.X         = varargin{i};
                        case 'Y'
                            Grd.Y         = varargin{i};
                        case 'Enclosure'
                            Grd.Enclosure = varargin{i};
                        case 'FileName'
                            filename      = varargin{i};
                    end
                elseif isempty(filename) && ~isempty(varargin{i})
                    filename=varargin{i};
                else
                    error('Cannot interpret: %s',varargin{i})
                end
        end
    elseif isnumeric(varargin{i})
        switch nparset
            case 0
                Grd.X         = varargin{i};
            case 1
                Grd.Y         = varargin{i};
            case 2
                Grd.Enclosure = varargin{i};
            case 3
                error('Unexpected numerical argument.')
        end
        nparset = nparset+1;
    else
        Grd = varargin{i};
        if ~isfield(Grd,'CoordinateSystem')
            Grd.CoordinateSystem='Cartesian';
            if isfield(Grd,'CoordSys')
                Grd.CoordinateSystem = Grd.CoordSys;
                Grd = rmfield(Grd,'CoordSys');
            end
        end
    end
    i=i+1;
end

% allow to supply orthogonal grid defined with two 1D sticks
%  while making sure the resulting grid is RIGHT-HANDED
%  (when m is locally defined as x, n should locally point in positive y direction)
sz.x = size(Grd.X);
sz.y = size(Grd.Y);
if ~isequal(sz.x,sz.y)
    if min(sz.x(1:2))==1 && min(sz.y(1:2))==1
        [Grd.X,Grd.Y] = ndgrid(Grd.X(:),Grd.Y(:)); % NOT meshgrid, which leads not to righthanded coordinate system.
    else
        error(['X and Y should have same size:',num2str(sz.x),' vs. ',num2str(sz.y)])
    end
end

% enclosure
if ~isfield(Grd,'Enclosure')
    Grd.Enclosure=[];
end
if ~isfield(Grd,'MissingValue')
    Grd.MissingValue=0;
end

%
j = ustrcmpi(fileformat,{'oldrgf','newrgf','swangrid'});
switch j
    case 1
        Format='%11.3f';
    case 2
        Format='%25.17e';
    case 3
        Format='%15.8f';
    otherwise
        Format=fileformat;
        fileformat='newrgf';
end
if isfield(Grd,'Format')
    Format=Grd.Format;
end

%
if Format(end)=='f'
    ndec       = sscanf(Format,'%%%*i.%i');
    coord_eps  = 10^(-ndec);
else
    coord_eps  = eps*max(1,abs(Grd.MissingValue));
end

% move any points that have x-coordinate/y-coordinate equal to missing
% value coordinate
Grd.X(Grd.X==Grd.MissingValue)=Grd.MissingValue+coord_eps;
Grd.Y(Grd.Y==Grd.MissingValue)=Grd.MissingValue+coord_eps;
Idx=isnan(Grd.X.*Grd.Y);                % change indicator of grid point exclusion
Grd.X(Idx)=Grd.MissingValue;            % from NaN in Matlab to (0,0) in grid file.
Grd.Y(Idx)=Grd.MissingValue;

% detect extension
[path,name,ext]=fileparts(filename);
if isempty(ext)
    ext='.grd';
end % default .grd file
filename=fullfile(path,[name ext]);
basename=fullfile(path,name);

if ~isempty(Grd.Enclosure),
    fid=fopen([basename '.enc'],'w');
    if fid<0
        error('* Could not open output file.')
    end
    fprintf(fid,'%5i%5i\n',Grd.Enclosure');
    fclose(fid);
end

% write
fid=fopen(filename,'w');
if strcmp(fileformat,'oldrgf') || strcmp(fileformat,'newrgf')
    if strcmp(fileformat,'oldrgf')
        fprintf(fid,'* MATLAB Version %s file created at %s.\n',version,datestr(now));
    elseif strcmp(fileformat,'newrgf')
        fprintf(fid,'Coordinate System = %s\n',Grd.CoordinateSystem);
        if isfield(Grd,'MissingValue') && Grd.MissingValue~=0
            fprintf(fid,['Missing Value = ' Format '\n'],Grd.MissingValue);
        end
    end
    fprintf(fid,'%8i%8i\n',size(Grd.X));
    fprintf(fid,'%8i%8i%8i\n',0,0,0);
    Format=cat(2,' ',Format);
    M=size(Grd.X,1);
    for fld={'X','Y'}
        fldx=fld{1};
        coords = Grd.(fldx);
        for j=1:size(coords,2)
            fprintf(fid,' ETA=%5i',j);
            fprintf(fid,Format,coords(1:min(5,M),j));
            if M>5
                fprintf(fid,cat(2,'\n          ',repmat(Format,1,5)),coords(6:M,j));
            end
            fprintf(fid,'\n');
        end
    end
elseif strcmp(fileformat,'swangrid')
    for fld={'X','Y'}
        fldx  =fld{1};
        fprintf(fid,'%s%s\n',lower(fldx),'-coordinates');
        coords = Grd.(fldx);
        Frmt  =repmat('%15.8f  ',[1 size(coords,1)]);
        k=8*3;
        Frmt((k-1):k:length(Frmt))='\';
        Frmt( k   :k:length(Frmt))='n';
        Frmt(end-1:end)='\n';
        fprintf(fid,Frmt,coords);
    end
end
%---------
fclose(fid);
OK=1;
