function dprintf(debug,varargin)
%DPRINTF   fprintf with option to write to screen
%
%   dprintf(fid,text)
%
% writes text to fid, and does nothing when fid=0.
% dprintf exists because fid=0 used to work with 
% FPRITNF in earlier Matlab versions.
%
% Useful for debug purposes: messages can be written to a log file or screen.
%
% Example:
%
%   dprintf( 0,'read such and so succesfully\n') % does nothing
%   dprintf( 1,'read such and so succesfully\n') % writes to screen
%   dprintf(>2,'read such and so succesfully\n') % write to file with identiefier
%
%See also: FPRINTF, WARNING, DISP

   if debug > 0
       fprintf(debug,varargin{:});
   end   

  % if debug
  %     fprintf(debug,varargin{:});
  %     fprinteol(debug);
  % else
  %     disp(varargin{:});
  % end
   

%% EOF
