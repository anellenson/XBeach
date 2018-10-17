function OK = isurl(str)
%ISURL   boolend whether char is url or not
%
% ok = isurl(string)
%
% See also: urlread

% $Id: isurl.m 7803 2012-12-07 11:19:30Z boer_g $
% $Date: 2012-12-07 12:19:30 +0100 (Fri, 07 Dec 2012) $
% $Author: boer_g $
% $Revision: 7803 $
% $HeadURL: https://svn.oss.deltares.nl/repos/openearthtools/trunk/matlab/general/string_fun/isurl.m $
% $Keywords$

OK = (length(str) >5 && strcmpi(str(1:6),'ftp://'  )) || ...
     (length(str) >6 && strcmpi(str(1:7),'http://' )) || ...
     (length(str) >7 && strcmpi(str(1:8),'https://'));

%% EOF
