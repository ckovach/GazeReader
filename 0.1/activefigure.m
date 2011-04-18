function varargout = activefigure(h)

%Makes the specified figure active without stealing focus.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin < 1
    h = figure(h);
    
elseif ~ishandle(h)
    
    h = figure(h);
    
else
    
    set(0,'CurrentFigure',h);
end

 if nargout > 0
     varargout{1} = h;
 end