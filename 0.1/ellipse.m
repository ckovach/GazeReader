function pl = ellipse(cent,trans,fig,varargin)


%Plots an ellipse centered at cent and with axes tranformed by matrix
%trans.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%Kovach 2008 


nplotpoint = 100;

if nargin < 2 | isempty(trans)
    
    trans = eye(2);
end

if nargin < 3 | isempty(fig)
    
    fig = figure;
else
    figure(fig)
end
% if nargin < 4
%     plotstring = {};
% elseif ~iscell(plotstring)
%     plotstring = {plotstring};
% end

ang = [0:nplotpoint]./nplotpoint*2*pi;

circ = cat(1,sin(ang),cos(ang));
cent = cent(:);
ellipse = trans*circ + cent(:,ones(1,nplotpoint+1));

% arm1 = trans*[-1 0;1 0]'; 
% arm2 = trans*[0 -1 ;0 1 ]'; 
[arms,len] = svd(trans); 


pl = plot(ellipse(1,:),ellipse(2,:),varargin{:});
hold on,
pl(2) = plot([-arms(1,1),arms(1,1)]*len(1) + cent(1),[-arms(2,1),arms(2,1)]*len(1)+ cent(2),varargin{:});
pl(3) = plot([-arms(1,2),arms(1,2)]*len(2,2) + cent(1),[-arms(2,2),arms(2,2)]*len(2,2)+ cent(2),varargin{:});

