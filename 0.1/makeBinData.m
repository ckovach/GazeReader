
function [binData, binfunctions] = makeBinData(varargin)

% binData = makeBinData(...)
%
%   Creates a sampling bin data structure with the following
%   fields:
% 
%     binData.groups - a structure array where each element definines an
%     bin or groups of bins.
% 
%     binData.precedence - a field indicating how membership in overlapping
%                   bins is to be resolved. This may be
%                        'nearest_to_center':  bin with the center nearest to the fixation point wins
%                        'array_order':   The bin specified earliest in the groups array wins.
% 
% 	binData = makeBinData(data,'type',bintype) creates a new set of bins
%       based on the information in data. The content of data differs for
%       differnt types:
%     
%         bintype = 'grid' creates a grid of rectangular bins
% 
%             data is a 2 element cell array with the first element a 4 element vector defining the 
%             dimensions of the grid as [xmin xmax ymin ymax] and data{2} a 2 element vector giving 
%             the number of divisions in each respective dimension. 
% 
%         bintype = 'rect' creates a set of rectangular bins
% 
%             data is a Nx4 matrix of rectangles.
%       
% 
%         bintype = 'ellipse' creates a set of elliptical bins
% 
%             If data is a Nx3 matrix of centers and radii,then bins are
%                circular within coordinates normalized to screen dimensions
%                (hence, oblate in typical screen coordinates). If data is a 
%                 Nx4 matrix of centers and distance for the x and y
%                axes, then data is a Nx4 matrix of centers and distance for the x and y
%                axes. If data is Nx5, then the 5th dimension gives the angle of rotation 
%                of one semi-axis counterclockwise with respect to the x axis.
% 
%         bintype = 'poly' creates a single polygonal bin
% 
%             data is a Mx2 matrix of vertices.
%
%         bintype = 'simplex' creates N simplicial bins 
%             data{1} is an Mx2 matrix of vertices.
%             data{2} is an Nx3 matrix of indices into data{1}.
% 
% Additional options:
% 
%       makeBinData(...,'label',binlabel) assigns a label to the group of
%       bins.
% 
%       makeBinData(...,'trials',trialvec) specifies which trials a bin is
%       active,  or the number of sequential trials if trialvec is
%       a scalar.
% 
%       makeBinData(binData,data,...) appends data to the binData
%       structure, binData.
% 
%       makeBinData(binData1,binData2,...) concatenates two binData structures.
%
%       binData.groups is a structure array with the following fields:
%           .type       -   type of bin (grid, simplex, rect, etc.)
%           .label      -   descriptive label
%           .binnums    -   bin numbers
%           .pos        -   array or structure defining bin position (depending on type)
%           .nbin       -   number of bins in group
%           .centers    -   bin centers [x y]
%           .isinside   -   function handle which takes xy coordinate as
%                          input and returns bin membership indicator
%           .plot       -  function for plotting bin boundaries
%           .patch      - function for plotting bin as a patch
%           .activeTrials - unused
%           .inputData  -  data argument passed to makeBinData
%           .code       - unique identifying code

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


if nargin > 0 && (~isstruct(varargin{1}) || isfield(varargin{1},'vert'))
    data = varargin{1};
    i = 2;   
    append = 0;
    concat= 0;
elseif nargin > 1 &&~isstruct(varargin{2})   
    binData = varargin{1};
    data = varargin{2};
    i = 3;
    append = 1;
    concat= 0;
elseif nargin > 1
    binData = varargin{1};
    binData2 = varargin{2};
%     data = varargin{3};
    i = 4;
    append = 1;
    concat= 1;
end

if concat

    newcode = [binData2.groups.code] + binData.codeincr;
    c = mat2cell(newcode,1,ones(1,length(newcode)));
    [binData2.groups.code] = c{:};    
    binData.groups = [binData.groups,binData2.groups];
    binData.codeincr= max([binData.groups.code]) ;
    
    return
    
end

type = 'undefined';
trials = [];
code = [];
ntrials = 0;
% units = 'normalized';
label = '';
% precedence = 'nearest_to_center';
precedence = 'group_order';
while i <= length(varargin)
   switch lower(varargin{i})
       case 'type'
          type = varargin{i+1};
          i = i+1;
       case 'trials'
          trials = varargin{i+1}; %1Xntrials logical vector with 1  whenever the to-be-assigned bin occurs
          i = i+1;
       case 'code'
          code = varargin{i+1};
          i = i+1;
       case 'gaps'  %codes of bingroups that are treated as gaps
          gaps = varargin{i+1};
          i = i+1;
       case 'ntrials'
          ntrials = varargin{i+1}; %Total number of trials
          i = i+1;
%        case 'screenunits'
%            units = 'screen';
       case 'label'  %Label for the bin group
           label  = varargin{i+1}; 
            i = i+1;
       case 'precedence'  %A string code used to decide how to resolve multiple assignment for overlapping bins.
           precedence  = varargin{i+1}; 
            i = i+1;
        otherwise

           error([varargin{i},' is not a valid option.']);
   end         
   i = i+1;
end


if ~append
        % binData = struct('groups',[],'activeTrials',[],'precedence','nearest_to_center');
    binData = struct('groups',BinGroupsTemplate,'precedence',precedence);
    binData.codeincr = 0;
    grpind = 1;
else
    grpind = length(binData.groups)+1;
end

binfunctions.isinside = @(varargin) isinside(varargin{:});
binfunctions.binarea = @(varargin) binarea(varargin{:});


binData.groups(grpind).label = label;
binData.groups(grpind).type= type;
binData.groups(grpind).nbin = 0;
if isempty(data)
    if isempty(varargin(2:end));
        binData.groups = binData.groups([]);
    else
        binData.groups(grpind).code = binData.codeincr+1;
        binData.codeincr = max([binData.groups(grpind).code]);
    end
    return
end


    % binData.groups.inputData = data;
if isempty(code)
    binData.groups(grpind).code = binData.codeincr+1;
else
    binData.groups(grpind).code = code;
end
% binData.codeincr = binData.codeincr+1;


if ~isempty(trials) && ntrials == 0;
    ntrials = max(trials);
end

switch lower(type)
    
    case 'grid'
        
        binData.groups(grpind) = makeGridbin(data);
        
    case 'poly'
        
        binData.groups(grpind) = makePolybin(data);
        
    case 'rect'
        
        binData.groups(grpind) = makeRectbin(data);

    case 'ellipse'                
        binData.groups(grpind) = makeEllipsebin(data);
 
    case {'tri','simplex'}                
        binData.groups(grpind) = makeTriBin(data);

    case 'circle'                 
        binData.groups(grpind) = makeEllipsebin(data);

    case 'undefined'
        if ~isempty(data)
            error('Type must be given for non-empty data.')
        end
    otherwise

       error([type,' is not a valid bin type.']);
end

% binData.groups(grpind).nbin = size(binData.groups(grpind).pos,1);

binData.groups(grpind).inputData = data;

binData.groups(grpind).label = label;
% multassign = @ (varargin) varargin{:};

if all(ismember(trials,[0 1]))
%     binData.groups(grpind).activeTrials = find(trials(:));
    binData.groups(grpind).activeTrials =  trials(:);
 
else 
%     binData.groups(grpind).activeTrials = trials;
        binData.groups(grpind).activeTrials = false( ntrials , 1 );
        binData.groups(grpind).activeTrials( trials ) = true;
end

if isempty(code)
    binData.groups(grpind).code = binData.codeincr + 1;
    binData.codeincr = binData.codeincr + 1;
else
    binData.groups(grpind).code = code;
    binData.codeincr = code;
end    
% if all(ismember(trials,[0 1]))
%     binData.activeTrials = trials(:);
% else 
%     binData.activeTrials = false( ntrials , 1 );
%     binData.activeTrials( trials ) = true;
% end




%%%%%%%%%%%%%%%%%%%%%%%
function BinGroups = BinGroupsTemplate
% initializes an empty BinGroups structure

BinGroups = struct('type',[],'label','','binnums',[],'pos',[],'nbin',[],'centers',[],'isinside',...
                    [],'plot',[],'patch',[],'area',[],'activeTrials',[],'inputData',[],'code',[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rect and Grid bins
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

function BinGroups = makeGridbin(data)


if ~iscell(data)    
    data = rect2grid(data); %assuming that if it's not a cell array it must be a rect matrix     
end    

rect = grid2rect(data);

BinGroups = makeRectbin(rect);
BinGroups.type = 'grid';
BinGroups.plot = @(BinGroups,varargin) PlotGridBin(BinGroups,varargin{:});
BinGroups.patch = @(BinGroups,varargin) DrawRectPatch(BinGroups,varargin{:});

BinGroups.nbin = size(BinGroups.pos,1);
%%%%%%%%%%%%%%%%%%%%%

function BinGroups = makeRectbin(rect)

%Create a set of rectangular bins

BinGroups = BinGroupsTemplate;
BinGroups.type = 'rect';

BinGroups.binnums = 1:size(rect,1); 

BinGroups.pos = rect;

BinGroups.centers = cat(2,sum(rect(:,1:2),2)./2,sum(rect(:,3:4),2)./2);


BinGroups.isinside = @(BinGroups,varargin) insideRect(BinGroups,varargin{:});
% BinGroups.isinside = @(BinGroups,varargin) isinside(BinGroups,varargin{:});
BinGroups.plot = @(BinGroups,varargin) PlotRectBin(BinGroups,varargin{:});
BinGroups.patch = @(BinGroups,varargin) DrawRectPatch(BinGroups,varargin{:});
BinGroups.area = @(BinGroups,varargin) areaRect(BinGroups,varargin{:});

BinGroups.nbin = size(BinGroups.pos,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = PlotRectBin(BinGroups,bins,varargin)

if ~ ( strcmp(BinGroups.type,'rect')) 
    error('This plot function is valid for bin''s of type rect only.')
end

if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

if isempty(bins) || bins == 0
    bins = 1:size(BinGroups.pos,1);
end

 
axlim = axis(ca);


ps = plot(axlim(2)*BinGroups.pos(bins,[1 1])', axlim(4)*BinGroups.pos(bins,[3 4])',varargin{:});

hold on,

ps = cat(2,ps,plot(axlim(2)*BinGroups.pos(bins,[2 2])', axlim(4)*BinGroups.pos(bins,[3 4])',varargin{:}));
ps = cat(2,ps,plot(axlim(2)*BinGroups.pos(bins,[1 2])', axlim(4)*BinGroups.pos(bins,[3 3])',varargin{:}));
ps = cat(2,ps,plot(axlim(2)*BinGroups.pos(bins,[1 2])', axlim(4)*BinGroups.pos(bins,[4 4])',varargin{:}));

axis(ca,axlim);

if nargout > 0
    varargout{1} = ps;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = DrawRectPatch(BinGroups,bins,Cdata,varargin)

% Draws a bin shaped patch
if ~ ( strcmp(BinGroups.type,'rect')) && ~ ( strcmp(BinGroups.type,'grid')) 
    error('This plot function is valid for bin''s of type rect and grid only.')
end

if nargin < 2 || isempty(bins) || max(bins) > size(BinGroups.pos,1)
    bins = 1:size(BinGroups.pos,1);
end

if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

axlim = axis(ca);

pos = BinGroups.pos;

%%%%
vx = cat(2,pos(bins,[1,2]),pos(bins ,[2,1]))';

vy = pos(bins ,[3,3,4,4])';

%%%%
% ph = patch(vx*screenres(1),vy*screenres(2),Cdata,varargin{:});
ph = patch(vx,vy,Cdata,varargin{:});


axis(ca,axlim);

if nargout > 0
    varargout{1} = ph;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X,Y] = rectVertices(BinGroups,bins,npts, varargin)
%Gets vertices of the ellipse bin

if ~ ( strcmp(BinGroups.type,'rect') )
    error('This plot function is valid for bin''s of type rect only.')
end

if nargin < 3
    npts = 500;
end

nu = linspace(0,1,npts/4)';

nxy = zeros(4*length(nu),2);
nxy(1:length(nu),1) =  nu;
nxy(1:length(nu),2) =  0;
nxy(length(nu) + (1:length(nu)),1) =  1;
nxy(length(nu) + (1:length(nu)),2) =  nu;
nxy(2*length(nu) + (1:length(nu)),1) =  nu(end:-1:1);
nxy(2*length(nu) + (1:length(nu)),2) =  1;
nxy(3*length(nu) + (1:length(nu)),1) =  0;
nxy(3*length(nu) + (1:length(nu)),2) =  nu(end:-1:1);

X = nxy(:,1)*BinGroups.pos(:,3)' + nxy(:,1)*BinGroups.pos(:,1)'; 
Y = nxy(:,2)*BinGroups.pos(:,4)' + nxy(:,2)*BinGroups.pos(:,2)'; 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = PlotGridBin(BinGroups,bins,varargin)


%This plots the entire grid -- bins does not affect anything and is
%included in the argument list for compatibility.

if ~ ( strcmp(BinGroups.type,'grid') )
    error('This plot function is valid for bin''s of type grid only.')
end

if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

axlim = axis(ca);

x = BinGroups.pos(:,1);
dx = diff(x);
dx = dx(find(dx,1));
if isempty(dx)    
    dx = find(BinGroups.pos(:,4)- BinGroups.pos(:,3),1);
end

y = BinGroups.pos(:,3);
dy = diff(y);
dy = dy(find(dy,1));
if isempty(dy)
    dy = find(BinGroups.pos(:,4)-BinGroups.pos(:,3),1);
end

gridticky = y(1):dy:y(end)+dy;
gridtickx = x(1):dx:x(end)+dx; 

% varargout{1} = []
% return

hold(ca,'on')
out = plot( ca,repmat(gridtickx*axlim(2),2,1) , repmat(axlim(4)*gridticky([1 end])',1,length(gridtickx) ),varargin{:});
out= cat(1,out,plot( ca,repmat(axlim(2)*gridtickx([1 end])',1,length(gridticky) ), repmat(axlim(4)*gridticky,2,1),varargin{:}));
hold(ca,'off')

axis(ca,axlim)
varargout{1} = out;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binmember = insideRect(BinGroups,FixLoc,tol)

if nargin < 3
    tol = 0;
end

rect = BinGroups.pos;
nbin = size(BinGroups.pos,1);
nfix = size(FixLoc,1);

binmember =    sparse(repmat( FixLoc(:,1) , 1 , nbin ) >= repmat( rect(:,1)'-tol , nfix , 1 )  & ...
               repmat( FixLoc(:,1) , 1 , nbin ) <  repmat( rect(:,2)'+tol , nfix , 1 )  & ...
               repmat( FixLoc(:,2) , 1 , nbin ) >= repmat( rect(:,3)'-tol , nfix , 1 )  & ...
               repmat( FixLoc(:,2) , 1 , nbin ) <  repmat( rect(:,4)' +tol, nfix , 1 )); 

           
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binarea = areaRect(BinGroups)


edgelength = diff(BinGroups.pos,[],2);

binarea = prod(edgelength(:,[1 3]),2);


   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ellipse bin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BinGroups = makeEllipsebin(data)
 
% if isnumeric(data)
%     data= {data(:,1:2),data(:,3)};
% end
if iscell(data)
    data= [data{1},data{2}];
end

if ~isempty(data) && size(data,2) < 4 %if data has 3 columns, the y-semi axis will be the same as the x
    data(:,4) = data(:,3);
end

if ~isempty(data) && size(data,2) < 5 %rotation parameter set to 0 if unspecified
    data(:,5) = 0;
end

BinGroups = BinGroupsTemplate;

BinGroups.type = 'ellipse';
    


BinGroups.centers = data(:,1:2);
BinGroups.pos = data;
BinGroups.binnums = 1:size(data,1);
BinGroups.nbin = size(data,1);


BinGroups.isinside = @(BinGroups,FixLoc,varargin) insideEllipseBin(BinGroups,FixLoc,varargin{:});
% BinGroups.isinside = @(BinGroups,varargin) isinside(BinGroups,varargin{:});   
BinGroups.plot = @(BinGroups,varargin) plotEllipseBin(BinGroups,varargin{:});
BinGroups.patch= @(BinGroups,varargin) DrawEllipsePatch(BinGroups,varargin{:});
BinGroups.area = @(BinGroups,varargin) areaEllipse(BinGroups,varargin{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binmember = insideEllipseBin(BinGroups,FixLoc,tol)

if nargin < 3 || isempty(tol)
     tol = 0;
end
if ~isa(FixLoc,'double')
    FixLoc = double(FixLoc);
end
radii = BinGroups.pos(:,3:4) + tol;
rotangle = BinGroups.pos(:,5);
nbin = size(BinGroups.pos,1);
nfix = size(FixLoc,1);

centers = BinGroups.pos(:,1:2);

%inverse scaling*rotation matrices
invrotmat = zeros(2,2*nbin*nfix);
invrotmat(1,1:2:end) =  kron( cos(rotangle')./radii(:,1)',ones(1,nfix));
invrotmat(1,2:2:end) =  kron( sin(rotangle')./radii(:,2)',ones(1,nfix));
invrotmat(2,1:2:end) =  kron( -sin(rotangle')./radii(:,1)',ones(1,nfix));
invrotmat(2,2:2:end) =  kron( cos(rotangle')./radii(:,2)',ones(1,nfix));

spinvr = sparseblock(invrotmat,2);

%distances to centers
fixds = sparseblock( repmat(FixLoc,nbin,1) - kron(centers,ones(nfix,1)),1,'transpose');

%normalized distances to centers
normd = sum( (fixds*spinvr).^2,2);

%A point is inside if it's normalized distance is <= 1
 binmember =   reshape( normd <= 1, nfix, nbin);
 
% binmember =    sparse( ( repmat( FixLoc(:,1) , 1 , nbin ) - repmat( centers(:,1)' , nfix , 1 ) ).^2  + ...
%                        ( repmat( FixLoc(:,2) , 1 , nbin ) - repmat( centers(:,2)' , nfix , 1 ) ).^2 ... 
%                         <= repmat( (radii(:)'+tol).^2 , nfix , 1 )  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binarea = areaEllipse(BinGroups)

%area of eliptical bins
radii = BinGroups.pos(:,3:4);

binarea = pi*prod(radii,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xy = ellipseVertices(BinGroups,bins,npts, varargin)
%Gets vertices of the ellipse bin

if ~ ( strcmp(BinGroups.type,'ellipse') )
    error('This plot function is valid for bin''s of type ellipse only.')
end

if nargin < 3
    npts = 500;
end

th = linspace(0,2*pi,npts);

centers = BinGroups.pos(bins,1:2);

if size(BinGroups.pos,2) < 4
    radii = BinGroups.pos(bins,[3 3]) + eps;
else
    radii = BinGroups.pos(bins,3:4) + eps;
end

if size(BinGroups.pos,2) < 5
    rotangle = 0+eps;
else
    rotangle =   BinGroups.pos(bins,5)+eps;
end


nbin = length(bins);
%scaling*rotation matrices
nth = length(th);

%Transformation matrix from ellipse normalized coordinates to screen
%coordinates
trmats = zeros(3,3*nbin);
trmats(1,1:3:end) =  cos(rotangle').*radii(:,1)';
trmats(1,2:3:end) =  -sin(rotangle').*radii(:,1)';
trmats(2,1:3:end) =  sin(rotangle').*radii(:,2)';
trmats(2,2:3:end) =  cos(rotangle').*radii(:,2)';
trmats(3,1:3:end) =  centers(:,1)';
trmats(3,2:3:end) =  centers(:,2)';
trmats(3,3:3:end) =  1;

% spr = sparseblock(rotmat,2);
sptr = sparseblock(trmats,3); % Block diagonal sparse matrix of transformations for each ellipse

xy = unsparsify(sparseblock(repmat([cos(th(:)),sin(th(:)),ones(length(th),1)],nbin,1),nth,'transpose')*sptr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = plotEllipseBin(BinGroups,bins,varargin)

if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

if nargin < 2 || isempty(bins) || all(bins == 0)
    bins = 1:size(BinGroups.pos,1);
end

axlim = axis(ca);

xy = ellipseVertices(BinGroups,bins);

out = plot(xy(:,1:3:end), xy(:,2:3:end),varargin{:});

axis(ca,axlim);

if nargout > 0
    varargout{1} = out;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = DrawEllipsePatch(binGroup,bins,Cdata,varargin)

if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

if nargin < 2 || isempty(bins) || max(bins) > size(binGroup.pos,1)
    bins = 1:size(binGroup.pos,1);
end

axlim = axis(ca);

xy = ellipseVertices(binGroup,bins);

if length(Cdata) ==3  %If Cdata has length 3, patch wrongly assumes it's RGB
 Cdata(end+1) = nan;
 xy(:,end+(1:3)) = nan;
end
out = patch(xy(:,1:3:end), xy(:,2:3:end),Cdata,varargin{:});

axis(ca,axlim);

if nargout > 0
    varargout{1} = out;
end


                    
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polygonal bins - these can only represent 1 bin at a time: each row is a
% vertex
%%%%%%%%%%%%%%%%%%%%%%%%%%

function BinGroups = makePolybin(data)



BinGroups = BinGroupsTemplate;

BinGroups.type = 'poly';

BinGroups.pos = data;

BinGroups.centers = mean(data);

BinGroups.binnums = 1;
BinGroups.nbin = 1;

BinGroups.isinside = @(BinGroups,FixLoc,tol)insidePolyBin(BinGroups,FixLoc,tol);
% BinGroups.isinside = @(BinGroups,varargin) isinside(BinGroups,varargin{:});
BinGroups.plot = @(BinGroups,varargin) plotPolyBin(BinGroups,varargin{:});
BinGroups.patch= @(BinGroups,varargin) DrawPolyPatch(BinGroups,varargin{:});
BinGroups.area = @(BinGroups,varargin) areaPoly(BinGroups,varargin{:});
%%%%%%%%%%%%%

function binmember = insidePolyBin(BinGroups,FixLoc,tol)

%tol doesn't do anything for polygonal bins

binmember = inpolygon(FixLoc(:,1),FixLoc(:,2),BinGroups.pos(:,1),BinGroups.pos(:,2));



function binarea = areaPoly(BinGroups)

%area of the polygon

binarea = polyarea(BinGroups.pos(:,1),BinGroups.pos(:,2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = plotPolyBin(BinGroups,bins,varargin)


if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

axlim = axis(ca);

varargout = cell(1,nargout);

[varargout{:}] = plot(axlim(2)*BinGroups.pos(:,1),axlim(2)*BinGroups.pos(:,2),varargin{:});

axis(ca,axlim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = DrawPolyPatch(BinGroups,bins,Cdata,varargin)


if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

axlim = axis(ca);

varargout = cell(1,nargout);

% [varargout{:}] = plot(axlim(2)*BinGroups.pos(:,1),axlim(2)*BinGroups.pos(:,2),varargin{:});
[varargout{:}] = patch(BinGroups.pos(:,1),BinGroups.pos(:,2),Cdata,varargin{:});

axis(ca,axlim);


%%%%%%%%%%%%%%%%%%%%%%
%Simplicial bins: bins are represented as a matrix of vertices and indices
%representing simplexes (triangles).
%%%%%%%%%%%%%%%%%%%%%


function BinGroups = makeTriBin(data)
 
% if isnumeric(data)
%     data= {data(:,1:2),data(:,3)};
% end
if iscell(data)
    vert = data{2};
    tri = data{1};
elseif isstruct(data)
    vert = data.vert;
    tri = data.tri;
end

BinGroups = BinGroupsTemplate;

BinGroups.type = 'simplex';

% BinGroups.pos.trirep = TriRep(tri,vert);


% dv(:,end) = 1;
v2 = vert;
v2(:,end+1) = 1;
nd = size(vert,2);
invTr = zeros(size(tri,1)*(nd+1),nd+1);

centers = zeros(size(tri,1),nd);
for i= 1:size(tri,1);    
    invTr((1:nd+1) + (i-1)*(nd+1),(1:nd+1) ) =eye(nd+1)/v2(tri(i,:),:)'; 
    centers(i,:) = mean(vert(tri(i,:),:));
end

BinGroups.pos.vert =vert;
BinGroups.pos.tri =tri;
BinGroups.pos.invTr = invTr;  %Inverse projection for each bin (used to determine whether a point lies within the bin).

BinGroups.centers = centers;
% BinGroups.pos.vert =vert;
% BinGroups.pos.tri =tri;

BinGroups.binnums = 1:size(tri,1);

BinGroups.nbin = size(tri,1);

BinGroups.isinside = @(BinGroups,FixLoc,varargin) insideTriBin(BinGroups,FixLoc,varargin{:});
% BinGroups.isinside = @(BinGroups,varargin) isinside(BinGroups,varargin{:});   
BinGroups.plot = @(BinGroups,varargin) plotTriBin(BinGroups,varargin{:});
BinGroups.patch= @(BinGroups,varargin) DrawTriPatch(BinGroups,varargin{:});
BinGroups.area = @(BinGroups,varargin) areaTriBin(BinGroups,varargin{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function binmember = insideTriBin(BinGroups,FixLoc,tol)


nd = size(BinGroups.pos.vert,2);
nbin = size(BinGroups.pos.tri,1);
nfix = size(FixLoc,1);

FixLoc(:,end+1) = 1;


trfix = reshape(BinGroups.pos.invTr*FixLoc',nd+1,nbin,nfix);


%A point is inside if it falls within the unit simplex
 binmember =  sparse(permute(sum(trfix(1:end-1,:,:)>=0) == nd & sum(trfix(1:end-1,:,:) <= 1) == nd & trfix(end-1,:,:) <= 1-sum(trfix(1:end-2,:,:),1),[3 2 1]));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = plotTriBin(BinGroups,bins,varargin)


if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

axlim = axis(ca);

varargout = cell(1,nargout);

pos = BinGroups.pos;
[varargout{:}] = TriPlot(pos.tri,axlim(2)*pos.vert(:,1),axlim(2)*pos.vert(:,2),varargin{:});

axis(ca,axlim);


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = DrawTriPatch(BinGroups,bins,Cdata,varargin)


if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

axlim = axis(ca);

varargout = cell(1,nargout);

vert = BinGroups.pos.vert;
tri = BinGroups.pos.tri;

v1 = vert(:,1);
v2 = vert(:,2);

% [varargout{:}] = plot(axlim(2)*BinGroups.pos(:,1),axlim(2)*BinGroups.pos(:,2),varargin{:});
[varargout{:}] = patch(v1(tri(bins,:))',v2(tri(bins,:))',Cdata,varargin{:});

axis(ca,axlim);


%%%%%%%%%%%%%%%%%%%%%%%%%

function binarea = areaTriBin(BinGroups)

%area (or volume) of the simplex

vert = BinGroups.pos.vert;
tri = BinGroups.pos.tri;
nd = size(vert,2);

binarea = zeros(1,size(tri,1));

for i = 1:length(tri)
    binarea(i) = 1./prod(1:3)*det(vert(tri(i,1:nd),:)-repmat(vert(tri(i,end),:),nd,1));
end






%%%%%%%%%%%%%%%%%%%%%%

function inside = isinside(BinGroups,points,varargin)

% function inside = isinside(BinGroups,points, gaps,varargin)
% %Determines if points are inside 
% if nargin < 2
%     gaps = [];
% end
% 
% if ~isempty(gaps)
%     insidegaps = isinside(gaps,points);
%     insidegaps = sum(cat(2,insidegaps{:}))>0;
% end

%decide which function to use to determine if points are inside
inside = cell(1,size(BinGroups,1));

for i = 1:length(BinGroups)
   
    type = BinGroups{i}.type;
    
    switch lower(type)

        case {'grid','rect'}

           inside{i}  = insideRect(BinGroups,points,varargin{:});

        case 'poly'

           inside{i}  = insidePolyBin(BinGroups,points,varargin{:});

        case {'ellipse','circle'}                
           inside{i}  = insideEllipseBin(BinGroups, points,varargin{:});

        case {'tri','simplex'}                
           inside{i}  = insideTriBin(BinGroups, points,varargin{:});

        case 'undefined'
            if ~isempty(BinGroups)
                error('Type must be given for non-empty data.')
            end
        otherwise

           error([type,' is not a valid bin type.']);
    end
%     insidepts{i}(insidegaps,:) = 0;
    
end



%%%%%%%%%%%%%%%%%%%

function area = binarea(BinGroups,gaps,varargin)

%compute bin area

if nargin < 2
    gaps = [];
end

area = cell(1,length(BinGroups));

for i = 1:length(BinGroups)
   
    type = BinGroups{i}.type;
    
    if isempty(gaps)
        switch lower(type)

            case {'grid','rect'}

               area{i}  = areaRect(BinGroups,varargin{:});

            case 'poly'

               area{i}  = areaPolyBin(BinGroups,varargin{:});

            case {'ellipse','circle'}                
               area{i}  = areaEllipseBin(BinGroups, varargin{:});
            case {'simplex','tri'}                
               area{i}  = areaTriBin(BinGroups, varargin{:});
            case 'undefined'
                if ~isempty(BinGroups)
                    error('Type must be given for non-empty data.')
                end
            otherwise

               error([type,' is not a valid bin type.']);
        end
    else
       
        binareaWithGaps(BinGroups,gaps,varargin{:})
        
    end
       
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function binareaWithGaps(BinGroups,gaps,varargin)
% 
% error('This function has not been finished')
%     
% for i = 1:length(BinGroups)
%     
%     
% end