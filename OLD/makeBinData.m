
function binData = makeBinData(data,varargin)

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
% Additional options:
% 
%       makeBinData(...,'label',binlabel) assigns a label to the group of
%       bins.
% 
%       makeBinData(...,'trials',trialvec) specifies which trials a bin is
%       active,  or the number of sequential trials if trialvec is
%       a scalar.

i = 1;
type = 'undefined';
trials = [];
ntrials = 0;
% units = 'normalized';
label = '';
precedence = 'nearest_to_center';
while i <= length(varargin)
   switch lower(varargin{i})
       case 'type'
          type = varargin{i+1};
          i = i+1;
       case 'trials'
          trials = varargin{i+1}; %1Xntrials logical vector with 1  whenever the to-be-assigned bin occurs
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


% binData = struct('groups',[],'activeTrials',[],'precedence','nearest_to_center');
binData = struct('groups',binGroupsTemplate,'precedence',precedence);
binData.groups.label = label;
binData.groups.type= type;
% binData.groups.inputData = data;

if isempty(data)
    return
end

if ~isempty(trials) && ntrials == 0;
    ntrials = max(trials);
end

switch lower(type)
    
    case 'grid'
        
        binData.groups = makeGridbin(data);
        
    case 'poly'
        
        binData.groups = makePolybin(data);
        
    case 'rect'
        
        binData.groups = makeRectbin(data);

    case 'ellipse'                
        binData.groups = makeEllipsebin(data);

    case 'circle'                 
        binData.groups = makeEllipsebin(data);

    case 'undefined'
        if ~isempty(data)
            error('Type must be given for non-empty data.')
        end
    otherwise

       error([type,' is not a valid bin type.']);
end

binData.groups.nbin = size(binData.groups.pos,1);

binData.groups.inputData = data;

binData.groups.label = label;
% multassign = @ (varargin) varargin{:};

if all(ismember(trials,[0 1]))
%     binData.groups.activeTrials = find(trials(:));
    binData.groups.activeTrials =  trials(:);
 
else 
%     binData.groups.activeTrials = trials;
        binData.groups.activeTrials = false( ntrials , 1 );
        binData.groups.activeTrials( trials ) = true;
end
% if all(ismember(trials,[0 1]))
%     binData.activeTrials = trials(:);
% else 
%     binData.activeTrials = false( ntrials , 1 );
%     binData.activeTrials( trials ) = true;
% end




%%%%%%%%%%%%%%%%%%%%%%%
function BinGroups = binGroupsTemplate
% initializes an empty BinGroups structure

BinGroups = struct('type',[],'label','','binnums',[],'pos',[],'nbin',[],'centers',[],'ismember',[],'plot',[],'patch',[],'activeTrials',[],'inputData',[]);


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

%%%%%%%%%%%%%%%%%%%%%

function BinGroups = makeRectbin(rect)

%Create a set of rectangular bins

BinGroups = binGroupsTemplate;
BinGroups.type = 'rect';

BinGroups.binnums = 1:size(rect,1); 

BinGroups.pos = rect;

BinGroups.centers = cat(2,sum(rect(:,1:2),2)./2,sum(rect(:,3:4),2)./2);


BinGroups.ismember = @(BinGroups,FixLoc,tol) insideRect(BinGroups,FixLoc,tol);
BinGroups.plot = @(BinGroups,varargin) PlotRectBin(BinGroups,varargin{:});
BinGroups.patch = @(BinGroups,varargin) DrawRectPatch(BinGroups,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = PlotRectBin(binGroups,bins,varargin)

if ~ ( strcmp(binGroups.type,'rect')) 
    error('This plot function is valid for bin''s of type rect only.')
end

if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

axlim = axis(ca);


ps = plot(axlim(2)*binGroups.pos(bins,[1 1])', axlim(4)*binGroups.pos(bins,[3 4])',varargin{:});

hold on,

ps = cat(2,ps,plot(axlim(2)*binGroups.pos(bins,[2 2])', axlim(4)*binGroups.pos(bins,[3 4])',varargin{:}));
ps = cat(2,ps,plot(axlim(2)*binGroups.pos(bins,[1 2])', axlim(4)*binGroups.pos(bins,[3 3])',varargin{:}));
ps = cat(2,ps,plot(axlim(2)*binGroups.pos(bins,[1 2])', axlim(4)*binGroups.pos(bins,[4 4])',varargin{:}));

axis(ca,axlim);

if nargout > 0
    varargout{1} = ps;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = DrawRectPatch(binGroups,bins,Cdata,varargin)

% Draws a bin shaped patch
if ~ ( strcmp(binGroups.type,'rect')) && ~ ( strcmp(binGroups.type,'grid')) 
    error('This plot function is valid for bin''s of type rect and grid only.')
end

if nargin < 2 || isempty(bins)
    bins = 1:length(binGroups.pos,1);
end

if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

axlim = axis(ca);

pos = binGroups.pos;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = PlotGridBin(binGroups,bins,varargin)


%This plots the entire grid -- bins does not affect anything and is
%included in the argument list for compatibility.

if ~ ( strcmp(binGroups.type,'grid') )
    error('This plot function is valid for bin''s of type grid only.')
end

if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

axlim = axis(ca);

x = binGroups.pos(:,1);
dx = diff(x);
dx = dx(find(dx,1));
if isempty(dx)    
    dx = find(binGroups.pos(:,4)- binGroups.pos(:,3),1);
end

y = binGroups.pos(:,3);
dy = diff(y);
dy = dy(find(dy,1));
if isempty(dy)
    dy = find(binGroups.pos(:,4)-binGroups.pos(:,3),1);
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

           
   


   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ellipse bin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function binGroups = makeEllipsebin(data)
 
% if isnumeric(data)
%     data= {data(:,1:2),data(:,3)};
% end
if iscell(data)
    data= [data{1},data{2}];
end

if ~isempty(data) && size(data,2) < 3 %if data has 3 columns, the y-semi axis will be the same as the x
    data(:,4) = data(:,3);
end

if ~isempty(data) && size(data,2) < 4 %rotation parameter set to 0 if unspecified
    data(:,5) = 0;
end

binGroups = binGroupsTemplate;

binGroups.type = 'ellipse';
    


binGroups.centers = data(:,1:2);
binGroups.pos = data;
binGroups.binnums = 1:size(data,1);


binGroups.ismember = @(binGroups,FixLoc,varargin) insideEllipseBin(binGroups,FixLoc,varargin{:});
    
binGroups.plot = @(binGroups,varargin) plotEllipseBin(binGroups,varargin{:});
binGroups.patch= @(binGroups,varargin) DrawEllipsePatch(binGroups,varargin{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binmember = insideEllipseBin(binGroups,FixLoc,tol)

if nargin < 3 || isempty(tol)
     tol = 0;
end

radii = binGroups.pos(:,3:4) + tol;
rotangle = binGroups.pos(:,5);
nbin = size(binGroups.pos,1);
nfix = size(FixLoc,1);

centers = binGroups.pos(:,1:2);

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

                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xy = ellipseVertices(binGroups,bins)
%Gets vertices of the ellipse bin

if ~ ( strcmp(binGroups.type,'ellipse') )
    error('This plot function is valid for bin''s of type ellipse only.')
end

th = 0:.01:2*pi;

centers = binGroups.pos(bins,1:2);

if size(binGroups.pos,2) < 4
    radii = binGroups.pos(bins,[3 3]) + eps;
else
    radii = binGroups.pos(bins,3:4) + eps;
end

if size(binGroups.pos,2) < 5
    rotangle = 0+eps;
else
    rotangle =   binGroups.pos(bins,5)+eps;
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
function varargout = plotEllipseBin(binGroups,bins,varargin)

if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

if nargin < 2 || isempty(bins)
    bins = 1:length(binGroup.pos,1);
end

axlim = axis(ca);

xy = ellipseVertices(binGroups,bins);

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

if nargin < 2 || isempty(bins)
    bins = 1:length(binGroup.pos,1);
end

axlim = axis(ca);

xy = ellipseVertices(binGroup,bins);

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



BinGroups = binGroupsTemplate;

BinGroups.type = 'poly';

BinGroups.pos = data;

BinGroups.centers = mean(data);

BinGroups.binnums = 1;

BinGroups.ismember = @(binGroups,FixLoc,tol)insidePolyBin(binGroups,FixLoc,tol);
BinGroups.plot = @(binGroups,varargin) plotPolyBin(binGroups,varargin{:});
BinGroups.patch= @(binGroups,varargin) DrawPolyPatch(binGroups,varargin{:});
%%%%%%%%%%%%%

function binmember = insidePolyBin(binGroups,FixLoc,tol)

%tol doesn't do anything for polygonal bins

binmember = inpolygon(FixLoc(:,1),FixLoc(:,2),binGroups.pos(:,1),binGroups.pos(:,2));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = plotPolyBin(binGroups,bins,varargin)

%bins does nothing and is only in the argument list for compatibility

if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

axlim = axis(ca);

varargout = cell(1,nargout);

[varargout{:}] = plot(axlim(2)*binGroups.pos(:,1),axlim(2)*binGroups.pos(:,2),varargin{:});

axis(ca,axlim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = DrawPolyPatch(binGroups,bins,Cdata,varargin)

%bins does nothing and is only in the argument list for compatibility

if any(strcmpi(varargin,'parent'))
    ca = varargin{find(strcmpi(varargin,'parent'))+1};
else
    ca = gca;
end

axlim = axis(ca);

varargout = cell(1,nargout);

% [varargout{:}] = plot(axlim(2)*binGroups.pos(:,1),axlim(2)*binGroups.pos(:,2),varargin{:});
[varargout{:}] = patch(binGroups.pos(:,1),binGroups.pos(:,2),Cdata,varargin{:});

axis(ca,axlim);
