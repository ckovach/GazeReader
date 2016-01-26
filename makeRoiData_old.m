
function binData = makeBinData(data,varargin)

% binData = makeBinData(...)
%
%   Creates an a region of interest data structure with the following
%   fields:
% 
%     binData.instances - a structure array where each element definines an bin or groups of bins.
% 
%     binData.trialIndices - a vector of inices into the array contained in the "instances" field. 
%                             The values of Indices(i) indicates that
%                             the indices(i)'th bin is active on the i'th trial
% 
%     binData.precedence - a field indicating how membership in overlapping
%                   bins is to be resolved. This may be
%                        'nearest_to_center':  bin with the center nearest to the fixation point wins
%                        'array_order':   The bin specified earliest in the instances array wins.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------





i = 1;
type = 'undefined';
trials = [];
ntrials = 0;
units = 'normalized';
label = '';
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
       case 'screenunits'
           units = 'screen';
       case 'label'
           label  = varargin{i+1}; %Total number of trials
            i = i+1;
        otherwise

           error([varargin{i},' is not a valid option.']);
   end         
   i = i+1;
end


binData = struct('instances',[],'trialIndices',[],'precedence','nearest_to_center');

if isempty(data)
    return
end

if ~isempty(trials) && ntrials == 0;
    ntrials = max(trials);
end

switch lower(type)
    
    case 'grid'
        
        binData.instances = makeGridbin(data);
        
    case 'poly'
        
        binData.instances = makePolybin(data);
        
    case 'rect'
        
        binData.instances = makeRectbin(data);

    case 'circle'        
        
        binData.instances = makeCiclebin(data);

    case 'undefined'
        if ~isempty(data)
            error('Type must be given for non-empty data.')
        end
    otherwise

       error([type,' is not a valid bin type.']);
end

binData.instances.nbin = size(binData.instances.pos,1);

binData.instances.label = label;
if all(ismember(trials,[0 1]))
    binData.trialIndices = trials(:);
else 
    binData.trialIndices = false( ntrials , 1 );
    binData.trialIndices( trials ) = true;
end




%%%%%%%%%%%%%%%%%%%%%%%
function BinInstances = binInstancesTemplate
% initializes an empty BinInstances structure

BinInstances = struct('type',[],'label',[],'binnums',[],'pos',[],'nbin',[],'centers',[],'ismember',[],'plot',[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rect and Grid bins
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%

function BinInstances = makeGridbin(data)


if ~iscell(data)
    error('Grid bin requires cell input')
end    
if length(data)~=2
    error('Input must be a cell array of length 2 for Grid bin')
end


lx = linspace(data{1}(1), data{1}(2), data{2}(1)+1);
lxoffset = (data{1}(2) - lx(end))./2;
lx = lx+lxoffset;

ly = linspace(data{1}(3), data{1}(4), data{2}(2)+1);
lyoffset = (data{1}(4) - ly(end))./2;
ly = ly+lyoffset;

[X,Y] = meshgrid(lx(1:end-1),ly(1:end-1));

[dx,dy] = meshgrid(diff(lx),diff(ly));


rect = cat(2 , X(:) , X(:) + dx(:) ,  Y(:), Y(:) + dy(:) );

BinInstances = MakeRectbin(rect);
BinInstances.type = 'grid';
BinInstances.plot = @(BinInstances,varargin) PlotGridBin(BinInstances,varargin{:});

%%%%%%%%%%%%%%%%%%%%%

function BinInstances = MakeRectbin(rect)

%Create a set of rectangular bins

BinInstances = binInstancesTemplate;
BinInstances.type = 'rect';

BinInstances.binnums = 1:size(rect,1); 

BinInstances.pos = rect;

BinInstances.centers = cat(2,sum(rect(:,1:2),2)./2,sum(rect(:,3:4),2)./2);


BinInstances.ismember = @(FixLoc) insideRect(BinInstances,FixLoc, rect);
BinInstances.plot = @(BinInstances,varargin) PlotRectBin(BinInstances,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotRectBin(binInstances,varargin)

if ~ ( strcmp(binInstances.type,'rect')) 
    error('This plot function is valid for bin''s of type rect only.')
end


plot(binInstances.pos(:,[1 1])', binInstances.pos(:,[3 4])',varargin{:});

hold on,

plot(binInstances.pos(:,[2 2])', binInstances.pos(:,[3 4])',varargin{:});
plot(binInstances.pos(:,[1 2])', binInstances.pos(:,[3 3])',varargin{:});
plot(binInstances.pos(:,[1 2])', binInstances.pos(:,[4 4])',varargin{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotGridBin(binInstances,varargin)

if ~ ( strcmp(binInstances.type,'grid') )
    error('This plot function is valid for bin''s of type grid only.')
end


x = binInstances.pos(:,1);
dx = diff(x);
dx = dx(find(dx,1));

y = binInstances.pos(:,3);
dy = diff(y);
dy = dy(find(dy,1));

gridticky = y(1):dy:y(end)+dy; 
gridtickx = x(1):dx:x(end)+dx; 

plot( repmat(gridtickx,2,1) , repmat([gridticky(1) gridticky(end)]',1,length(gridtickx) ));

hold on,

plot(  repmat([gridtickx(1) gridtickx(end)]',1,length(gridticky) ), repmat(gridticky,2,1));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binmember = insideRect(BinInstances,FixLoc)


rect = BinInsances.rect;
nbin = size(BinInstances.rect,1);
nfix = size(FixLoc,1);

binmember =    sparse(repmat( FixLoc(:,1) , 1 , nbin ) >= repmat( rect(:,1)' , nfix , 1 )  & ...
               repmat( FixLoc(:,1) , 1 , nbin ) <  repmat( rect(:,2)' , nfix , 1 )  & ...
               repmat( FixLoc(:,2) , 1 , nbin ) >= repmat( rect(:,3)' , nfix , 1 )  & ...
               repmat( FixLoc(:,3) , 1 , nbin ) <  repmat( rect(:,4)' , nfix , 1 ));

           
   


   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Circle bin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function binInstances = makeCiclebin(data);

    
    
binInstances = binInstancesTemplate;

binInstances.type = 'circle';
    

binInstances.centers = data{1};
binInstances.pos = data{1};
binInstances.binnums = 1:size(data{1},1);

radii = data{2};

try

    binInstances.pos(:,end+1) = radii(:,1);
    if size(radii,2) == 2
        binInstances.pos(:,end+1) = radii(:,2);
    end
catch
        error('Number of circle radii doesn''t match number of circle centers.')
    
end

binInstances.ismember = @(binInstances,FixLoc) insideCircleBin(binInstances,FixLoc);
    
binInstances.plot = @(binInstances,varargin) plotCircleBin(binInstances,varargin{:})


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function binmember = insideCircleBin(binInstances,FixLoc)


radii = binInstances.pos(:,3);

centers = binInstances.pos(:,1:2);


binmember =    sparse( ( repmat( FixLoc(:,1) , 1 , nbin ) - repmat( centers(:,1)' , nfix , 1 ) ).^2  + ...
                       ( repmat( FixLoc(:,2) , 1 , nbin ) - repmat( centers(:,2)' , nfix , 1 ) ).^2 ... 
                        <= repmat( radii(:)'.^2 , nfix , 1 )  );

                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotCircleBin(binInstances,varargin)

if ~ ( strcmp(binInstances.type,'circle') )
    error('This plot function is valid for bin''s of type circle only.')
end

th = 0:.01:2*pi;

centers = binInstances.pos(:,1:2);

radii = binInstances.pos(:,3:4);



for i = 1:size(centers,1);    
    plot(radii(i,1)*cos(th) + centers(i,1), radii(i,2)*sin(th) + centers(i,2));
    hold on
end


                    
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polygonal bins - these can only represent 1 bin at a time
%%%%%%%%%%%%%%%%%%%%%%%%%%

function binInstances = makePolybin(data)



binInstances = binInstancesTemplate;

BinInstances.type = 'polygon';

BinInstances.pos = data;

BinInstances.centers = mean(data);

BinInstances.binnums = 1;

BinInstances.ismember = insidePolyBin(binInstances,FixLoc);
BinInstances.plot = @(binInstances,varargin) plotPolyBin(binInstances,varargin{:});
%%%%%%%%%%%%%

function binmember = insidePolyBin(binInstances,FixLoc);

binmember = inpolygon(FixLoc(:,1),FixLoc(:,2),binInstances.pos(:,1),binInstances.pos(:,2));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotPolyBin(binInstances,varargin)

plot(binInstances.pos(:,1),binInstances.pos(:,2),varargin{:});

