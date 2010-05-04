function RX = interaction(R1,R2,varargin)

% function RX = interaction(R1,R2)
%
%   Creates a regressor of interaction terms between R1 and R2. 
% 
% function RX = interaction(R,intxnord)
%
%   Creates interactions among all the regressors within the regressor
%   array, R, up to order intxnord.
%
% function RX = interaction(R)
%
%   Creates 2nd order interactions among all the regressors within the regressor
%   array, R.
%
% See also MAKEREGRESSOR


i=1;
while i <= length(varargin)
    
    switch lower(varargin{i})
    otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;
    
end



if nargin > 1 && ~isnumeric(R2)
  
    if any(R1.noptions ~= R2.noptions)
        error('Blocks do not appear to be identical for all regressors')
    end
    RX = interaction2(R1,R2);    

    return
    
elseif length(R1) == 2
    
    RX = interaction(R1(1),R1(2));
    return
    
elseif length(R1) == 1 && (nargin == 1  || isempty(R2) )
    
    RX = R1;
    return
    
end

if nargin == 1 || isempty(R2)
    intxnord = 2;
else 
    intxnord = R2;
end

indx = 1;
for intxn = 2:intxnord
    B = chooseperm(length(R1),intxn);
    
    for i = 1:size(B,1)
        RX(indx) = interaction(R1(B(i,1)), interaction( R1(B(i,2:end)) ) );
        indx = indx+1;
    end
    
end

% if nargin < 3
%     ixnlev = 2;
% end

try
    fid = fopen([mfilename,'.m'],'r');
    R.info.COMMAND =  fread(fid);
    fclose(fid);
catch
    warning('Failed to record the script used to generate regressor %s poly',R.label);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RX = interaction2(R1,R2)
rxind = 1;


RX = makeregressor([]);

    
for i = 1:length(R1)
    if strcmp(R1.info.form,'sparse')
        Vi = unsparsify(R1.value,'transpose');
    else
        Vi = R1.value;
    end        
    
    for j = 1:length(R2)
        if strcmp(R2.info.form,'sparse')
            Vj = unsparsify(R2.value,'transpose');
        else
            Vj = R2.value;
        end        
    
        
%         RX(rxind).value = repmat(RX(i).value,size( RX(j).value,2)).*kron(RX(i).value,ones(1,size( RX(i).value,2)));
        RX(rxind).value = sparseblock(repmat( Vi,1,size( Vj,2)).*kron( Vj ,ones(1,size( Vi,2))), R1(i).noptions,'transpose');        
        RX(rxind).normconst = R1(i).normconst*R2(j).normconst; 
        RX(rxind).label = sprintf('%s * %s',R1(i).label,R2(j).label) ;
        RX(rxind).info.label = sprintf('%s * %s',R1(i).label,R2(j).label) ;
        RX(rxind).info.intxnord = R1(i).info.intxnord + R2(j).info.intxnord;
        RX(rxind).info.parent = cat(1,R1(i).info,R2(j).info);
        RX(rxind).info.hashcode = R1(i).info.hashcode./2 + R2(j).info.hashcode;
        RX(rxind).info.form = 'sparse';
        RX(rxind).noptions= R1(i).noptions;
        RX(rxind).Npar = R1(i).Npar*R2(j).Npar;
        
        rxind = rxind+1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = chooseperm(N, k, item, rows)

% 
%    Returns indices for every unordered combination of k items from a
%    population of N using a recursive algorithm.


if nargin < 3
    rows = 1;
    item = 1;
end

B = [];

if item >= k
    
    B = (rows:N)';
    
else
    
    for strow = rows:N

       b = chooseperm( N, k, item+1, strow+1);
    
       if ~isempty(b)
           
           B = cat(1,B,cat(2,ones(size(b,1),1)*strow,b));
           
       else
           
           return
           
       end
    end    
end

