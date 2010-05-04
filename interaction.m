function RX = interaction(varargin)

% function RX = interaction(R1,R2)
%
%   Creates a regressor of interaction terms between R1 and R2. 
% 
% function RX = interaction(R,'intxnord',intxnord)
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

intxnord = 2;

if isstruct(varargin{1})
    R1 = varargin{1};
    varargin(1)=[];
end


if ~isempty(varargin) && isstruct(varargin{1})
    R2 = varargin{1};
    varargin(1)=[];
elseif ~isempty(varargin) && isnumeric(varargin{1})
%     R2 = varargin{1};
    intxnord = varargin{1};
else
    R2 = 0;
end



i=1;
codeincr = 0;
label = [];
maxord = Inf;
while i <= length(varargin)
    
    switch lower(varargin{i})
        case 'intxnord' %order of interaction 
            intxnord = varargin{i+1};
            i = i+1;            
        case 'codeincr'
            codeincr = varargin{i+1};
            i = i+1;            
        case 'maxord'  %Returns only interactions of order less than or equal to this term, as determined from the info.intxnord field.
            maxord = varargin{i+1};
            i = i+1;            
        case 'label'
            label = varargin{i+1};
            i = i+1;            
        otherwise
            error('%s is not a valid keyword.',varargin{i})
    end
    i = i+1;
    
end

if intxnord == 1 %interaction order 1 implies no terms other than the ones already in the input
    RX = R1([]);
    return
end


if nargin > 1 && isstruct(R2)
    %For two regressor arrays 2nd order interactions between groups
    indx = 1;
    for i = 1:length(R1)
        for j = 1:length(R2)

            if any(R1(i).noptions(:) ~= R2(j).noptions(:))
                error('Blocks do not appear to be identical for all regressors')
            end
            int = interaction2(R1(i), R2(j),maxord);
            if ~isempty(int)
                RX(indx) = int;
                RX(indx).code = codeincr + indx;
                RX(indx).codevec(1:RX(indx).Npar) = codeincr + indx;
                if ~isempty(label)
                    RX(indx).label = sprintf('%s %i',label,indx);
                end
                indx = indx+1;
            end
        end
    end
    
    return
    
elseif length(R1) == 2
    
    RX = interaction(R1(1),R1(2),'codeincr',codeincr);
    RX.code = codeincr + 1;
    if ~isempty(label)
         RX.label = label;
    end
    return
    
elseif length(R1) == 1 && (nargin == 1  || isempty(R2) )
    
    RX = R1;
    return
    
end

%For one array, all interactions among members are computed up to specified
%interaction order.



indx = 1;
for intxn = 2:intxnord
    B = chooseperm(length(R1),intxn);
    
    for i = 1:size(B,1)
        RX(indx) = interaction(R1(B(i,1)), interaction( R1(B(i,2:end)) ) );
        RX(indx).code = codeincr + indx;
        RX(indx).codevec(1:RX(indx).Npar) = codeincr + indx;
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

function RX = interaction2(R1,R2,maxord)


rxind = 1;


RX = makeregressor([]);
info = RX.info;

RX = RX([]);
    
for i = 1:length(R1)
    if strcmp(R1.info.form,'sparse')
        Vi = unsparsify(R1.value,'transpose');
    else
        Vi = R1.value;
    end        
    
    for j = 1:length(R2)
        
        if R1(i).info.intxnord + R2(j).info.intxnord > maxord  %only return regressors of order less than maxord
            continue
        end
        RX(rxind).info = info;
        if strcmp(R2.info.form,'sparse')
            Vj = unsparsify(R2.value,'transpose');
        else
            Vj = R2.value;
        end        
    
        
%         RX(rxind).value = sparseblock(repmat( Vi,1,size( Vj,2)).*kron( Vj ,ones(1,size( Vi,2))), R1(i).noptions,'transpose');     
%         RX(rxind).levmat = cat(1,repmat( R1.levmat,1,size( Vj,2)), kron( R2.levmat ,ones(1,size( Vi,2))));
%         RX(rxind).factmat = cat(1,repmat( R1.factmat,1,R2.Npar), kron( R2.factmat + 1,ones(1,R1.Npar)));

        if strcmp(R2.info.form,'sparse') ||strcmp(R1.info.form,'sparse')  
            RX(rxind).value = sparseblock(kron( Vi,ones(1,size( Vj,2))).*repmat( Vj ,1,size( Vi,2)), R1(i).noptions,'transpose');     
        elseif ~isempty(Vi) && ~isempty(Vj)
            RX(rxind).value = kron( Vi,ones(1,size( Vj,2))).*repmat( Vj ,1,size( Vi,2));     
        else
            RX(rxind).value =[];
            RX(rxind).function = makemefun(R1.function,R2.function,R1.Npar,R2.Npar);
%             continue
        end
        RX(rxind).levmat = cat(1,kron( R1.levmat,ones(1,R2.Npar)), repmat( R2.levmat ,1,R1.Npar));
        RX(rxind).factmat = cat(1,kron( R1.factmat,ones(1,R2.Npar)), repmat( R2.factmat,1,R1.Npar));
        RX(rxind).normconst = R1(i).normconst*R2(j).normconst; 
        RX(rxind).label = sprintf('%s * %s',R1(i).label,R2(j).label) ;
        RX(rxind).info.label = sprintf('%s * %s',R1(i).label,R2(j).label) ;
        RX(rxind).info.intxnord = R1(i).info.intxnord + R2(j).info.intxnord;
        RX(rxind).info.parent = cat(1,R1(i).info,R2(j).info);
%         RX(rxind).info.hashcode = R1(i).info.hashcode./2 + R2(j).info.hashcode;
        if strcmp(R2.info.form,'sparse') ||strcmp(R1.info.form,'sparse')  
            RX(rxind).info.form = 'sparse';
        else
            RX(rxind).info.form = R1.info.form ;
        end
        RX(rxind).noptions= R1(i).noptions;
        RX(rxind).Npar = R1(i).Npar*R2(j).Npar;
        RX(rxind).fixed = zeros(1,RX(rxind).Npar);
        RX(rxind).function = makemefun(R1.function,R2.function,R1.Npar,R2.Npar);
 
        rxind = rxind+1;
    end
end

%%%%%%%%%%%%%%%%%

function mkfun = makemefun(R1fun,R2fun,npar1,npar2)


if isempty(R1fun) || isempty(R2fun)
    mkfun = [];
    return
end

nargs1 = nargin(R1fun);
nargs2 = nargin(R2fun);

argliststr = strcat('X',cellfun(@num2str,num2cell(1:nargs1+nargs2),'uniformOutput',0));

funstr = sprintf('@(%s) kron( R1fun(%s),ones(1,npar2) ).*repmat( R2fun(%s) , 1 , npar1 );',...
            [sprintf('%s,',argliststr{1:end-1}),argliststr{end}],...
            [sprintf('%s,',argliststr{1:nargs1-1}),argliststr{nargs1}],...
            [sprintf('%s,',argliststr{(1:nargs2-1)+nargs1}),argliststr{nargs1+nargs2}]);

mkfun = eval(funstr);



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

