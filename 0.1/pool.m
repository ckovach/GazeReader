function Rout = pool(varargin)

% function Rout = pool(R,grp)
% 
%   Pools members of a regressor array R into blocks as determined by grp, where grp is a 1xB
%   cell or vector with the number of sequentially ordered members in the fomer case and
%   the vectors of members in the latter case. By default all of R is gouped.
% 
% See also SPLIT, MAKEREGRESSOR

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


i=1;
codeincr = 0;
poolinput=false;
label = [];
while i <= length(varargin)
    
    if ischar(varargin{i})
        switch lower(varargin{i})
            case 'codeincr'
                codeincr = varargin{i+1};
%                 i = i+1;
                varargin(i:i+1) = [];
                i = i-1;
            case 'label'
                label = varargin{i+1};
%                 i = i+1;
                varargin(i:i+1) = [];
                 i = i-1;
%             case 'poolinput'        %function takes only one input
%                 poolinput = varargin{i+1};
%                 i = i+1;
            otherwise
                error('%s is not a valid keyword.',varargin{i})
        end
%         varargin(i) = [];
    end
    i = i+1;
    
end

if ~isnumeric(varargin{end}) && ~iscell(varargin{end})
    R = cat(2,varargin{1:end});
    
    grp = {1:length(R)};
    
elseif ~iscell(varargin{end})
    
    grp = varargin{end};
    R = cat(1,varargin{1:end-1});

    stindex = cumsum([0,grp(1:end-1)]);
    for i = 1:length(grp)
        cellgrp{i} = stindex(i)+ (1:grp(i));
    end
    grp = cellgrp;
elseif iscell(varargin{end})
        R = cat(1,varargin{1:end-1});
        grp = varargin{end};
else
        R = cat(1,varargin{1:end-1});
end


COMMAND= '';
try
    fid = fopen([mfilename,'.m'],'r');
    COMMAND =  fread(fid);
    fclose(fid);
catch
    warning('Failed to record the script used to generate regressor %s poly',R.label);
end

Rout = makeregressor([]);

for i = 1:length(grp)   
    noptions = R(grp{i}(1)).noptions;
    V = [];
    sparseform = 0;
    for j = 1:length(grp{i})
        if ~isequal(R(grp{i}(j)).noptions,noptions) 
            error('Not all regressors share the same trial structure!')
        end
        
        if strcmp(R(grp{i}(j)).info.form, 'sparse')
            V = cat(2,V,unsparsify(R(grp{i}(j)).value,'transpose'));
            sparseform = true;
        else
            V = cat(2,V,R(grp{i}(j)).value);
        end
    end
    
    if sparseform
        Rout(i).value = sparseblock(V,noptions,'transpose');
        Rout(i).info.form = 'sparse';
    else
        Rout(i).value = V;
        Rout(i).info.form = R(grp{i}(1)).info.form;
    end
    Rout(i).noptions = noptions;
    
    Rout(i).normconst = unique([R(grp{i}).normconst]); 
    
    if isempty(label)
        Rout(i).label = sprintf('pooled');
        Rout(i).info.label = sprintf('pooled');
    else 
        Rout(i).label = label;
    end
    
    Rout(i).info.COMMAND = COMMAND;
    try
        Rout(i).info.parent = cat(1,R(grp{i}).info);
    catch
    end
    Rout(i).info.hashcode = 0;
    Rout(i).Npar = sum([R(grp{i}).Npar]);
    
    emptyfixed = cellfun(@isempty,{R(grp{i}).fixed});
    if any(emptyfixed)
        Rout(i).fixed = [];
    else        
        Rout(i).fixed= [R(grp{i}).fixed];
    end
    
    try
%         Rout(i).function = makefunction(R(grp{i}).function);      
        [Rout(i).function,terms,Rout(i).deriv] = makePooledFunction(R(grp{i}));      
        Rout(i).info.functionInputCodes = terms;
catch
        warning('Unable to generate a function for %s',Rout(i).label)
    end
%     try
%         Rout(i).deriv = makefunction(R(grp{i}).deriv);      
%     catch
%         warning('Unable to generate a function for %s',Rout(i).label)
%     end
    strow = 0;
     Rout(i).info.parentindex=[];
    for j = 1:length(grp{i})
        
        xord = size(R(grp{i}(j)).levmat,1);
        Rout(i).levmat(1:xord,end+(1:R(grp{i}(j)).Npar)) = R(grp{i}(j)).levmat ;
        Rout(i).factmat(1:xord,end+(1:R(grp{i}(j)).Npar)) = R(grp{i}(j)).factmat;

        Rout(i).info.contrasts{j}  = zeros(Rout(i).Npar,R(grp{i}(j)).Npar);
        Rout(i).info.contrasts{j}(strow + ( 1 : R( grp{i}(j) ).Npar ) ,:) = eye( R( grp{i}(j) ).Npar );
%         Rout(i).info.hashcode =  bitxor(Rout(i).info.hashcode , bitcmp(R(grp{i}(j)).info.hashcode));
        Rout(i).info.parentindex(1,end+(1:R(grp{i}(j)).Npar)) = grp{i}(j)*ones(1,R(grp{i}(j)).Npar);
        
%         pooled{j}  = R(grp{i}(j)).info.pooled_labels;
        
        strow = strow+R(grp{i}(j)).Npar;
    end
    
    Rout(i).code = codeincr+1;
    Rout(i).codevec(1:Rout(i).Npar) = codeincr+1;
    codeincr = codeincr+1;
%     Rout(i).info.pooled_labels = cat(2,pooled{:}, { R(grp{i}).label } );
    
    
end

%---------------------------------------    
function fun = makefunction(varargin)

if any(cellfun(@isempty, varargin))
    fun = [];
    return
end

nargs = cellfun(@nargin,varargin);
csnargs = [0,cumsum(nargs)];

funstr = ['@(X1',sprintf(',X%i',2:sum(nargs)),') cat(2'];
for i = 1:length(varargin)
    if nargs(i) >=2
        str = sprintf(',X%i',csnargs(i)+(2:nargs(i)));
    else
        str = '';
    end
    funstr = cat(2,funstr,sprintf(',varargin{%i}(X%i%s)',i,csnargs(i)+1,str));
end
funstr = [funstr,');'];

fun = eval(funstr);

%--------------------------------------
function [pooled_function,terms,dpooled] = makePooledFunction(Rs)



%Create a function which distributes inputs to the appropriate arguments


% Get regressor codes
rcodes = [Rs.code];
rfunctions = {Rs.function};
factmats = {Rs.factmat};
dfunctions = {Rs.deriv};


terms= zeros(2,0);
for i = 1:length(Rs)
    [unq,q,polyterm] = unique(factmats{i}(:,1),'stable');
     rinputs{i} = [];
     for j = 1:length(unq)
         for k = 1:sum(polyterm==j)
             
             input = find(terms(1,:) == unq(j) & terms(2,:) == k);
             if isempty(input)
                 terms(:,end+1) = [unq(j),k];
                 rinputs{i}(end+1) = size(terms,2);
             else
                 rinputs{i}(end+1) = input;
             end
                 
         end
     end
end             

nargs = max([rinputs{:}]);
args = cellfun(@(n) sprintf('X%i',n),num2cell(1:nargs),'uniformoutput',false);
args2 = cellfun(@(n) sprintf('%i',n),num2cell(1:nargs),'uniformoutput',false);
trim = @(a) a(1:end-1);
genarg = @(X)trim(sprintf('%s,',X{:}));
funs = cellfun(@(inp,i) sprintf('rfunctions{%i}(%s)',i,genarg(args(inp))),rinputs,num2cell(1:length(rinputs)),'uniformoutput',false);

dfuns = cellfun(@(inp,i) sprintf('dfunctions{%i}(%s,n([%s]))',i,genarg(args(inp)),genarg(args2(inp))),rinputs,num2cell(1:length(rinputs)),'uniformoutput',false);

pooled_function = str2func( sprintf('@(%s) cat(2, %s )',genarg(args),genarg(funs)));

try
    dpooled = str2func( sprintf('@(%s,n) cat(2, %s )',genarg(args),genarg(dfuns)));
catch
    dpooled = [];
end
