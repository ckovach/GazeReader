function Rout = pool(varargin)

% function Rout = pool(R,grp)
% 
%   Pools members of a regressor array R into blocks as determined by grp, where grp is a 1xB
%   cell or vector with the number of sequentially ordered members in the fomer case and
%   the vectors of members in the latter case. By default all of R is gouped.
% 
% See also SPLIT, MAKEREGRESSOR

if ~isnumeric(varargin{end}) || ~iscell(varargin{end})
    R = cat(2,varargin{1:end});
    
    grp = {1:length(R)};
    
elseif ~iscell(varargin{end})
    
    grp = varargin{end};
    R = cat(1,varargin{1:end-1});

    stindex = cumsum([0,grp(1:end-1)]);
    for i = 1:lengh(grp)
        cellgrp{i} = stindex(i)+ (1:grp(i));
    end
    grp = cellgrp;
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
    noptions = R(1).noptions;
    V = [];
    for j = 1:length(grp{i})
        if any(R(j).noptions ~= noptions) 
            error('Not all regressors share the same trial structure!')
        end
        
        if strcmp(R(grp{i}(j)).info.form, 'sparse')
            V = cat(2,V,unsparsify(R(grp{i}(j)).value,'transpose'));        
        else
            V = cat(2,V,R(grp{i}(j)).value);
        end
    end
    
    Rout(i).value = sparseblock(V,noptions,'transpose');
    Rout(i).noptions = noptions;
    
    Rout(i).normconst = unique([R(grp{i}).normconst]); 
        
    Rout(i).label = sprintf('pooled');
    Rout(i).info.label = sprintf('pooled');
            
    Rout(i).info.COMMAND = COMMAND;
    Rout(i).info.parent = cat(1,R(grp{i}).info);
    Rout(i).info.hashcode = 0;
    Rout(i).Npar = sum([R(grp{i}).Npar]);
    Rout(i).info.form = 'sparse';
    
    strow = 0;
    for j = 1:length(grp{i})
        
        Rout(i).info.contrasts{j}  = zeros(Rout(i).Npar,R(grp{i}(j)).Npar);
        Rout(i).info.contrasts{j}(strow + ( 1 : R( grp{i}(j) ).Npar ) ,:) = eye( R( grp{i}(j) ).Npar );
        Rout(i).info.hashcode =  bitxor(Rout(i).info.hashcode , bitcmp(R(grp{i}(j)).info.hashcode));
        
        pooled{j}  = R(grp{i}(j)).info.pooled_labels;
        
        strow = strow+R(grp{i}(j)).Npar;
    end
    
    
    Rout(i).info.pooled_labels = cat(2,pooled{:}, { R(grp{i}).label } );
    
    
end

    
    