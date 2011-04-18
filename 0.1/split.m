function Rout = split(R,grps,codeincr)

% function Rout = split(R,grp)
% 
%   Splits the block R.value into separate blocks of regressors, where
%   block i has grp(i) parameters.
%
% See also POOL, MAKEREGRESSOR

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


Rout = makeregressor([]);


if nargin < 3
%     codeincr = [];
    codeincr = max([R.code]);
end

COMMAND= '';
try
    fid = fopen([mfilename,'.m'],'r');
    COMMAND =  fread(fid);
    fclose(fid);
catch
    warning('Failed to record the command.');
end


indx = 1;
for k = 1:length(R)
    
    

    if nargin < 2
        grp = [];
    elseif iscell(grps)
        grp = grps{k};
    else
        grp = grps;
    end
    
    if  isempty(grp) 
        grp = ones(1,R(k).Npar);
    end
    
    stindex = cumsum([0,grp(1:end-1)]);

    
    if strcmp(R(k).info.form,'sparse')
        V = unsparsify(R(k).value,'transpose');
    else
        V = R(k).value;
    end
        
    for i = 1:length(grp)   

        grpvec = stindex(i)+ (1:grp(i));
        
        Rout(indx) = makeregressor([],'codeincr',i-1+codeincr,'noptions',R(k).noptions);

        if strcmp(R(k).info.form,'sparse')
            Rout(indx).value = sparseblock(V(:,grpvec),R(k).noptions, 'transpose');
        else
            Rout(indx).value =V(:,grpvec);
        end
        
        
        Rout(indx).normconst = R(k).normconst; 

         Rout(indx).levmat = R(k).levmat(:,grpvec);
         Rout(indx).factmat = R(k).factmat(:,grpvec);
         
        if length(grp) == 1 && grp == R(k).Npar;
            Rout(indx).label = R(k).label;
        else
            Rout(indx).label = sprintf('%s_sub_%i',R(k).label,i);
        end  
        
%         Rout(indx).info = R(1).info;
        Rout(indx).info.label = sprintf('%s_sub_%i',R(k).label,i);

        Rout(indx).info.COMMAND = COMMAND;
        Rout(indx).info.parent = R(k).info;
        Rout(indx).info.form = R(k).info.form;
        Rout(indx).info.intxnord = R(k).info.intxnord;
%         Rout(indx).noptions = R(k).noptions;
        Rout(indx).Npar = grp(i);
%         Rout(indx).code = codeincr+i;
        Rout(indx).codevec = Rout(indx).code*ones(1,grp(i));
        
        Rout(indx).info.hashcode = 0;
        for j = 1:length(grpvec )
            
           reseed; 
           r32 = typecast(rand,'uint32');
           
           newhash = r32(1);
%            if newhash == 0 
%                newhash = Rout(indx).info.hashcode;
%            end

           Rout(indx).info.hashcode = newhash;
                          
        end

        indx = indx+1;
    end
end
    
