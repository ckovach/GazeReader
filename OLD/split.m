function Rout = split(R,grps)

% function Rout = split(R,grp)
% 
%   Splits the block R.value into separate blocks of regressors, where
%   block i has grp(i) parameters.
%


Rout = makeregressor([]);



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
    

    if iscell(grps)
        grp = grps{k};
    elseif nargin < 2
        grp = [];
    else
        grp = grps;
    end
    
    if  isempty(grp) 
        grp = ones(1,R(k).Npar);
    end
    
    stindex = cumsum([0,grp(1:end-1)]);

    
    if strcmp(R(k).info.form,'sparse')
        V = unsparsify(R(k).value,'transpose');
    end

        
    for i = 1:length(grp)   

        grpvec = stindex(i)+ (1:grp(i));

        Rout(indx).value = sparseblock(V(:,grpvec),R(k).noptions, 'transpose');
        Rout(indx).normconst = R(k).normconst; 
        
        if length(grp) == 1 && grp == R(k).Npar;
            Rout(indx).label = R(k).label;
        else
            Rout(indx).label = sprintf('%s_sub_%i',R(k).label,i);
        end  
        
        Rout(indx).info.label = sprintf('%s_sub_%i',R(k).label,i);

        Rout(indx).info.COMMAND = COMMAND;
        Rout(indx).info.parent = R(k).info;
        Rout(indx).info.form = R(k).info.form;
        Rout(indx).info.intxnord = R(k).info.intxnord;
        Rout(indx).noptions = R(k).noptions;
        Rout(indx).Npar = grp(i);

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
    
