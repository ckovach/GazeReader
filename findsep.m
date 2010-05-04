function FS = findsep(h,model,data_set)

% function FS = findsep(h,model,data_set)
% Finds separation in the fit for a group of models and returns a structure with the index of the regressor
% group sorted by weighting along each dimension in the null space
% of the information matrix, as well as the weigting: FS.Sepmat, FS.usep,
% FS.member. Separation is detected when the ratio of the minimum to
% maximum eigenvalue of the information matrix falls below a secified
% tolerance.

tol = 1e-10; %tol

regData = getappdata(h,'regData');
modelData = getappdata(h,'modelData');

if nargin < 3
    data_set = 1:length(modelData);
end


FS = struct('model',[]);
FS(1) = [];
for i=data_set 

    if nargin < 2
        model = 1:length(modelData(i).models);
    end

    for j = model
        
        I = modelData(i).models(j).fit(1).I;
        
        [u,d] = svd(I);

        d = diag(d)./d(1);

        sep = d < tol;

        if ~any(sep)    
            continue
        end
        
       fprintf('\nSeparation detected in data set %i, model %i....',i,j)
        
       rgcodes= [regData(i).regressors.code];
           
        modreg = modelData(i).models(j).regressors;
        
        [q,rgind] = ismember(modreg,rgcodes);
        
        mrg =  regData(i).regressors(rgind);
        rp = pool(mrg);
        
        spb = sparseblock(ones(1,rp.Npar),[mrg.Npar]);
        member = unsparsify(cumsum(spb,2).*spb);
        
        
        [usep,srtind] = sort(abs(u(:,sep)),'descend');
        
        FS(i).model(j).sepmat = rgind(rp.info.parentindex(srtind));
        FS(i).model(j).usep = usep;
        FS(i).model(j).member = member(rp.info.parentindex(srtind)); 
        
    end
    
end

 

