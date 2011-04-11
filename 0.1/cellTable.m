function Vcell = cellTable(model, getmodels,subjlabel)

% Prints a text file summarizing model fitting results
%

if nargin < 3
    dosubjlabel = 0;
    subjlabel = 0;
end


if nargin < 2 || getmodels == 0;
    getmodels = 1:length(model(1).models);
end

if length(subjlabel) ==1 && subjlabel
    subjlabel = num2cell(1:length(model(1).models(1).fit));
    dosubjlabel = 1;
elseif ~isequal(subjlabel,0) 
    dosubjlabel = 1;
end

normalized=0;
if nargin < 2
    
    fnstr = 'stats_summary_%i.txt';
    i = 1;
    while exist(sprintf(fnstr,i),'file')
        i = i+1;
    end
    
    fname = sprintf(fnstr,i);
end



rlbl = model(1).models(getmodels(1)).regLabels;


            V = {};
         for q = 1: length(model(1).models(getmodels(1)).fit)
            fname = regexprep(sprintf('%s',model(1).models(getmodels(1)).fit(q).label),'_','\\_');
%             Header = {' ','-2*LLR','','Pvalue','df','norm. $\rho^2$','$\Delta$AIC','$\Delta$BIC'};
            Header = {' ','-2*LLR','','Pvalue','df','Pseudo-R^2','dAIC','dBIC'};
%             Header = {' ','-2*LLR','','Pvalue','df','Pseudo-R$^2$','$\Delta$AIC','$\Delta$BIC'};
                
            V(end+1,:) = mat2cell(nan(size(Header)),1,ones(1,length(Header)));
            V{end,1} = fname;
            for j = 1:length(model)        
                    if length(getmodels)>1
                        getmodel = getmodels(j);
                    else
                        getmodel = getmodels;
                    end
                if q == 1
                    bln = 2*model(j).models(getmodel).fit(1).npar;
                    LL0 = model(j).models(getmodel).fit(1).LL-abs(model(j).models(getmodel).fit(1).LLR);
                    dAIC=model(j).models(getmodel).fit(1).AIC+2*LL0;
                    dBIC=model(j).models(getmodel).fit(1).BIC+2*LL0;
                    model(j).models(getmodel).fit(q).LLR = -abs(model(j).models(getmodel).fit(q).LLR);
                else
                    bln = model(j).models(getmodel).fit(1).npar;
                    dAIC=diff([model(j).models(getmodel).fit([q 1]).AIC]);
                    dBIC = diff([model(j).models(getmodel).fit([q 1]).BIC]);
                end
                
                
                if dosubjlabel 
                    lbl = subjlabel{j};
                else
                    lbl = nan;
                end
                
                if ~normalized
                V = cat(1,V,{lbl,-2*model(j).models(getmodel).fit(q).LLR,[],...
                    model(j).models(getmodel).fit(q).llrpval, bln-model(j).models(getmodel).fit(q).npar,model(j).models(getmodel).fit(q).R, dAIC, dBIC});
                else
                    V = cat(1,V,{lbl,-2*model(j).models(getmodel).fit(q).LLR,[],...
                    model(j).models(getmodel).fit(q).llrpval, bln-model(j).models(getmodel).fit(q).npar,model(j).models(getmodel).fit(q).R./model(j).models(getmodel).fit(1).R, dAIC, dBIC});
                end
            end
            
    end
ps = [V{:,4}];
V(ps>.05,3 ) = strcat({'    '},V(ps>.05 ,3));
V(ps<.05 & ps>.01,3) = strcat({'*   '},V(ps<.05 & ps>.01,3));
V(ps<.01 & ps>.001,3) = strcat({'**  '},V(ps<.01 & ps>.001,3));
V(ps<.001,3) = strcat({'*** '},V(ps<.001,3));


if nargout > 0
    Vcell = cat(1,Header(:)',V);
    Vcell(cellfun(@(X)any(isnan(X)),Vcell(:))) = {''};
end

