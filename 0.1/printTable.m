function printTable(model,texfname)

% Prints a text file summarizing model fitting results

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

getmodels = 10;
% getmodels =[ 1     2     1     3     2     2     3]
if nargin < 2
    
    fnstr = 'stats_summary_%i.txt';
    i = 1;
    while exist(sprintf(fnstr,i),'file')
        i = i+1;
    end
    
    fname = sprintf(fnstr,i);
end

% if isnumeric(fname)
%     fid = fname;
% else
% %     fid = fopen(fname,'w');
% end


rlbl = model(1).models(getmodels(1)).regLabels;
%     for i = 1:length(model(getmodel).models)
%        if isempty(model(j).models(getmodel).fit)
%             error('Model does not appear to have been fit to data yet ("fit" field is empty).')
%        end
        
       
%          fprintf(fid,'\n\nModel: %s',model(1).models(getmodel).label);

%          fprintf(fid,'\nRegressors:%s',sprintf('\n\t%s',rlbl{:}));
%         for j = 1:length(model)
%            if ~isequal(model(j).models(getmodel).regLabels,rlbl), error('Models do not appear to be the same across subjects');end
% 
%              fprintf(fid,'\n\nFull:\n\tLL: %i\n\tNpar: %i\n\tAIC: %i\n\tBIC: %i',model(j).models(getmodel).fit(1).LL,...
%              model(j).models(getmodel).fit(1).npar,model(j).models(getmodel).fit(1).AIC,model(j).models(getmodel).fit(1).BIC);
% 
%             par = zeros(1,2*length(model(j).models(getmodel).fit(1).parest));
% 
%              par(1:2:end) = model(j).models(getmodel).fit(1).parest;
%              par(2:2:end) = sqrt(diag(model(j).models(getmodel).fit(1).I^-1));
%             fprintf(fid,'\nParameter Est: ');
%             fprintf(fid,'\t%1.2g(%1.2g)',par);
%         end


%          fprintf(fid,'\n\nSubmodels:');
            V = {};
         for q = 1: length(model(1).models(getmodels(1)).fit)
            fname = regexprep(sprintf('%s',model(1).models(getmodels(1)).fit(q).label),'_','\\_');
            Header = {' ','-2*LLR','','Pvalue','df','$\Delta$AIC','$\Delta$BIC'};
                
            V(end+1,:) = mat2cell(nan(size(Header)),1,ones(1,length(Header)));
            V{end,1} = fname;
            for j = 1:length(model)        
                    if length(getmodels)>1
                        getmodel = getmodels(j);
                    else
                        getmodel = getmodels;
                    end
%                 V = cat(1,V,{nan,-2*diff([model(j).models(getmodel).fit([1 q]).LL]),[],...
%                     model(j).models(getmodel).fit(q).llrpval, diff([model(j).models(getmodel).fit([q 1 ]).npar]),diff([model(j).models(getmodel).fit([q 1]).AIC]), diff([model(j).models(getmodel).fit([q 1]).BIC])});
                if q == 1
                    bln = 0;
                    LL0 = model(j).models(getmodel).fit(1).LL-abs(model(j).models(getmodel).fit(1).LLR);
                    dAIC=model(j).models(getmodel).fit(1).AIC+2*LL0;
                    dBIC=model(j).models(getmodel).fit(1).BIC+2*LL0;
                    model(j).models(getmodel).fit(q).LLR = -abs(model(j).models(getmodel).fit(q).LLR);
                else
                    bln = model(j).models(getmodel).fit(1).npar;
                    dAIC=diff([model(j).models(getmodel).fit([q 1]).AIC]);
                    dBIC = diff([model(j).models(getmodel).fit([q 1]).BIC]);
                end
                V = cat(1,V,{nan,-2*model(j).models(getmodel).fit(q).LLR,[],...
                    model(j).models(getmodel).fit(q).llrpval, model(j).models(getmodel).fit(q).npar-bln, dAIC, dBIC});
                
            end
            
%             eval(sprintf('!latex "%s"',fname))
%             eval(sprintf('!latex2html ''%s''',fname))
%             eval(sprintf('!dvipdf "%s"',fname))
            
    %         par = zeros(1,2*length(model(i).fit(j).parest));
    %         
    %         fprintf(fid,'\t\t(%i)',2*diff([model(i).fit([1 j]).LL]),...
    %                 -model(i).fit(j).llrpval, diff([model(i).fit([j 1 ]).npar]),diff([model(i).fit([j 1]).AIC]), diff([model(i).fit([j 1]).BIC]));
    % %         fprintf(fid,'\t\t\tLLR: %2.2f\t\tPvalue: %2.2f\t\tdf: %i\t\tAIC: %2.2f\t\tBIC: %2.2f',2*[model(i).fit([j 1]).LL],...
    %                 model(i).fit(j).llrpval, diff([model(i).fit([1 j ]).npar]),model(i).fit(j).AIC, model(i).fit(j).BIC);
         end
            ps = [V{:,4}];
            V(ps>.05,3 ) = strcat({'    '},V(ps>.05 ,3));
            V(ps<.05 & ps>.01,3) = strcat({'*   '},V(ps<.05 & ps>.01,3));
            V(ps<.01 & ps>.001,3) = strcat({'**  '},V(ps<.01 & ps>.001,3));
            V(ps<.001,3) = strcat({'*** '},V(ps<.001,3));
            
            MakeTable(V,'head',Header','filename',texfname,'fontsize','small'); 

%          fprintf(fid,'\n\n-----------------\n');
%     end
% % end
% if ~isnumeric(fname)
%     fclose(fid);
% end