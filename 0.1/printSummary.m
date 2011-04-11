function printSummary(model,fname)

% Prints a text file summarizing model fitting results

if nargin < 2
    
    fnstr = 'stats_summary_%i.txt';
    i = 1;
    while exist(sprintf(fnstr,i),'file')
        i = i+1;
    end
    
    fname = sprintf(fnstr,i);
end

if isnumeric(fname)
    fid = fname;
else
    fid = fopen(fname,'w');
end

if isempty(model.fit)
    error('Model does not appear to have been fit to data yet ("fit" field is empty).')
end

for i = 1:length(model)
  
     fprintf(fid,'\n\nModel: %s',model(i).label);
      
     fprintf(fid,'\nRegressors:%s',sprintf('\n\t%s',model(i).regLabels{:}));
     fprintf(fid,'\n\nFull:\n\tLL: %i\n\tNpar: %i\n\tAIC: %i\n\tBIC: %i',model(i).fit(1).LL,...
         model(i).fit(1).npar,model(i).fit(1).AIC,model(i).fit(1).BIC);

     par = zeros(1,2*length(model(i).fit(1).parest));
     
     par(1:2:end) = model(i).fit(1).parest;
     par(2:2:end) = sqrt(diag(model(i).fit(1).I^-1));
    fprintf(fid,'\nParameter Est: ');
    fprintf(fid,'\t%1.2g(%1.2g)',par);
     
     

     fprintf(fid,'\n\nSubmodels:');
     fprintf(fid,'\n\t\t\tLLR\t\t\tPvalue\t\t\tdf\t\tDelta AIC\t\t\tDelta BIC');
     for j = 2: length(model.fit)
        fprintf(fid,'\n\n%s:\n',model(i).fit(j).label);
        fprintf(fid,'\t\t%2.2f\t\t%7.2g\t\t\t%i\t\t%8.4g\t\t%.4g',-2*diff([model(i).fit([1 j]).LL]),...
                model(i).fit(j).llrpval, diff([model(i).fit([j 1 ]).npar]),diff([model(i).fit([j 1]).AIC]), diff([model(i).fit([j 1]).BIC]));
%         par = zeros(1,2*length(model(i).fit(j).parest));
%         
%         fprintf(fid,'\t\t(%i)',2*diff([model(i).fit([1 j]).LL]),...
%                 -model(i).fit(j).llrpval, diff([model(i).fit([j 1 ]).npar]),diff([model(i).fit([j 1]).AIC]), diff([model(i).fit([j 1]).BIC]));
% %         fprintf(fid,'\t\t\tLLR: %2.2f\t\tPvalue: %2.2f\t\tdf: %i\t\tAIC: %2.2f\t\tBIC: %2.2f',2*[model(i).fit([j 1]).LL],...
%                 model(i).fit(j).llrpval, diff([model(i).fit([1 j ]).npar]),model(i).fit(j).AIC, model(i).fit(j).BIC);
     end
     
     fprintf(fid,'\n\n-----------------\n');
end
    
if ~isnumeric(fname)
    fclose(fid);
end