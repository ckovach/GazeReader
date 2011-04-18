function RD = concatRD(lastfxns,lasttrs, varargin)

%RDcat = concatRD(nfixations,RD1,RD2,RD3,...)
%Concatenates multiple regData structures into a single one. This assumes
%that the corresponding regressors are the same in each case.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

RD = makeRegData([]);
trind = [];
lasttr = 0;
fxind = [];
lastfx = 0;
add_indicator = true; %Adds a regressor that indicates the original subject

if length(varargin)>1
    RDin = cat(2,varargin{:});
else
    RDin = varargin{1};
end

indr = [];
for i = 1:length(RDin)
        if isfield(RDin(i),'trialIndex')
            trind = cat(1,trind(:),RDin(i).trialIndex + lasttr);
            lasttr = lasttr + lasttrs(i);
            RD.trialIndex = trind;
        end
        if isfield(RDin(i),'fixationIndex')
            fxind = cat(1,fxind(:),RDin(i).fixationIndex + lastfx);
            lastfx = lastfx + lastfxns(i);
            RD.fixationIndex = fxind;
        end

        if i == 1
            RD= RDin(i);
        else


             for k = 1:length(RDin(i).regressors)   

                 if ~isequal(RD.regressors(k).label,RDin(i).regressors(k).label)
                     warning('Labels do not match for all regressors. Are they the same quantities?')
                 end

                 RD.regressors(k).value = cat(1,RD.regressors(k).value,RDin(i).regressors(k).value);
                 RD.regressors(k).noptions = cat(1,RD.regressors(k).noptions(:),RDin(i).regressors(k).noptions(:));
             end
        end
        if    add_indicator
            indr = cat(1,indr,ones(size(RDin(i).regressors(1).value,1),1)*i);
        end
end

if   add_indicator
    RD.regressors = appendReg(RD.regressors,indr,[],'label','data_set');
end