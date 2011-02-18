function varargout = appendReg(h,R2,dataset,varargin)

% appendReg(h,R,dataset)
%
% h - GazeReader handle.
% R - array of regressors.
%
% Appends a regressor to regData.regressors for the specified data
% set,  checking for consistency of the design matrix size and noptions 
% field. If dataset is not specified, the currently active set is used
% If noptions is empty, then it is assigned from an existing regressor.
% Otherwise, if the size of R or noptions is not consistent with existing 
% regressors, an error is returned.
% 
% Rout = appendReg(R1,R2,dataset)
% 
% If the first argument is a regressor structure rather than the GazeReader
% handle, R2 is appended to R1 and returned in Rout.
%
% appendReg automatically assigns a unique code to the appended regressor.
%
%
% appendReg(h,X,dataset,...)
%
% Creates a regressor from the matrix X and appends it to regressors in regData
% (or to R1). All arguments after data set are passed to makeRegressor.
%
% X must either have one row for each bin in for each observation, or it must
% have length equal to the number of observations or the number of trials.
% In the latter cases element elements in X are assigned to every bin within the the
% corresponding interval.
%
% See also MAKEREGRESSOR

% C Kovach 2011



if ishandle(h)
   
    if nargin < 3 || isempty(dataset)
        dataset = getappdata(h,'CurrentDataSet');
    end

    regData = getappdata(h,'regData');
    R1 = regData(dataset).regressors;

    set_regData = true;
    trialData = getappdata(h,'trialData');
    ntr = length(trialData(dataset).trials);


elseif isstruct(h)
    
    R1 = h;
    
    set_regData = false;

    ntr = 0;
else
    
    error('First argument is not recognized as a handle or regressor structure.')
    
end

if isnumeric(R2)
    
%     R2 = makeRegressor(R2,varargin{:});
    R2 = buildpolyreg(R2,1,'subtractdc',false,varargin{:});

end

noptions = R1(1).noptions;      % Number of bins for each observation
nobs = length(R1(1).noptions); % Number of observations
nrows = size(R1(1).value,1);   % Number of rows in the design matrix



nrows2 = size(R2.value,1);

noptions = R1(1).noptions;
nbin = unique(noptions);
if length(nbin) > 1
    nbin = nan;
end

if ishandle(h)
    trialnbins = [trialData(dataset).trials.nbin];
else
    trialnbins = -1;
end

%%% for numeric inputs, create regressors of the appropriate size
switch nrows2
    case  sum( noptions )
    % OK
    case length( noptions )  
        % If the input matches the number of observations, in is applied to
        % all bins on a given observation
        R2.value = R2.value(regData(dataset).fixationIndex,:);
    
    case ntr
    % If it matches the number of trials, it is applied to all bins within a given trials
        
        R2.value = R2.value(regData(dataset).trialIndex,:);        
        
    case nbin
        % If number of bins is constant input equals bin size, then input
        % is assumed to apply to each bin. 
        
        R2.value = repmat(R2.value,length(noptions),1);        

    case sum(trialnbins)
        
        R2.value = R2.value(regData(dataset).binIndex(:,1));
        
    
    otherwise %oops
                
        error(sprintf('The number of rows (=%i) is not consistent with noptions.\nInput should have %i, %i or %i rows.',nrows2, nrows, length(noptions), ntr))
    
    
end


%%% Check consistency of regressors: before doing anything else we will
%%% verify that the regressors are consistent with each other.


%%%%%% Check that regressors are consistent in their dimensions %%%%%
for i = 2:length(R1)
    
    if length(R1(i).noptions) ~= nobs || any(R1(i).noptions ~= noptions)
        error('noptions field is not consistent for all elements of the first regressor group.')
    end
    
    if size(R1(i).value,1) ~= nrows
        error('Number of rows is not consistent for all elements of the first regressor group.')
    end
    
end

%%%  Check that codes are unique, assign new codes if not.
if length(unique([R1.code])) ~= length(R1)
    warning('Codes for the first group are not unique. They will be reassigned.')
    for i = 1:length(R1)
        R1(i).factmat(R1(i).factmat==R1(i),code) = i;
        R1(i).code = i;
        R1(i).codevec(:) = i;
        R1(i).codevec(:) = i;
        
    end
end

for i = 1:length(R1)
    
    if any(R1(i).noptions ~= noptions)
        error('noptions field is not consistent for all elements of R1.')
    end
    
    if size(R1(i).value,1) ~= nrows
        error('Number of rows is not consistent for all elements of R1.')
    end
end

codeincr = R1(end).code;
codes1 = [R1.code];
for i = 1:length(R2)
    
    if ~isempty(R2(i).noptions) && any(R2(i).noptions ~=noptions)        
        error('noptions for R2 is non-empty but doesn''t match that for R1.')
    end
    
    if size(R2(i).value,1) ~= nrows
        error('Number of rows in the design for R2 doesn''t match that for R1.')
    end
    if    isempty(R2(i).code) || ismember(R2(i).code,[codes1,R2(1:i-1).code])
%         updatefact = (R2(i).factmat==R2(i).code);
%         if isempty(updatefact), updatefact = ones(size(R2.levmat));end
        R2(i).factmat(R2(i).factmat==R2(i).code) = codeincr+i;
        R2(i).info.functionInputCodes(1,R2(i).info.functionInputCodes(1,:)==R2(i).code) = codeincr+i;
        R2(i).code = codeincr + i;         
    end
    
    R2(i).noptions = noptions;
    R2(i).codevec = ones(1,R2.Npar)*R2(i).code;
%     R2(i).factmat = ones(1,R2.Npar)*R2(i).code;
end



Rout = [R1,R2];

if set_regData
    
    regData(dataset).regressors = Rout;
    regData(dataset).codeincr = max([Rout.code]);
   
    % reassign the appended regressors
    setappdata(h,'regData',regData)
end

if nargout > 0
    varargout{1} = Rout;
end



    
    





    
    