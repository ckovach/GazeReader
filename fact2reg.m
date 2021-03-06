      
function RF  = fact2reg(F,varargin)

% function RF  = fact2reg(F)
%
% Builds Regressor from a factor matrix, F.
%
% Each level in each factor generates a column in the regressor R.value
% which contains a dummy variable for the corresponding factor x level.
% RF is a regressor structure.
%
% Note that the columns of RF.value can be mapped onto the corresponding 
% factor x level interactions based on the colums of RF.factmat (factor) and
% RF.levmat (level). Which give the factors and levels for all terms in the
% interaction.                     
%
%
% function RF  = fact2reg(F,'option',value)
%
%  Create RF with the following options:
%         'intxnord'      -      Specify order of interaction (default  =1)
%         'fullintxn'     -      model with all interaction orders (takes
%                                 no argument)
%         'center'        -     Center regressors if true (default == true)
%                                   and decorrelated with the intercept.
%         'noptions'
%         'ignore'        - logical index into F: ignores entries in F for 
%                             which ignore is true. fact2reg does not create
%                             an indictor variable for unique entries
%                             contained only in ignored rows.
%
%   For remaining options see also MAKEREGRESSOR 

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

% C Kovach 2008

intxnord = 1;
center = true;
noptions = [];
labels = repmat({[]},1,size(F,2));
% cellmeans = false;
codeincr = 0;
postmult = 1;
i = 1;
obsfreq = ones(size(F,1),1);

ignore = false(size(F,1),1);

varargin{end+1} = 'finis';
while i <= length(varargin)
    switch lower(varargin{i})
       case 'intxnord'                  %order of interaction
            intxnord = varargin{i+1};
            i = i+1;
       case 'fullintxn'                 %mode with all interaction orders
            intxnord = size(F,2);
       case 'center'
            center = varargin{i+1};
            i = i+1;
       case 'noptions'
            noptions = varargin{i+1};
            i = i+1;
       case {'labels','label'}
            labels = varargin{i+1};
            i = i+1;
%         case 'cellmeans'
%             cellmeans = varargin{i+1};
%             i = i+1;
        case 'codeincr'
            codeincr = varargin{i+1};
            i = i+1;
        case 'ignore'       % logical index into F: ignores entries in F for 
                            % which ignore is true. fact2reg does not create
                            % an indictor variable for unique entries
                            % contained only in ignored rows.
            ignore = varargin{i+1};
            i = i+1;
        case 'postmultiply'
            postmult = varargin{i+1};
            i = i+1;
        case 'obsfreq'
            obsfreq = varargin{i+1};
            i = i+1;
      case 'finis'
         
      otherwise
            error([varargin{i},' is not a valid keyword.'])
    end
    i = i+1;
end 

if ~iscell(labels)
    labels = {labels};
end

if iscell(F)  %Allow for cell arrays as well as numeric vectors
    eqfun = @(a,b) cellfun(@(a) isequal(a, b{1}),a);
else
    eqfun = @(a,b) a==b;
end

for f = 1:size(F,2)

%     ufs = unique(F(:,f));
      ufs = unique(F(~ignore,f));  % Get unique values in the factor

    XF = zeros(size(F,1),length(ufs));
    
    for l = 1:length(ufs)
        
        XF(:,l) = eqfun(F(:,f) , ufs(l) ).*~ignore;

    end


    if center
        %remove mean
        if ~all(obsfreq==1)
    %         u = ones(1,size(XF,1))'./sqrt(size(XF,1));
            spbl = sparseblock(ones(1,size(XF,1)),noptions.*ones(size(obsfreq)));
        else 
            spbl = 1; 
        end
        
        u1 = (spbl'*obsfreq).*~ignore;
        u1 = u1./sum(u1);
        u2 = ~ignore;
        m = u2*(u1'*XF);
        XF = XF-m;
        XF = XF(:,2:end);
        fun= mkfun(ufs,m(find(~ignore,1),:));
    else
        fun= mkfun(ufs);
    end

    RF(f) = makeregressor(XF,'noptions',noptions,'label',labels{f},'codeincr',codeincr,'postmultiply',postmult); %#ok<*AGROW>
    RF(f).info.factorlabels = ufs;
    RF(f).function =fun;
end

if intxnord>1
    RF = cat(2,RF,interaction(RF,'intxnord',intxnord));
end

%%%%%%%%%%%
function fn = mkfun(ufs,m)

%%% The function is nested here to avoid having to save the workspace to disk

    if nargin > 1 %center
        fn = @(F)mkcfun(F,ufs,m);
    else
        fn = @(F)mkcfun(F,ufs);
        
    end


%%%%%%%%%
function XF = mkcfun(F,ufs,m)

if iscell(F)  %Allow for cell arrays as well as numeric vectors
    eqfun = @(a,b) cellfun(@(a) isequal(a, b{1}),a);
else
    eqfun = @(a,b) a==b;
end
XF = zeros(size(F,1),length(ufs));

for l = 1:length(ufs)

    XF(:,l) = eqfun(F , ufs(l) );

end
if nargin > 2
    XF = XF-repmat(m,size(XF,1),1);
    XF = XF(:,2:end);
end



