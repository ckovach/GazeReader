      
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

% C Kovach 2008

intxnord = 1;
center = true;
noptions = [];
labels = repmat({[]},1,size(F,2));
cellmeans = false;
codeincr = 0;
postmult = 1;
i = 1;

ignore = false(size(F));

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
      case 'finis'
         
      otherwise
            error([varargin{i},' is not a valid keyword.'])
    end
    i = i+1;
end 

if ~iscell(labels)
    labels = {labels};
end

for f = 1:size(F,2)

%     ufs = unique(F(:,f));
      ufs = unique(F(~ignore,f));

    XF = zeros(size(F,1),length(ufs));

    for l = 1:length(ufs)

        XF(:,l) = F(:,f) == ufs(l);

    end


    if center
        %remove mean
        u = ones(1,size(XF,1))'./sqrt(size(XF,1));
        XF = XF-u*(u'*XF);
        XF = XF(:,2:end);
    end

    RF(f) = makeregressor(XF,'noptions',noptions,'label',labels{f},'codeincr',codeincr,'postmultiply',postmult);

end

if intxnord>1
    RF = cat(2,RF,interaction(RF,'intxnord',intxnord));
end

