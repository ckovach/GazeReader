      
function RF  = fact2reg(F,varargin)

%Builds Regressor from a factor matrix.
%Each level in each factor generates a regressor
%factlev is a matrix of factors (top row) and levels (bottom row) associated
% with corresponding columns of XF. contrasts is the contrast matrix needed
% to compute F statistics for each factor.

% C Kovah 2008

intxnord = 1;
center = true;
noptions = [];
labels = repmat({[]},1,size(F,2));
cellmeans = false;
codeincr = [];
i = 1;

varargin{end+1} = 'finis';
while i <= length(varargin)
    switch lower(varargin{i})
       case 'intxnord'
            intxnord = varargin{i+1};
            i = i+1;
       case 'fullintxn'
            intxnord = size(F,2);
       case 'center'
            center = varargin{i+1};
            i = i+1;
       case 'noptions'
            noptions = varargin{i+1};
            i = i+1;
       case 'labels'
            labels = varargin{i+1};
            i = i+1;
        case 'cellmeans'
            cellmeans = varargin{i+1};
            i = i+1;
        case 'codeincr'
            codeincr = varargin{i+1};
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

    ufs = unique(F(:,f));

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

    RF(f) = makeregressor(XF,'noptions',noptions,'label',labels{f},'codeincr',codeincr);

end

if intxnord>1
    RF = cat(2,RF,interaction(RF,'intxnord',intxnord));
end

