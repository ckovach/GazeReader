
function imstruct = makeImageStruct(pths,imgs)

imstruct = struct('filename',[],'path',[],'info',[],'position',[],'screenres',[],'xySize',[],'code',[],'activeTrials',[]);

if nargin == 0
    return
elseif nargin == 1
    imgs = pths;
    pths = '';
end

if ~iscell(imgs), imgs = {imgs};end
if ~iscell(pths), pths = {pths}; end
if length(pths) == 1
    pths = repmat(pths,1,length(imgs));
end

imstruct = repmat(imstruct,1,length(imgs));

multassign = @ (varargin) varargin{:};
    
[imstruct.filename] = multassign(imgs{:});
[imstruct.path] = multassign(pths{:});
