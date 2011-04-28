function binData  = importIAS(fnames,codeincr,varargin)

% Function to import an eyelink IAS (region of interest) file and convert
% into bin groups.
%
% Usage:
%
%  binData  = importIAS(filename(s),codeincr)
%           
%       Arguments are filenames and codeincr value.
%
%  binData  = importIAS(filename(s),gr_handle)
%
%        Appends bin data to existing instance of gazeReader.
%
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


screendim = [1 1];
i = 1;
while i <= length(varargin)
   switch lower(varargin{i})
       case 'type'
          type = varargin{i+1};
          i = i+1;
       case 'screendim'
          screendim = varargin{i+1};
          i = i+1;
       otherwise
           error([varargin{i},' is not a valid option.']);
   end         
   i = i+1;
end

    
           

if ischar(fnames)
    fnames = {fnames};
end

if nargin < 2
    codeincr = 0;
end

if ishandle(codeincr) && isstruct(getappdata(codeincr,'expEventData')) 
    h = codeincr;
    binData = getappdata(h,'binData');
    if isempty(binData)
        binData = makeBinData([]);
    end
    codeincr = binData.codeincr;
else
    binData = makeBinData([],'codeincr',codeincr);
    h = nan;
end
    
    
for i = 1:length(fnames)
    
    fid = fopen(fnames{i},'r');
    
    while ~feof(fid)
        
        line = fgets(fid);
        
        type = regexp(line,'[A-Z]*','match','once');
        
        
        switch type
            case 'RECTANGLE'
                xydata = regexp(line,'(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*(\w*)','tokens','once');
                lbl = xydata{end};
                xydata = cellfun(@str2num,xydata,'uniformoutput',false);
                data = [xydata{[2 4 3 5]}]./screendim([1 1 2 2]);
                bintype = 'rect';
                
%                 iasdata = regexp(line,'(?<id>\d+)\s+(?<left>\d+)\s+(?<top>\d+)\s+(?<right>\d+)\s+(?<bottom>\d+)\s+(?<label>\w*','tokens');
                id = xydata{1};
            case 'ELLIPSE'
%                 iasdata = regexp(line,'(?<id>\d+)\s+(?<left>\d+)\s+(?<top>\d+)\s+(?<right>\d+)\s+(?<bottom>\d+)\s+(?<label>\w*)','tokens');
                xydata = regexp(line,'(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*(\w*)','tokens','once');
                lbl = xydata{end};
                xydata = cellfun(@str2num,xydata,'uniformoutput',false);
                xydata = [xydata{2:end}];
                eldat = xydata([ 1 2 3 4])./screendim([1 2 1 2]);
                data = [mean(eldat([1 3;2 4]')), diff(eldat([1 3;2 4]'))/2];
                bintype = 'ellipse';
            case 'FREEHAND'
                id = regexp(line,'(\d+)','tokens','once');
                coord = regexp(line,'(\d+,\d+)','tokens');
                coord = cellfun(@str2num,[coord{:}],'uniformoutput',false);
                data = cat(1,coord{:})*diag(screendim.^-1);
                lbl = regexp(line,'(\w+)\s*$','tokens','stringanchors');
                lbl = lbl{1}{1};
                bintype = 'poly';
            otherwise
                continue
        end
        
        [pth,fn] = fileparts(fnames{i});
        label =  sprintf('%s: %s',fn,lbl);
        binData = makeBinData(binData,data,'type',bintype,'label',label,'codeincr',codeincr);
        
                
    end        
end

if ishandle(h) && isstruct(getappdata(h,'expEventData'))   
    setBinData(h,binData);
end

