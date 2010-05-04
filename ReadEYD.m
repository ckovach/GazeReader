function ASLdat = ReadEYD( fname, varargin)

%function ASLdat = ReadEYD( fname );
%Reads data from data file into a matlab structure.
%
%function ASLdat = ReadEYD( fname , segment );
%Reads data from specified segment.
%

%Modified 10/21 to handle xml parts

headerOnly = false;
segment = 0;

i = 1;
while i <= length(varargin)
   switch lower(varargin{i})
       case 'headeronly'
            headerOnly  = true;
       case 'segment'
            segment  = varargin{i+1};
            i = i+1;
       otherwise

           error([varargin{i},' is not a valid option.']);
   end         
   i = i+1;
end

xdatMask = 255; %Mask of 255 assumes that only the first 8 bits of the xdat are used in event markers.

ASLdat = struct('fileInfo',struct,'allHeaderData',struct,'seg',struct);



if nargin < 1 || isempty(fname)
    [fname,fpath] = uigetfile({'*.eyd','ASL data file'});
else
    [fpath,fname,ext] = fileparts(fname);
    fname = [fname,ext];
end
           
           
if isempty(fpath)
    fpath = pwd;
end

fid = fopen(deblank(fullfile(fpath,fname)));

ASLdat.filename = deblank(fullfile(fpath,fname));

if fid <= 0 
    error('Failed to open file')
end

formatLabels = {'VT_UI1','VT_UI2','VT_I1','VT_I2','byte','uint16','int8','int16'};
precisnMap = {'uint8','uint16','int8','int16','uint8','uint16','int8','int16'};


ASLdat.allHeaderData = struct;
ASLdat.fileInfo.name = fname;
ASLdat.fileInfo.path = fpath;
RecDat = [];

line = '';

getXML = false;

while ~strcmp(line,'[Segment_Data]')

    line = deblank(line);
    while isempty(line) && ~feof(fid)
        line = deblank(fgetl(fid));
    end
    
    % Newer versions contain xml segment -- get this and save it to another
    % file for parsing later
    isxml = regexp(line,'.xml]','once');
    if ~isempty(isxml)
        getXML = true;
        fxml = fopen('temp.xml','w');
        
        while isempty(regexp(line,'</Ete>','once'))
            line = fgetl(fid);
            fprintf(fxml,'%s\n',line);
        end
        fclose(fxml);
        line = fgetl(fid);
        continue
    end
    
    q = regexp(line,'\[[\w\s-]*\]','match');
    line = regexprep(q{end},'[\s&()\.-]','_');
                 
    currBranchLbl = line(2:end-1);
    
    switch line

        case '[File_Description]'
            line = fgetl(fid);
            while ~isempty(strtrim(line))  
                  [a,b] = strtok(line,':');                  
                  a = strtrim(a); b = strtrim(b(2:end));
                  a = regexprep(a,'[\s&()]','_');
                  ASLdat.fileInfo.(a) = b;              
                  
                  
                  ASLdat.allHeaderData.(currBranchLbl).(a) = b;
                  line = fgetl(fid);
            end
        case '[System_Data_Items]'
            line = fgetl(fid);
            line = fgetl(fid);
            i = 1;
            while ~isempty(strtrim(line))  
                  
                  [lbl,pos,fmt,nb] = strread(line,'%s%f%s%f');
                  RecDat(end+1).lbl = lbl{1};
                  RecDat(end).pos = pos;
                  RecDat(end).fmt = precisnMap{ strcmpi( fmt{1},formatLabels)};
                  RecDat(end).nbytes = nb;
                  RecDat(end).scalef = 1;
                  i = i+1;
                  
                  ASLdat.allHeaderData.(currBranchLbl).(lbl{1}) = sprintf('%f\t%s\t%f',pos,fmt{1},nb);
                  line = fgetl(fid);
            end
        case '[User_Description]'
            
             line = fgetl(fid);
            while ~isempty(strtrim(line))  
                  [a,b] = strtok(line,':');                  
                  a = strtrim(a); b = strtrim(b(2:end));
                  a = regexprep(a,'[\s&()]','_');
                 
                  ASLdat.allHeaderData.(currBranchLbl).(a) = b;
                  line = fgetl(fid);
            end
            
        case '[Data_Items_Selected_by_User]'
            line = fgetl(fid);
            line = fgetl(fid);
            i = 1;
            while ~isempty(strtrim(line))  
                  
                  [lbl,pos,fmt,nb,sf] = strread(line,'%s%f%s%f%f');
                  RecDat(end+1).lbl = lbl{1};
                  RecDat(end).pos = pos;
                  RecDat(end).fmt = precisnMap{ strcmpi( fmt{1},formatLabels)};
                  RecDat(end).nbytes = nb;
                  RecDat(end).scalef = sf;
                  i = i+1;
                  
                  ASLdat.allHeaderData.(currBranchLbl).(lbl{1}) = sprintf('%f\t%s\t%f\t%f',pos,fmt{1},nb,sf);
                  line = fgetl(fid);
            end
        
        case '[Segment_Data]'

        otherwise
                          
              line = fgetl(fid);  
             while ~isempty(strtrim(line)) && line(1) ~= '['  
                  [a,b] = strtok(line,': =');
                  a = strtrim(a); b = strtrim(b(2:end));
                  a = regexprep(a,'[\s&()\.-]','_');
                  
                  if ~isempty(regexp(a(1),'[1-9]'))
                      a = ['x',a];
                  end
                  
                  
                 ASLdat.allHeaderData.(currBranchLbl).(a) = b;
                 line = fgetl(fid);
              end
            
    end
end        

for i = 1:length(RecDat)
    
    ASLdat.fileInfo.RecordFormat.(RecDat(i).lbl) = RecDat(i);
    
end

if ~headerOnly
    
    nsegments = str2double(ASLdat.allHeaderData.Segment_Variables.User_Recorded_Segments);
    segDirAddr = str2double(ASLdat.allHeaderData.Segment_Variables.Segment_Directory_Start_Address);

    if ASLdat.allHeaderData.File_Description.File_Version(1) == '5'
        asl5000Format = true;
    elseif ASLdat.allHeaderData.File_Description.File_Version(1) == '6'

        asl5000Format = false;
    else
        error('Unrecognized Version')
    end



    recByteN = [RecDat.nbytes];
    recBytePos = [RecDat.pos];
    BytesPerRec = sum(recByteN);

    startPos = ftell(fid);

    NRec = (segDirAddr - startPos)./BytesPerRec ;


    segPosDat = struct; %Segment positions data
    fseek(fid,segDirAddr,-1);
    s = 1;
    nPseudoSegs = 0;
    tmp = 0;
    if asl5000Format 
            cases = [255,254,253]; %Segment, Pseudosegment, and end of segment data
            nget = [3,2,2];        %Number of bytes for each type 
            offset = 0;
    else
            cases = [251,-1,252]; %I don't know indicator of pseudosegment for asl6000...
            nget = [7,2,2];
            offset = 1;
    end

    while tmp ~= cases(3)

        tmp = fread(fid,1,'*uint8');

        switch tmp

            case cases(1)
                inp = fread(fid,nget(1),'*uint32');
                segPosDat(s).startByte = inp(1+offset);
                segPosDat(s).startTime = inp(2+offset);
                segPosDat(s).endTime = inp(3+offset);
                if ~asl5000Format 
                    segPosDat(s).nRec =inp(end);
                end
            case cases(2)
                fread(fid,nget(2),'*uint32');
                nPseudoSegs = nPseudoSegs+1; 
            case cases(3)
                lastByte = fread(fid,nget(3),'*uint32');

            otherwise 
                error('Error reading segment metadata!');
        end
        s = s+1;
    end

    SegNRec = diff([segPosDat(:).startByte,segDirAddr])./BytesPerRec ;

    if asl5000Format 
        for i = 1:length(SegNRec)
            segPosDat(i).nRec = SegNRec(i);
        end
    end


    if segment == 0
        segment = 1:length(segPosDat);
    end

    if max(segment) > length(segPosDat) 
        error('There aren''t that many segments!')
    end

    ASLdat.seg = struct('horz',[],'vert',[],'pupil',[],'xdat',[]);
    reclbls = {RecDat(:).lbl};
    Ncycles = length(segPosDat)*length(RecDat);
    % Ncycles = length(segPosDat);
    cyc = 0;
    h =waitbar(cyc/Ncycles,sprintf('Reading EYD file %s ...',fname));
    
    etfs = str2num(ASLdat.fileInfo.Update_Rate_Hz_);
    for s = 1:length(segment) 

        for i = 1:length(RecDat)
            fseek(fid,double(recBytePos(i) - 1 + segPosDat(segment(s)).startByte),-1);

            prec = RecDat(i).fmt;

            ASLdat.seg(s).otherDat.(RecDat(i).lbl) = fread(fid,double(segPosDat(segment(s)).nRec ),['*',prec],BytesPerRec - RecDat(i).nbytes );

            cyc = cyc +1;
            waitbar(cyc/Ncycles,h,sprintf('Reading segment %i of %i ...',s,length(segPosDat)));
            drawnow  
        end

        %Correcting lag introduced by overtimes
        overtime = ASLdat.seg(s).otherDat.overtime_count;
        ovtCorrecn = cumsum(double(overtime));
        corrtimes = 0;
        corrtimes([1:length(overtime)]' + ovtCorrecn(:)) = 1;
        corrtimes = cumsum(corrtimes); %Indices adjsuted to reflect lost data points in overtime 
                                       %Lost points are replaced with the last good point. 

        ASLdat.seg(s).horz = RecDat(strcmp(reclbls,'horz_gaze_coord')).scalef .* double(ASLdat.seg(s).otherDat.horz_gaze_coord(corrtimes))';
        ASLdat.seg(s).vert = RecDat(strcmp(reclbls,'vert_gaze_coord')).scalef .* double(ASLdat.seg(s).otherDat.vert_gaze_coord(corrtimes))';
        ASLdat.seg(s).pupil = RecDat(strcmp(reclbls,'pupil_diam')).scalef .* double(ASLdat.seg(s).otherDat.pupil_diam(corrtimes))';
        ASLdat.seg(s).xdat = bitand( double( ASLdat.seg(s).otherDat.XDAT(corrtimes) ) , xdatMask )';
        ASLdat.seg(s).otherDat = rmfield(ASLdat.seg(s).otherDat,{'XDAT','vert_gaze_coord','pupil_diam','start_of_record'});
        ASLdat.seg(s).otherDat.otCorrecter = corrtimes;  % Fields in otherDat have to be corrected for overtime 
                                                              % with ..otherDat.fieldname(ASLdat.seg(s).otherDat.otCorrecter ) 
        ASLdat.seg(s).time =  [1:length( ASLdat.seg(s).horz)]./etfs;
        
        %Extract event markers from the xdat column
        ASLdat.seg(s).EventMark.T = find(diff([0, ASLdat.seg(s).xdat]) > 0);
        if any(ASLdat.seg(s).EventMark.T == 1)
            warning('\nXDAT signal may be erroneous at %i in segment %i', find(diff(ASLdat.seg(s).EventMark.T) == 1),s)
        end
        ASLdat.seg(s).EventMark.id = ASLdat.seg(s).xdat( ASLdat.seg(s).EventMark.T);


    end


end
    
if getXML
    fprintf('Parsing XML...');
    
    xmlstruct = parseXML('temp.xml');
    ASLdat.ETsettings = xmlstruct.Ete.m_eteEye.EteEye;
    ASLdat.xml_all = xmlstruct.Ete;
    fprintf('Done\n');
    
end
    
fclose(fid);

if ~headerOnly
    delete(h)  
end

    
            


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following functions are copied and modified
% from the xmlread help in matlab 7.4.0.287
%
function theStruct = parseXML(filename)
% PARSEXML Convert XML file to a MATLAB structure.
try
   tree = xmlread(filename);
catch
   error('Failed to read XML file %s.',filename);
end

% Recurse over child nodes. This could run into problems 
% with very deeply nested trees.
try
   theStruct = parseChildNodes(tree);
catch
   error('Unable to parse XML file %s.',filename);
end


% ----- Subfunction PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
   childNodes = theNode.getChildNodes;
   numChildNodes = childNodes.getLength;
   allocCell = cell(1, numChildNodes);

%    children = struct(             ...
%       'Name', allocCell, 'Attributes', allocCell,    ...
%       'Data', allocCell, 'Children', allocCell);

       
      
        for count = 1:numChildNodes
            ch = childNodes.item(count-1);
            nname = char(ch.getNodeName);
            if ~strcmp(nname,'#text')
                if isfield(children,nname)
                    
                    children.(nname) = cat(2,children.(nname), parseChildNodes(ch));
                else
                    children.(nname) = parseChildNodes(ch);
                end
            elseif numChildNodes == 1 && any(strcmp(methods(ch), 'getData'))            
               children = char(ch.getData); 
            end                
            
         end
%         theChild = childNodes.item(count-1);
%         children(count) = makeStructFromNode(theChild);
 end



% ----- Subfunction MAKESTRUCTFROMNODE -----
% function nodeStruct = makeStructFromNode(theNode)
% % Create structure of node info.
% 
% nodeStruct = struct(                        ...
%    'Name', char(theNode.getNodeName),       ...
%    'Attributes', parseAttributes(theNode),  ...
%    'Data', '',                              ...
%    'Children', parseChildNodes(theNode));
% 
% if any(strcmp(methods(theNode), 'getData'))
%    nodeStruct.Data = char(theNode.getData); 
% else
%    nodeStruct.Data = '';
% end
% 
% % ----- Subfunction PARSEATTRIBUTES -----
% function attributes = parseAttributes(theNode)
% % Create attributes structure.
% 
% attributes = [];
% if theNode.hasAttributes
%    theAttributes = theNode.getAttributes;
%    numAttributes = theAttributes.getLength;
%    allocCell = cell(1, numAttributes);
%    attributes = struct('Name', allocCell, 'Value', ...
%                        allocCell);
% 
%    for count = 1:numAttributes
%       attrib = theAttributes.item(count-1);
%       attributes(count).Name = char(attrib.getName);
%       attributes(count).Value = char(attrib.getValue);
%    end
% end
% 
% 
% 
