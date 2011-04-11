function FX = readFIX(fname,varargin)

%Reads data from an ASL FIX file

% modified 6/7 to allow for error in the sequential position of the xdat signal


i = 1;
oldstyle = 0;
bitmask = 255;  %By default, masking xdat with 255
pulse = true;  %Xdat values are assumed to be pulses rather than steady values
nhist = 1;      %number of points in xdat history to include in xdhist (exlcuding the current xdat value)
headerOnly=false;
segment = [];
while i <= length(varargin)
   switch lower(varargin{i})
       case 'bitmask'
            bitmask = varargin{i+1};
            i = i+1;
       case 'pulse'
            pulse = varargin{i+1};
            i = i+1;
       case 'impulse'
            impulse = 1;
       case 'oldstyle'
            oldstyle = 1;
       case 'headeronly'
           headerOnly=true;
       case 'segment'
           segment =varargin{i+1};
            i = i+1;
       otherwise

           error([varargin{i},' is not a valid option.']);
   end         
   i = i+1;
end



FX = struct('fileInfo',struct,'allHeaderData',struct,'seg',struct([]));
FX.seg = struct('info',struct,'allHeaderData',struct,'fix',struct([]),'xdat',struct([]),'etfs',[],'xdatCodes',[],'summary',[]);

if nargin < 1 | isempty(fname)
    [fname,fpath] = uigetfile({'*.fix','ASL FIX file'});
else
    [fpath,fname,ext] = fileparts(fname);
    fname = [fname,ext];
end

if isempty(fpath)
    fpath = pwd;
end

fid = fopen(deblank(fullfile(fpath,fname)));

FX.filename = deblank(fullfile(fpath,fname));

if fid <= 0 
    error('Failed to open file')
end

% formatLabels = {'VT_UI1','VT_UI2','VT_I1','VT_I2'};
% precisnMap = {'uint8','uint16','int8','int16'};


FX.allHeaderData = struct;
FX.fileInfo.name = fname;
FX.fileInfo.path = fpath;
RecDat = [];
xdhist = zeros(1,nhist+1);

line = '';
segnum = 1;
xdt = -Inf;
oldseg = 0;
while ~feof(fid) && ~headerOnly


    while isempty(line) && ~feof(fid)
        line = fgetl(fid);
    end
    q = regexp(line,'\[[\w\s-]*\]','match');
    line = regexprep(q{end},'[\s&()\.-]','_');

    currBranchLbl = line(2:end-1);


    switch regexprep(line,'\d*','#')

        case '[Fixation_File_Header]'
            line = fgetl(fid);
            while ~feof(fid) & ~isempty(strtrim(line))  
                  [a,b] = strtok(line,':=');                  
                  a = strtrim(a); b = strtrim(b(2:end));
                  a = regexprep(a,'[\s&()]','_');
                  FX.fileInfo.(a) = b;              


                  FX.allHeaderData.(currBranchLbl).(a) = b;
                  line = fgetl(fid);
            end
        case '[Fix_Segment_Header_#]'
            temp = regexp(line,'\d*','match');
            segnum = str2double(temp{1});
            line = fgetl(fid);
            if ~isempty(segment) && ~ismember(segnum,segment)
                segnum = 0;
%                 continue
            elseif ismember(segnum,segment)
                FX.seg(segnum == segment).segment = segnum;              
                segnum = find(segnum == segment);
            end

            while ~feof(fid) & ~isempty(strtrim(line))  
                  [a,b] = strtok(line,':=');                  
                  a = strtrim(a); b = strtrim(b(2:end));
                  a = regexprep(a,'[\s&()]','_');
                  if segnum > 0
                    FX.seg(segnum).info.(a) = b;              
                    FX.seg(segnum).allHeaderData.(currBranchLbl).(a) = b;
                  end
                  line = fgetl(fid);
            end


        case '[Fix_Segment_Summary_#]'
            temp = regexp(line,'\d*','match');
            segnum = str2double(temp{1});
            line = fgetl(fid);
            
            if ~isempty(segment) && ~ismember(segnum,segment)
                segnum = 0;
%                 continue
            elseif ismember(segnum,segment)

                FX.seg(segnum == segment).segment = segnum;              
                segnum = find(segnum == segment);
                
            end
            
            while ~feof(fid) & ~isempty(strtrim(line))  
                  [a,b] = strtok(line,':=');                  
                  a = strtrim(a); b = strtrim(b(2:end));
                  a = regexprep(a,'[\s&()]','_');
                  if segnum > 0
                    FX.seg(segnum).summary.(a) = b;              
                  end
                  line = fgetl(fid);
            end



        case '[Fix_Data_#]'

             temp = regexp(line,'\d*','match');
             segnum = str2double(temp{1});
             line = fgetl(fid);
           
            if ~isempty(segment) && ~ismember(segnum,segment)
                segnum = 0;
%                 continue
            elseif ismember(segnum,segment)
                FX.seg(segnum == segment).segment = segnum;              
                segnum = find(segnum == segment);
                
            end

             temp = regexp(line,'\w*','match');

        %      FX.seg(segnum).colLabel = temp(2:end);
             line = fgetl(fid);

             t0 = 24*3600*datenum(regexp(line,'\d*:\d*:[\d.]*','match'));
             xdat = 0;

             fn = 1;
             xdn = 1;
             while ~feof(fid) & ~isempty(strtrim(line))  


                sep = regexp(line,'[\w\.:]*','match');
                t = 24*3600*datenum(regexp(line,'\d*:\d*:[\d.]*','match')) - t0;

                if strcmp(sep{2},'x')
                    xdtemp = str2double(sep{3});
                    if ~isempty(bitmask)
                            xdtemp = bitand(xdtemp,bitmask);
                    end
                    if ~pulse || xdtemp ~= 0
                        xdat = xdtemp;
                        xdhist(1:end-1) = xdhist(2:end);
                        xdhist(end) = xdat;
%                         xdatT = t;
                        if segnum > 0
                            FX.seg(segnum).xdat(xdn).startT = t;
                            FX.seg(segnum).xdat(xdn).id = xdat;
                            FX.seg(segnum).xdat(xdn).code = num2str(xdat);
                        end
                        
                        xdn = xdn+1;
                        xdt = t;
                    end
                    
                elseif segnum > 0
                    
                    FX.seg(segnum).fix(fn).meanPos = [str2double(sep{9}),str2double(sep{10})];
                    FX.seg(segnum).fix(fn).startT = t;
                    FX.seg(segnum).fix(fn).num = str2double(sep{3});
                    FX.seg(segnum).fix(fn).dur = str2double(sep{6});
                    FX.seg(segnum).fix(fn).endT = t + str2double(sep{6});
%                   if oldstyle
                        if t > xdt
                            FX.seg(segnum).fix(fn).xdat = xdat;
                        end
%                   end                    
                    FX.seg(segnum).fix(fn).xdhist = xdhist;
                    FX.seg(segnum).fix(fn).xdindex = xdn-1;
                    FX.seg(segnum).fix(fn).interfixt = str2double(sep{7}); %time between the end of the last fixation and start of the 
                                                                        % current one, as reported by ASL
                    FX.seg(segnum).fix(fn).pupil = str2double(sep{11});
                    FX.seg(segnum).fix(fn).eye = str2double(sep{12});
                    FX.seg(segnum).fix(fn).interfixdeg = str2double(sep{8});
                    FX.seg(segnum).fix(fn).pln = str2double(sep{4});
                    FX.seg(segnum).fix(fn).scn = str2double(sep{13});
                    FX.seg(segnum).fix(fn).flags = str2double(sep{14});
                    FX.seg(segnum).fix(fn).loss = str2double(sep{15});
                    if fn>1
                        FX.seg(segnum).fix(fn).shiftvec = diff(cat(1,FX.seg(segnum).fix(fn + [-1:0]).meanPos));
                        FX.seg(segnum).fix(fn).dt = diff(cat(1,FX.seg(segnum).fix(fn + [-1:0]).startT)); % difference between fixation onsets of 
                                                                                                         %the current and previous fixation
                    else
                        FX.seg(segnum).fix(fn).shiftvec =[nan nan];
                         FX.seg(segnum).fix(fn).dt =nan;
                   end
                        
                    fn = fn+1;
                end
                 line = fgetl(fid);

             end     


                 line = fgetl(fid);



        otherwise

              line = fgetl(fid);  
             while ~feof(fid) & ~isempty(strtrim(line)) & line(1) ~= '['  
                  [a,b] = strtok(line,': =');
                  a = strtrim(a); b = strtrim(b(2:end));
                  a = regexprep(a,'[\s&()\.-]','_');

                  if ~isempty(regexp(a(1),'[1-9]'))
                      a = ['x',a];
                  end


                 FX.allHeaderData.(currBranchLbl).(a) = b;
                 line = fgetl(fid);
             end

    end
    if segnum > 0
        if isfield(FX.fileInfo,'exact_data_rate')   
            FX.seg(segnum).etfs = str2double(FX.fileInfo.exact_data_rate);
        elseif isfield(FX.fileInfo,'data_rate')
            FX.seg(segnum).etfs = str2double(FX.fileInfo.data_rate);
        else
            FX.seg(segnum).etfs = str2double( FX.allHeaderData.File_Description.Update_Rate_Hz_);
        end    
    end
    
    
    if segnum > 0 && segnum ~= oldseg
        if ~isempty(FX.seg(segnum).xdat)
            xdts = [FX.seg(segnum).xdat.startT];
        end

        
        if ~oldstyle && ~isempty(FX.seg(segnum).fix) && ~isempty(FX.seg(segnum).xdat)   

            for fn = 1:length(fxts)
                xdi = max(find(fxts(fn)>xdts));
                if ~isempty(xdi)
                    currxdat = FX.seg(segnum).xdat( xdi ).id;
                else
                    currxdat = 0;
                end
                FX.seg(segnum).fix(fn).xdat = currxdat;
            end

        end
        oldseg = segnum;
        
    end

end
if ~isempty(FX.seg(segnum).fix)
    fxts = [FX.seg(segnum).fix.startT];
    xds = [FX.seg(segnum).fix.xdat];
    FX.seg(segnum).xdatCodes = unique(xds(xds~=0));
else
    FX.seg(segnum).xdatCodes = [];
end


FX.units = 'xy Eyetracker';
fclose(fid);
