function RWD = concatRWD(lastendts,varargin)

%RWDcat = concatRWD(RWD1,RWD2,RWD3,...)
%Concatenates multiple rawGazeData structures into a single one.


if length(varargin)>1
    RWDin = cat(2,varargin{:});
else
    RWDin = varargin{1};
end


lastendt = 0;

fields = {'horz','vert','pupil','xdat'};
for i = 1:length(RWDin)
        
        if i == 1
            RWD= RWDin(i);
        else
            
            
            for k = 1:length(fields)
             RWD.(fields{k}) = cat(2,RWD.(fields{k}),RWDin(i).(fields{k}));
            end
             if isfield(RWD,'degConversion')
                 RWD.degConversion = cat(1,RWD.degConversion,RWDin(i).degConversion);
             end
             RWD.time = cat(2,RWD.time,RWDin(i).time+double(lastendt));
                                
          
        end

        lastendt = lastendts(i);
end




        