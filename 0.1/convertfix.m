function FIX = convertfix(FIX, et2scrn);

%converts data in FIX to screen units from eyetracker and vice versa

default_starting_units = 'xy Eyetracker';
if nargin < 2
    et2scrn = FIX.et2scrn;
end

for i = 1:length(FIX)
    
        if strcmp(FIX(i).units,'xy Eyetracker');
            E = et2scrn;
            FIX(i).units = 'xy Screen';
        elseif strcmp(FIX(i).units,'xy Screen');
    %         E = et2scrn^-1;
    %         FIX(i).units = 'xy Eyetracker';
            E = eye(size(et2scrn));
        end

        if isfield(FIX(i),'units')
            units = FIX(i).units;
        else
            units = default_starting_units;
        end

      for s = 1:length(FIX(i).seg)



             for j = 1:length(FIX(i).seg(s).fix)         
                 newpos = [FIX(i).seg(s).fix(j).meanPos, 1]*E';
                 FIX(i).seg(s).fix(j).meanPos = newpos(1:2);
                 newshift = [FIX(i).seg(s).fix(j).shiftvec, 0]*E';
                 FIX(i).seg(s).fix(j).shiftvec = newshift(1:2);
                 if isfield(FIX(i).seg(s).fix(j),'posdetrend')
                     newpos = [FIX(i).seg(s).fix(j).posdetrend, 1]*E';
                     FIX(i).seg(s).fix(j).posdetrend = newpos(1:2);
                 end
             end

    end

end