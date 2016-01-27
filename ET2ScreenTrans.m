
function ET2SCRN = ET2ScreenTrans(eydfile, calibPts,varargin)

%function ET2SCRN = ET2ScreenTrans(eydfile, calibPts)
% Returns an affine transformation matrix to convert between raw eyetracker
% coordinates and screen coordinates. CalibPts can be an Nx2 list of
% calibration points in the order used by the eye tracker. If CalibPts
% is a 2X2 vector then it is assumed to represent fixed x and y spacing (top row)
% and screen center (bottom row) for a rectangular 9 point calibration with 
% ordering top left to bottom right.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------




flipy = false;

i = 1;
while i <= length(varargin)
   switch lower(varargin{i})
       case 'flipy' %flips about the mean for the y of the calibration points 
            flipy  = true;
            
       otherwise

           error([varargin{i},' is not a valid option.']);
   end         
   i = i+1;
end


if size(calibPts,1) == 2
    calibPts = [kron(calibPts(1,1)*[1 1 1]',[-1 0 1]') + calibPts(2,1),kron(calibPts(1,2)*[-1 0 1]',[1 1 1]')+ calibPts(2,2)];
end
    
    
eydDat = ReadEYD(eydfile,'headerOnly');

for i = 1:size(calibPts,1)
    if isfield(eydDat.allHeaderData,'Calibration_Values')
        ETCalibPts(i,1) = str2double(eydDat.allHeaderData.Calibration_Values.(sprintf('htgt_data_%i',i)));
        ETCalibPts(i,2) = str2double(eydDat.allHeaderData.Calibration_Values.(sprintf('vtgt_data_%i',i)));
    elseif isfield(eydDat,'ETsettings')
        ETCalibPts(i,1) = str2double(eydDat.ETsettings(1).m_eteCalibration.target_points.Point(i).X);
        ETCalibPts(i,2) = str2double(eydDat.ETsettings(1).m_eteCalibration.target_points.Point(i).Y);
    else
        error('Unable to find calibration data.')
    end
        
        
end

if flipy
    calibPts(:,2) = 2*mean(calibPts(:,2))-calibPts(:,2);
end

scrncal = cat(2,calibPts,ones(size(calibPts,1),1));
etcal = cat(2,ETCalibPts,ones(size(ETCalibPts,1),1));
%Transformation matrix from eyetracker to screen pixel coordinates
%ET2SCRN = (calibPts'*etcal)*(etcal'*etcal)^-1; 
ET2SCRN = (scrncal'*etcal)*(etcal'*etcal)^-1; 
