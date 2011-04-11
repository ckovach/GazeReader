function CheckFix(ASL,trials)

bl = 1;
s =1 ;
trialidfield = 'xdat';
durfield = 'dur';
et2screenfield = 'et2scrn';
% trialidfield = 'trial';
% durfield = 'durT';

%durfield = 'durT';
showCorrection = true;


if nargin < 2
   trials =  1:length(ASL(s).segData(bl).arrayFixRoi)
end

circw = 0:.1:2*pi;
rad = 200; %Radius of circle in pixels per second

if showCorrection & isfield(ASL(s).segData.arrayFixData,'trendP')
    trendP = ASL(s).segData.arrayFixData.trendP    
else
    trendP = [0 0]';
end
%Check fixations by plotting them
              
    etdatH =ASL(s).segData(bl).EYDdata.horz;
    etdatV =ASL(s).segData(bl).EYDdata.vert;
    otcorr = cumsum(double(ASL(s).segData(bl).EYDdata.otherDat.overtime_count)); %OVertime correction
    
    screenDat = (ASL(s).segData(bl).arrayFixData.(et2screenfield )*cat(2,etdatH,etdatV,ones(length(etdatH),1))')';
 %   ArrayStartTs(ASL(s).segData(bl).trialData.Array.id) = ASL(s).segData(bl).trialData.Array.T; 
     ArrayStartTs = ASL(s).segData(bl).arrayStartTs; 
    f1 = figure;
    a1 = gca;
    a2 = axes;
    %f2 = figure;
    %ArrayStartTs = ASL.segData(bl).arrayFixRoi.startT;  
    %ArrayStartTs = [ASL(s).segData(bl).trialData.Array.T];
    etfs = ASL(s).segData.arrayFixData.etfs
    if isfield(ASL(s),'tunits') & strcmp(ASL(s).tunits,'sec')
%         tnorm = str2num(ASL(s).fixfileinfo.exact_data_rate);
%         etfs = tnorm;
        tnorm = etfs;
    else
        tnorm = 1;
%         etfs = 60;
    end
    
    t = [1:length(etdatH)]./tnorm;
    trend = [polyval(trendP(1,:),t); polyval(trendP(2,:),t);polyval(0,t)]';
    screenDatDetrend = screenDat-trend;
    
    stimFileOrder = ASL(s).respData.stimset
    
    for i = trials

        if ~isempty(ASL(s).segData(bl).arrayFixRoi(i).roi)

%            X = imread(sprintf(stimPicNames{s},i));
            
            X = imread(sprintf(ASL(s).stimPicNames{ASL(s).respData(bl).stimset},i));
    %         subplot(2,1,1)
            %figure(f1)
            axes(a1)
            X = X(end:-1:1,:,:);
            imagesc(X),axis('xy'), shading flat, colormap gray, axis image
            axis([0  1024 0 768])
            tpoints = round([0:4*etfs] + ASL(s).segData(bl).arrayStartTs(i)*tnorm);


            hold on 

            af = ASL(s).segData(bl).arrayFixData.fix([ASL(s).segData(bl).arrayFixData.fix.(trialidfield)] == i );
            trd = ASL(s).segData(bl).arrayFixRoi(i);
            drawt = tpoints(1);
            for k = 1:length(af);

    %             subplot(2,1,1)
%               figure(f1)
              axes(a1)
              tsub = drawt:round(tnorm*(af(k).startT + af(k).(durfield)));
              tsub = tsub + otcorr(tsub)'; 
                if ~any(trend(:))
                    plot(screenDat(tsub,1),screenDat(tsub,2))
                else
                    plot(screenDat(tsub,1),screenDat(tsub,2),':')
                    plot(screenDatDetrend(tsub,1),screenDatDetrend(tsub,2))
                end    %             subplot(2,1,2), hold on
%                 figure(f2), hold on
                  axes(a2), hold on
                plot(tpoints./etfs - ASL(s).segData(bl).arrayStartTs(i),[etdatH( tpoints),etdatV( tpoints),ones(length( tpoints),1)]*ASL(s).segData(bl).arrayFixData.(et2screenfield )(1:2,:)');

                drawt = round(tnorm*(af(k).startT + af(k).(durfield)));

                tstring = sprintf('Trial %i\tFixation %i\tEm: %i \t Inv: %i \t Targ: %i \t',i,k,trd.em(k),(trd.roi(k) > 3) - 2*(trd.roi(k) == 0) ,((trd.roi(k) <= 3) == ASL(s).respData.targupright) - 2*(trd.em(k) == 0 ) );
                %if af(k).suspect == 1, tstring = [tstring,sprintf('\tSUSPECT')];end

                fxtimepts = round(tnorm*([1./etfs:1./etfs:af(k).(durfield)]+af(k).startT));
               fxtimepts = fxtimepts + otcorr(fxtimepts)';
    %             subplot(2,1,1)
%                 figure(f1)
                axes(a1)
                if ~any(trend(:))
                    p= plot(screenDat(fxtimepts,1),screenDat(fxtimepts,2),'r');
                else
                    dur = af(k).dur;
                    p(1)= plot(screenDat(fxtimepts,1),screenDat(fxtimepts,2),'r:');
                    p(2)= plot(screenDatDetrend(fxtimepts,1),screenDatDetrend(fxtimepts,2),'r');
                    p(3)= plot(rad*dur*cos(circw)+ af(k).posdetrend(1),rad*dur*sin(circw) + af(k).posdetrend(2),'r');
                    p(4)= plot(af(k).posdetrend(1),af(k).posdetrend(2),'r+');
                end
    %             subplot(2,1,2)
%                 figure(f2)
                axes(a2)
               
                p2= plot(double(fxtimepts)./etfs  - ASL(s).segData(bl).arrayStartTs(i) ,[etdatH( fxtimepts),etdatV( fxtimepts),ones(length( fxtimepts),1)]*ASL(s).segData.arrayFixData.(et2screenfield )(1:2,:)','r');
                axis([0 4 0 1024])
                set(a1,'position',[.05 .25 .9 .7])
                set(a2,'position',[.05 .05 .9 .15])
                
                drawnow, shg
%             figure(f1)
            axes(a1)
                title(tstring)
                pause
                delete(p(1:2))
                delete(p2)
            end
    %         subplot(2,1,2)
%             figure(f2),
            axes(a2)
            cla
        end
    end
% close(f1)
% close(f2)

