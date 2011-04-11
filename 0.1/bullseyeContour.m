function c = bullseyeContour(x,y,z, rs, ths, varargin)

%Plot polar patch data with a supermiposed target


plotargs = {'k'};
textargs = {};
i = 1;
fig = [];
textth = ths;
while i <= length(varargin)
   switch lower(varargin{i})
       case 'plotargs'
            plotargs = varargin{i+1};
            i = i+1;
       case 'textargs'
            textargs = varargin{i+1};
            i = i+1;
       case 'figure'
            fig = varargin{i+1};
            i = i+1;
       case 'textth'
            textth = varargin{i+1};
            i = i+1;
            
       otherwise
            break
%            error([varargin{i},' is not a valid option.']);
   end         
   i = i+1;
end

if nargin < 6
    varargin = {'k'};
end

s = [0:.01:2*pi]';

if isempty(fig)
    fig = figure;
else
    figure(fig)
end
[gz,gx,gy] = patch2surf(x,y,z);

c = contour(gx,gy,gz,varargin{i:end});

shading interp
axis image

hold on

rlbl = {};
rloc = [];
dr = mean(diff(rs));
for r = rs 
   plot(r*cos(s),r*sin(s),plotargs{:})
   rlbl{end+1} = sprintf('%1.2f',r);
   
   rloc(:,end+1) = [0 (r+.1*dr)]; 
end

thloc = [0 0]';
thlbl = {''};
for th = ths
    
   plot(r*sin(th)*[0 1],  r*cos(th)*[0 1],plotargs{:})
  
end

for th = textth
    if abs(th) > .01 & abs(th - 2*pi) > .01
       thlbl{end+1} =   sprintf('%1.2f',th*180./pi);
       thloc(:,end+1) = [r*sin(th),  r*cos(th)];
    end
end

text(rloc(1,:),rloc(2,:), rlbl,textargs{:});
 
text(thloc(1,:),thloc(2,:), thlbl,textargs{:});
   
   