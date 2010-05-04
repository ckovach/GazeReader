function [f,frq] = makefourier(nfrq,varargin)

% [f,frq] = makefourier(nfrq,varargin)
%
%   f is a handle to a function which returns the value of each sinusoid at
%   the specified point.  nfrq is a vector containing the number frequency 
%     components to retain for each dimension.

%C Kovach 2007


varargin{end+1} = 'finis';

circular = 1;
keepdc = 0;
freqconst = 2*pi; %Fits to the circle (toroid). 2*pi assumes function to be continuous across opposite ends. 
i = 1;
useindex = 1:length(nfrq); %indices of input to use in computing the function value
while i <= length(varargin)
    switch lower(varargin{i})
       case 'circular' %Filters circular region of frequency space
            circular = 1;
       case 'square' %Filters square region of frequency space
            circular = 0;            
       case 'keepdc' 
            keepdc = 1;            
       case 'freqconst' %set freqconst to some other frequency unit, eg pi (default is 2*pi)
            freqconst = varargin{i+1};
            i = i+1;
       case 'useindex' 
            useindex  = varargin{i+1};
            i = i+1;            
       case 'finis'
       otherwise
            error([varargin{i},' is not a valid keyword.'])
    end
    i = i+1;
end


frq = [-nfrq(1):nfrq(1)]';

for i = 2:length(nfrq)
    partfrq = [-nfrq(i):nfrq(i)]';
    frq= cat(2,kron(frq,ones(length(partfrq ),1)),kron(ones(length(frq),1),partfrq));
end

%rows that are the additive inverse of another row are redundant
[redundant,redi] = ismember(frq,-frq,'rows');
redmat = spalloc(size(frq,1),size(frq,1),size(frq,1));
redmat((find(redi)-1)*size(frq,1) + redi(find(redi))) = 1;

%vector of redundant frequency components to discard
discfrq = sum(tril(redmat))' > 0;


%Discard frequency components that lie outside of ellipse bounded by elements of nfreq
%in each dimension
if circular
    frqdist = sqrt(sum((frq./nfrq(ones(size(frq,1),1),:)).^2 ,2)) ;
    discfrq = discfrq  | frqdist > 1;
end

frq = frq(~discfrq,:);
frqdist = frqdist(~discfrq,:);

%Sort rows in terms of ascending frequencies
[srtf,srti] = sort(frqdist);
frq= frq(srti,:);
frq = kron(frq,[1 1]');
ph = zeros(size(frq,1),1);
ph(2:2:end) = -pi./2;
if keepdc
   frq(2:end+1,:) = frq;
   frq(1,:) = 0;
   
   frqdist(2:end+1) = frqdist;
   frqdist(1,:) = 0;
   if isempty(ph)
       ph = 0;
   else
        ph(2:end+1) = ph;
   end
end

%This is the function handle
%f = @(txy,th) (th*[exp(j*freqconst*frq*txy(:,useindex)')])';
f = @(txy,th) (th*[cos(freqconst*frq*txy(:,useindex)' + ph(:,ones(size(txy,1),1)))]);


