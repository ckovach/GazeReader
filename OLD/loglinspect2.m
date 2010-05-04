
function LL = loglinspect(N,T,boxord)

% Computes the log linear spectrum with likelihood maximization

% sinord = 100; %Number of points in the 
if nargin == 1 || isempty(T)
    T = find(N);
    N = length(N);
end

if max(T) > N
    error('Maximum event time must be less than ore equal to N')
end



specstart = zeros(1,N);

Y = zeros(1,N);

Y(T) = 1;

tol = 1e-6;
derr = Inf;

Tdiffs = -repmat(T',1,length(T)) + repmat(T,length(T),1);
Tdiffs = Tdiffs(find(Tdiffs > 0)); % Only past events are considered

unqTdiffs = unique(Tdiffs);

TdiffsN = sum( repmat( Tdiffs, 1 , length(unqTdiffs) ) == repmat(unqTdiffs',length(Tdiffs),1));

DY = zeros(1,N);

DY(unqTdiffs) = TdiffsN;




fY = fft(Y);
fDY = fft(DY);

zeta = abs(fY).^2;

toeind = toeplitz(uint16(1:N));

w = [0:floor(N./2),-ceil(N./2-1):-1];

if nargin < 3
    boxord = 20;
end
    
A = specstart;
Anew = A;
alpha = 0;
alphanew = 0;

LL = -Inf;

Q=[];
lls = [];
damp = 0;
dampadd = .5;
while derr > tol
   
  
    
    
    frho = Anew .* conj(fY);  %Convolution of gamma and observed vector
    
    frho(1) = alphanew;
    rho = real(ifft(frho));
    
  
    LLnew = rho*Y' - sum(exp(rho));

    
    if LLnew >= LL
        
        A = Anew;
        alpha = alphanew;
        lls(end+1) = LL;
        damp = damp+dampadd ;
        LL = LLnew;

        hold on, 
        plot(rho)
        shg
        pause

    else
        damp = damp-dampadd;
    end
    
    fL = fft( exp ( rho ) );
    
%     Q(:,end+1) = fL;
    
    DL = zeta - fL.*conj(fY);

    dalpha = (sum(Y) - fL(1))./(fL(1));

    
    %The second Derivative matrix is the circular convolution matrix (Toeplitz) of fL
    
%     D2L = fL(toeind);
    
%     del = DL*D2L^-1;
    del = DL./(fL(1)*abs(fY(1)).^2);
    
    del(abs(w) > boxord) = 0;
    Anew = A + del*N*exp(damp);
    
      gnew = ifft(Anew);
    gnew(1:floor(N./2)) = 0;    
    Anew = fft(gnew);
    
    alphanew = alpha + dalpha*N*exp(damp);
    
%     Anew(2:end) = (Anew(2:end) + fliplr(conj(Anew(2:end))))./2;
    
    derr = abs(sum((Anew-A).^2) + dalpha^2)
    
end
    
    
keyboard   
    
    