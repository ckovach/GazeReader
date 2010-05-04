
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

fY = fft(Y);
zeta = abs(fY).^2;

toeind = toeplitz(uint16(1:N));

w = [0:floor(N./2),-ceil(N./2-1):-1];

% boxord = 5;

    
A = specstart;
Anew = A;
LL = -Inf;

Q=[];
lls = [];
damp = 0;
dampadd = .5;
while derr > tol
%     hold on, 
%     plot(real(ifft(A)))
%     shg
%     pause
    
    rho = real(ifft(Anew));
    
    LLnew = rho*Y' - sum(exp(rho));

    
    if LLnew > LL
        
        A = Anew;
        lls(end+1) = LL;
        damp = damp+dampadd ;
        LL = LLnew;

    else
        damp = damp-dampadd;
    end
    
    fL = fft( exp ( rho ) );
    
%     Q(:,end+1) = fL;
    
    DL = fY - fL;

    %The second Derivative matrix is the circular convolution matrix (Toeplitz) of fL
    
%     D2L = fL(toeind);
    
%     del = DL*D2L^-1;
    del = DL./fL(1);
    
    del(abs(w) > boxord) = 0;
    Anew = A + del*N*exp(damp);
    
    Anew(2:end) = (Anew(2:end) + fliplr(conj(Anew(2:end))))./2;
    
    derr = abs(sum(del.^2))
    
end
    
    
% keyboard    
    
    