function reseed

%
% RESEED seeds the uniform and normal random generators, RAND and RANDN
% with the current time obtained with NOW.
%
% See also: RAND, RANDN, NOW.

umethod = 'twister';
nmethod = 'state';

s = typecast(now,'uint32');
rand(umethod, double( s(1) )); %Using the least significant 4 bytes of the value returned by now

randn( nmethod, double( swapbytes( s(1) ) ) ); %using swapbytes so that the randn and rand seeds are
                                             %at least somewhat arbitrary
                                             %with respect to each other.
                                             


