
function B = chooseperm(N, k, item, rows)

% function  B = chooseperm(N, k)
%    Returns indices for every unordered combination of k items from a
%    population of N using a recursive algorithm.


if nargin < 3
    rows = 1;
    item = 1;
end

B = [];

if item >= k
    
    B = (rows:N)';
    
else
    
    for strow = rows:N

       b = chooseperm( N, k, item+1, strow+1);
    
       if ~isempty(b)
           
           B = cat(1,B,cat(2,ones(size(b,1),1)*strow,b));
           
       else
           
           return
           
       end
    end    
end

