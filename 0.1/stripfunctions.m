
function R = stripfunctions(R)

% Removes all function handles from a structure array, descending through all fields. 
% One possible reason for wanting to do this is that saving a function handle to a .mat file 
% entails saving the entire workspace visible to the function, resulting in
% unnecessarily large files.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

%C. Kovach 2011


if isempty(R)
    return
end


if ~isstruct(R)
    error('Input must be a structure or structure array, or else empty.')
end

fieldns = fieldnames(R);


for k = 1:length(fieldns)
    
    fn = fieldns{k};
    if isa(R(1).(fn),'function_handle')  %If field contains a function handle, remove it
    
        for i = 1:length(R)
            R(i).(fn) = [];  
%             R(i).(fn)(1) = [];          % Retain the field as a zero-length function handle
                                        % to avoid possible inconsistencies later.
        end
        
    elseif isa( R(1).(fn) ,'struct')  %If field contains a structure run stripfunctions recursiveley
        
        for i = 1:length(R)
            R(i).(fn) = stripfunctions( R(i).(fn) );            
        end
        
    end 
    
end



