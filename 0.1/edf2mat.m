function edf2mat(dirfile)

% edf2mat
%
%     Extract data from all edf files using readEDF and save to mat file. 
%
% edf2mat(dirname) 
%
%     Extract all files in specified directory.
%
% edf2mat(edffile)
%
%     Extract data in specified edf file. 
%
% See also READEDF



% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin <1 
    dirfile = cd;
end

if ~isempty(regexp(dirfile,'\.edf','once'));
    edfs = {dirfile};
else
    currdir = cd;
    cd(dirfile);
%     edfs = dir(spritnf('%s%s%s',dirfile,filesep,'*.edf'));
    edfs = dir('*.edf');

    edfs = {edfs.name};

end

err = {};

for i = 1:length(edfs)
    savefname = sprintf('edf_%s',regexprep(edfs{i},'.edf','.mat'));
   
    fprintf('\nConverting %s\n',savefname)
    
    if exist(savefname)
        fprintf('%s already exists in this folder. Skipping...\n',savefname)
        continue;
    end
    
    try
        [FIX,RAW] = readEDF(edfs{i});
    catch
        warning('Error loading %s...skipping',savefname);
        err{end+1} = lasterror;
        continue
    end
    save(savefname,'FIX','RAW');
    
end
cd(currdir);