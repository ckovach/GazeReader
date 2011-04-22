function gr_present = GRcheck

% This function checks if GazeReader (GR) is present in the current path. If
% not, it offers to attempt a fresh installation through svn. It returns
% a value of 1 if GR is present or successfully installed and 0 otherwise.

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------



if ~exist('GazeReader.m','file')
    beep
    
    grurl = 'https://saccade.neurosurgery.uiowa.edu/GazeReader/0.1';
   fprintf(['\n--------MISSING FILES---------\nNecessary files appear not to be in the current path.\n\nPlease install GazeReader and add to path\n',...
            '%s\n\n',grurl]) 
    inp = 'x';    
    while ~ismember(lower(inp),{'y','n'})
        inp = input(sprintf(['\nDo you want me to try to install GazeReader now?\n(This requires',...
                    ' command line svn, and username/passwd access)\nY/N:']),'s');
    end
    
    if isequal(lower(inp),'y')
        com = sprintf('svn co %s GazeReader',grurl);
        [stat,res] = system(com);
        if stat >0
            beep
            if ispc                  
                    url = sprintf('download a free version at\n\n\t\thttp://www.sliksvn.com/en/download');
            elseif ismac 
                    url = sprintf('download a free version at\n\n\t\thttp://www.open.collab.net/downloads/community/');
            elseif isunix
                    url = sprintf('install it from your \nrepository (eg run  ''sudo apt-get install svn'')');
            else
                    url = sprintf('find\na version suitable for your system at\n\n\thttp://subversion.apache.org/packages.html');                    
            end
            fprintf(['\n Installation failed with error \n\n\t%s.\n\nIf this happened because you don''t have a working',...
                     ' copy of SVN, you can %s\n\n'],res,url)
        end
        
    else
        stat = 1;
    end
    
    if stat == 0
        addpath(fullfile(cd,'GazeReader'))
    end
    
else
    pth = fileparts(which('GazeReader'));
    fprintf('\nGazeReader installed at %s\n',pth)
    stat = 0;
end

gr_present = stat==0;