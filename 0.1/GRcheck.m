function gr_present = GRcheck

% This function checks if GazeReader (GR) is present in the current path. If
% not, it offers to attempt a fresh installation through svn. Installation only works
% if a command-line implementation of svn is installed on the system. It returns
% a value of 1 if GR is present or successfully installed and 0 otherwise.

% ----------- SVN REVISION INFO ------------------
% $URL: file:///home/svn/Tasks/Neuroecon/intracranial/AttentionValue/Incentive_Salience/GRcheck.m $
% $Revision: 35 $
% $Date: 2011-04-20 13:13:10 -0500 (Wed, 20 Apr 2011) $
% $Author: ckovach $
% ------------------------------------------------

%C Kovach 2011

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
            fprintf(['\n Installation failed with error \n\n\t%s.\nIf this happened because you don''t have a working',...
                     'copy of SVN, download at \n\thttp://subversion.apache.org/packages.html\n\n'],res)
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