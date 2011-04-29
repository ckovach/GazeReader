function gr_present = GRcheck

% This function checks if GazeReader (GR) is present in the current path. If
% not, it offers to attempt a fresh installation through svn. It returns
% a value of 1 if GR is present or successfully installed and 0 otherwise.
%
% Subversion can be downloaded at http://subversion.apache.org/packages.html
%

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------


%%% URL for current GazeReader repository
grurl = 'https://saccade.neurosurgery.uiowa.edu/GazeReader/0.1';

if ~exist('GazeReader.m','file')
    beep
    
    
   fprintf(['\n--------MISSING FILES---------\nNecessary files appear not to be in the current path.\n\nPlease install GazeReader and add to path\n',...
            '%s\n\n',grurl]) 
    inp = 'x';    
    while ~ismember(lower(inp),{'y','n'})
        inp = input(sprintf(['\nDo you want me to try to install GazeReader now?\n(This requires',...
                    ' command line svn, and username/passwd access)\nY/N:']),'s');
    end
    
    if isequal(lower(inp),'y')
        
        
        fprintf('\nChoose a place to install..')
        installdir = '';
        while exist(fullfile(installdir,'.svn'),'dir') > 0
           installdir = fullfile('..',installdir); 
        end
                
        installdir = uigetdir(installdir,'Choose where to install GazeReader.');
        while exist(fullfile(installdir,'.svn'),'dir') > 0  && ~exist(fullfile(installdir,'GazeReader.m'),'file')
            warndlg('Selected directory must not be a working copy of a subversion repository. Please choose another.')
            installdir = uigetdir(installdir,'Selected directory must not be a subversion repository.');        
            if installdir == 0
                return
            end
        end    
        
        if exist(installdir,'dir') > 0 && exist(fullfile(installdir,'GazeReader.m'),'file')
            fprintf('\nAdding existing version to the matlab path.')
            addpath(installdir)
            gr_present = true;
            return
        elseif exist(installdir,'dir') > 0
            grpath = fullfile(installdir,'GazeReader');
        else
            grpath = installdir;
        end
        
        com = sprintf('svn co %s %s',grurl,grpath);
        
        %%% Attempt to checkout using subversion at the system command-line
        [stat,res] = system(com);  
        
        if stat >0  %%% If command returned an error
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

        else
            addpath(grpath)
        end
        
    else
        stat = 1;
    end
    
    
else
    pth = fileparts(which('GazeReader'));
    fprintf('\nGazeReader installed at %s\n',pth)
    stat = 0;
end

gr_present = stat==0;