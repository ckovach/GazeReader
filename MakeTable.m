
function MakeTable(varargin)

%Generates a latex table

X = varargin{1};
Head = {};
FileName = 'Table.tex';
fontsize = 'normalsize';
appnum = 1;
while exist(FileName)
    appnum = appnum+1;
    FileName = sprintf('Table%i.tex',appnum);
end


varargin{end+1} = 'finis';
i = 2;
while i <= length(varargin)
    switch lower(varargin{i})
       case 'filename'
            FileName = varargin{i+1};
            i = i+1;
       case 'head'
            Head = varargin{i+1};
            i = i+1;
       case 'fontsize'
            fontsize = varargin{i+1};
            i = i+1;
       case 'finis'
       otherwise
            error([varargin{i},' is not a valid keyword.'])
    end
    i = i+1;
end
            
FileName = [strtok(FileName,'.'),'.tex'];
fid = fopen(FileName,'w');
fprintf(fid,'\\documentclass{article}\n');
fprintf(fid,'\\addtolength{\\hoffset}{-1in}\n');
fprintf(fid,'\\addtolength{\\marginparwidth}{-1in}\n');
fprintf(fid,'\\addtolength{\\voffset}{-1in}\n');
fprintf(fid,'\n\\begin{document}');
fprintf(fid,'\n\\pagestyle{empty}');
fprintf(fid,'\n\\begin{%s}',fontsize);
fprintf(fid,'\n\t\\begin{tabular}{%s}',108*ones(1,size(X,2))); %Ignore warning
fprintf(fid,'\n\t');
for i = 1:length(Head)
    if i == 1
        fprintf(fid,'%s',Head{i});
    else
        fprintf(fid,'& %s',Head{i});
    end
end
if ~isempty(Head)
    fprintf(fid,'\\\\\n\\hline')  ;
end

for i = 1:size(X,1)
  
    for j = 1:size(X,2)        
        x = X(i,j);
        if iscell(x)
            x = x{1};
        end
        if isnan(x)
            x = [];
        end
        if isstr(x)
            prstr = '%s';
        else
            prstr = '%0.4g';
        end
        if j == 1
           fprintf(fid,['\n\t',prstr],x);
        else
           	fprintf(fid,[' & ',prstr],x);
        end
        
    end
        
    fprintf(fid,'\\\\');
end
   
fprintf(fid,'\n\t\\end{tabular}')  ;
fprintf(fid,'\n\\end{%s}',fontsize);
fprintf(fid,'\n\\end{document}');

fclose(fid)    