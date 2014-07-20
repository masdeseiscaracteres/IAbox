% This function generates standalone packages for all the tools that
% can be found in this toolbox
% It uses the deppkg function by Jit Sarkar

% To-Do:
% - Copy duplicate files just once (duplicates are already removed from the
% printed list but still copied)
% - Show an error when a file does not exist (look for it with all different
% extensions)
% - Modify to copy the help files associated to MEX files
% - Extend to deal with different architecture MEX files (e.g. using
% exist('spqr_solve_tbb.mexmaci64'))
% - Be case-sensitive
function makePackage(package_dir,varargin)

package_dir=fullfile(package_dir,' '); % Convert path separator to filesep and ensure we add a trailing filesep
package_dir=package_dir(1:end-2); % Remove trailing filesep and space
[~, package_name]=fileparts(package_dir);

package_files=varargin(1:end);

fprintf('Building package %s in %s\n',package_name,package_dir);

%Remove redundant files
file_list=[];
for iMainFile=1:length(package_files)
    file_list=[file_list; deppkg(package_files{iMainFile},package_dir,true)];
end
file_list=unique(file_list);

fprintf('The following files have been INCLUDED in the package:\n');
for iFile=1:length(file_list)
    fprintf('\t%s\n',file_list{iFile});
end
fprintf('\n');
end

%%
% Function to find and consolidate all file dependencies for an m-file
% script/function or all such files in a given directory
% Anything that falls under the matlab root path will be ommitted
% (e.g. toolboxes, and standard functions).
%
% Jit Sarkar
% MPL/SIO
% 2011/04/20
% Modified by Óscar González
% 2014/03/02
%
%
%	Copy_List	=	deppkg(SRC,DEST)
%
%
% Input Parameters:
%	SRC		{str}(n)			:	Source m-file script/function to parse
%	DEST	{srt}(m)			:	Destination directory to copy
%									dependencies to
%
%
% Output Parameters:
%	Copy_List	{cell}{str}(Nc)	:
%

function	Copy_List	=	deppkg(SRC, DEST, IS_recursive)
%	IS_recursive is only used internally to trigger actual file copying,
%	otherwise the recursion would copy your given SRC file(s) as well.

%%	Input parsing
if	~exist(SRC, 'file')
    error('Input file/folder does not exist');
end

if	~exist('IS_recursive','var')
    IS_recursive	=	false;
end


%%	If SRC is a directory find top level files
%	will NOT process nested directories, because that can get messy if the
%	DEST folder is a subfolder in the SRC, and will take unecessarily long
if	isdir(SRC)
    File_List	=	dir(SRC);
    Copy_List	=	[];
    for	nn	=	1:length(File_List)
        if	~File_List(nn).isdir
            file_name	=	fullfile(SRC,File_List(nn).name);
            Copy_List	=	[Copy_List; deppkg(file_name,DEST, IS_recursive)]; %#ok<AGROW>
        end
    end
    return;
end


%%	Otherwise find all top-level dependencies for current file

try
    s=getcallinfo(which(SRC));
    list1=s(1).calls.fcnCalls.names(:); %Missing packages
    [Dep_List_Local, idxa, idxb]=unique(list1);
    Lines_Local=s(1).calls.fcnCalls.lines(idxa);
    % Dep_List={which(SRC)};
    Notfound_List={};
    for iFile=1:length(Dep_List_Local)
        switch exist(Dep_List_Local{iFile})
            case {2,3,6} %.m,.mex or .p
                %             full_filename = which(Dep_List_Local{iFile});
                %             IN_matlab	=	strncmp(matlabroot, full_filename, length(matlabroot));
                %             if not(IN_matlab)
                %             Dep_List=[Dep_List; full_filename];
                %             end
            case 5 %built-in
            case 0 %not-found
                Notfound_List=[Notfound_List; Dep_List_Local{iFile}];
                %fprintf('%s NOT FOUND at %s, line %d\n',Dep_List_Local{iFile},which(SRC),Lines_Local(iFile))
                fprintf('%s NOT FOUND at <a href="matlab: opentoline(which(''%s''),%d)">%s, line %d</a>\n',Dep_List_Local{iFile},which(SRC),Lines_Local(iFile),SRC,Lines_Local(iFile))
                
            otherwise
                warning('This kind of file is not considered');
        end
    end
end

Dep_List	=	depfun(SRC, '-quiet', '-toponly'); %Missing not found files
% %	Determine which dependencies are matlab bundled functions/toolboxes
IN_matlab	=	strncmp(matlabroot, Dep_List, length(matlabroot));
% %	Remove them from the list
Dep_List	=	Dep_List(~IN_matlab);

%	First item is always the current file itself
%	Only copy the file, if this is a recursive step
Copy_List	=	[];
if	IS_recursive
    if	~isdir(DEST);	mkdir(DEST);	end;
    copyfile(Dep_List{1},DEST);
    Copy_List	=	Dep_List(1);
end

% Find MEX files
for	nn	=	2:length(Dep_List)
    %   Find MEX files
    if exist(Dep_List{nn},'file')==3
        [pathstr,name,ext] = fileparts(Dep_List{nn});
        exts=mexext('all');
        %         exts(end+1).ext='m'; %TODO: copy also associated .m file
        for iExt=1:length(exts)
            arch_filename=fullfile(pathstr,[name '.' exts(iExt).ext]);
            if exist(arch_filename,'file')
                Dep_List=[Dep_List; arch_filename];
            end
        end
    end
end

%	Loop through reduced dependency list recursively finding all other
%	required files
for		nn	=	2:length(Dep_List)
    %Dep_List{nn}
    Copy_List	=	[Copy_List; deppkg(Dep_List{nn},DEST, true)]; %#ok<AGROW>
end



%%	End of Function
return;
end