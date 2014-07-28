function make(command)
switch lower(command)
    case 'install'
        addpath(genpath(fullfile(pwd,'example'))); % Usage examples
        addpath(genpath(fullfile(pwd,'src'))); % Core source code
        addpath(genpath(fullfile(pwd,'test'))); % Scripts testing proper functioning of files in "src"
        fprintf('IAbox has been added to search path (until the end of this session)\n');
    case 'uninstall'
        rmpath(genpath(fullfile(pwd,'example'))); % Usage examples
        rmpath(genpath(fullfile(pwd,'src'))); % Core source code
        rmpath(genpath(fullfile(pwd,'test'))); % Scripts testing proper functioning of files in "src"
        fprintf('IAbox has been removed from search path\n');
    case 'packages'
        make install;
        % Define your own packages here
        makePackage StandalonePackages/ClosedForm3Users ClosedForm3Users;
        makePackage StandalonePackages/GaussNewtonMinLeakage test_GaussNewtonMinLeakage example_NumberSolutionsViaGaussNewton;
        %
    otherwise
        error('Unrecognized option');
end