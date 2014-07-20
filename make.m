function make(command)
switch lower(command)
    case 'install'
        fprintf('IAbox has been added to search path (until the end of this session)\n');
        addpath(genpath('example')); % Usage examples
        addpath(genpath('src')); % Core source code
        addpath(genpath('test')); % Scripts testing proper functioning of files in "src"
    case 'uninstall'
        rmpath(genpath('example')); % Usage examples
        rmpath(genpath('src')); % Core source code
        rmpath(genpath('test')); % Scripts testing proper functioning of files in "src"
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