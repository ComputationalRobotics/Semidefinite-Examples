function install_mex(recompile)
%% need matlab version >= R2016x
    if (nargin == 0); recompile = 0; end
    
    os = computer;
    if (contains(os, 'MAC'))
        mexcmd = 'mex -largeArrayDims -output ';
    else
        mexcmd = 'mex -O -largeArrayDims -output ';
    end
    %% go to source file
    eval(['cd ', 'mexfiles']);
    mexname{1} = 'mexeig';
    mexname{2} = 'mexFnorm';
    mexname{3} = 'mexMatvec';
    mexname{4} = 'mexnnz';
    mexname{5} = 'mexsmat';
    mexname{6} = 'mexsvec';
    mexname{7} = 'mextriang';
    mexname{8} = 'mextriangsp';
    %% compile
    for k = 1:length(mexname)
        if (recompile)
            cmd([mexcmd, mexname{k}, ' ', mexname{k}, '.c'])
        end
    end
end

function cmd(s)
    fprintf('\n %s', s);
    eval(s);
end