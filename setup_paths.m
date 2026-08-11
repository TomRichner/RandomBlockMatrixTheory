function setup_paths()
% SETUP_PATHS Add the RandomBlockMatrixTheory source tree to the MATLAB path.
%
% Run once per session. Idempotent, and never calls savepath.
%
% The function derives everything from its own location, so it resolves from a
% cold session as long as the working directory is the repo root:
%
%   >> setup_paths
%
% After that, anything in src/, tests/ and examples/ runs from any cwd.

this_dir = fileparts(mfilename('fullpath'));

folders = { ...
    fullfile(this_dir, 'src'), ...
    fullfile(this_dir, 'tests'), ...
    fullfile(this_dir, 'examples') ...
    };

for k = 1:numel(folders)
    if exist(folders{k}, 'dir')
        addpath(genpath(folders{k}));
    end
end

end
