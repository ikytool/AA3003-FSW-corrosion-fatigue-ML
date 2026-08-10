%RUN_ALL Regenerate every figure that the current data allow.
%   Stages s1-s3 always run (published constants). Stages s4-s6 need the
%   digitized dataset and are skipped with a message until it exists.
%   (Stage names are written out literally because the stage scripts
%   clear the workspace.)

cd(fileparts(mfilename('fullpath')));

s1_paris_baseline;
s2_arrhenius_crosscheck;
s3_inhibitor_efficiency;

try, s4_ml_hybrid;         catch err, fprintf('--- s4_ml_hybrid skipped: %s\n', err.message); end
try, s5_life_integration;  catch err, fprintf('--- s5_life_integration skipped: %s\n', err.message); end
try, s6_probabilistic;     catch err, fprintf('--- s6_probabilistic skipped: %s\n', err.message); end

fprintf('Done. Figures are in ../figures.\n');
