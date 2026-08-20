function comsol61_matlab_smoke(outputDir, mliDir, port)
% COMSOL 6.1 LiveLink for MATLAB smoke test. The MPHTXT is written only by
% COMSOL's mesh.export API; MATLAB never parses, converts, or rewrites it.
if nargin < 1 || isempty(outputDir), outputDir = 'D:\AI agent\codex\DDM_Schur\data\generated\two_by_three_comsol61\smoke'; end
if nargin < 2 || isempty(mliDir), mliDir = 'D:\Program Files\COMSOL\COMSOL61\Multiphysics\mli'; end
if nargin < 3 || isempty(port), port = 20361; end
if ~isfolder(outputDir), mkdir(outputDir); end
addpath(mliDir);
mphtxtPath = fullfile(outputDir, 'matlab_smoke.mphtxt');
mphPath = fullfile(outputDir, 'matlab_smoke.mph');
if isfile(mphtxtPath), delete(mphtxtPath); end
if isfile(mphPath), delete(mphPath); end
fprintf('MATLAB_RELEASE=%s\n', version('-release'));
fprintf('LIVELINK_DIR=%s\n', mliDir);
fprintf('COMSOL_SERVER_PORT=%d\n', port);
mphstart('localhost', port);
fprintf('COMSOL_VERSION=%s\n', mphversion());
import com.comsol.model.*
import com.comsol.model.util.*
model = ModelUtil.create('matlab_smoke');
model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);
model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('size', {'1','1','1'});
model.component('comp1').geom('geom1').run;
model.component('comp1').mesh.create('mesh1');
model.component('comp1').mesh('mesh1').create('ftet1', 'FreeTet');
model.component('comp1').mesh('mesh1').run;
model.component('comp1').mesh('mesh1').export(mphtxtPath);
mphsave(model, mphPath);
if ~isfile(mphtxtPath) || dir(mphtxtPath).bytes <= 0, error('COMSOL did not create a nonempty MPHTXT file.'); end
if ~isfile(mphPath) || dir(mphPath).bytes <= 0, error('COMSOL did not create a nonempty MPH file.'); end
stats = mphmeshstats(model, 'mesh1');
fprintf('MPHTXT_PATH=%s\n', mphtxtPath);
fprintf('MPHTXT_BYTES=%d\n', dir(mphtxtPath).bytes);
fprintf('MPH_PATH=%s\n', mphPath);
fprintf('MPH_BYTES=%d\n', dir(mphPath).bytes);
disp(stats);
fprintf('COMSOL61_MATLAB_SMOKE_OK\n');
ModelUtil.remove('matlab_smoke');
end
