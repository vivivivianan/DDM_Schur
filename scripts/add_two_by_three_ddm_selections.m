function add_two_by_three_ddm_selections(mphFile)
% Add the regression model's explicit, model-owned DDM partition selections.
% This is intentionally the only MPH edited in this selection-based regression.
if ~isfile(mphFile), error('MPH file does not exist: %s',mphFile); end
addpath('D:\Program Files\COMSOL\COMSOL61\Multiphysics\mli'); import com.comsol.model.util.*
model=mphload(mphFile,'two_by_three_add_ddm_selections'); cleanup=onCleanup(@()ModelUtil.remove(model.tag));
c=model.component('comp1'); assignments={[1 2],[3 5],[4 6]};
for q=1:3
 tag=sprintf('ddm_sd_%03d',q); if any(strcmp(cellstr(char(c.selection.tags)),tag)), c.selection.remove(tag); end
 c.selection.create(tag,'Explicit'); c.selection(tag).label(tag); c.selection(tag).geom('geom1',3); c.selection(tag).set(assignments{q});
 S=mphgetselection(c.selection(tag));
 if double(S.dimension)~=3||~isequal(sort(double(S.entities(:))'),assignments{q}), error('Could not create expected Domain selection %s.',tag); end
 fprintf('%s: Explicit Domain IDs %s\n',tag,mat2str(double(S.entities(:))'));
end
mphsave(model,mphFile); clear cleanup
end
