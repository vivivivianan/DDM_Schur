function inspect_comsol_physics_inventory(mphFile,outputFile)
% Read-only inventory of active physics features without assuming properties.
import com.comsol.model.util.*
m=mphload(mphFile,['physics_inventory_' regexprep(char(java.util.UUID.randomUUID),'-','')]); c=onCleanup(@()ModelUtil.remove(m.tag));
f=fopen(outputFile,'w'); cp=m.component('comp1');
for p=cellstr(char(cp.physics.tags))'
  ph=cp.physics(p{1}); fprintf(f,'PHYSICS tag=%s type=%s\n',p{1},char(ph.getType()));
  for q=cellstr(char(ph.feature.tags))'
    x=ph.feature(q{1}); fprintf(f,'FEATURE tag=%s type=%s active=%d selection_count=%d\n',q{1},char(x.getType()),x.isActive(),numel(double(x.selection.entities)));
  end
end
fclose(f); clear c
end
