function dump_comsol_model_raw(mphFile,outputFile)
% Read-only diagnostic dump used when supported-feature inspection stops early.
import com.comsol.model.util.*
m=mphload(mphFile,['rawdump_' regexprep(char(java.util.UUID.randomUUID),'-','')]); c=onCleanup(@()ModelUtil.remove(m.tag)); f=fopen(outputFile,'w');
for ct=cellstr(char(m.component.tags))', comp=ct{1}; cp=m.component(comp); fprintf(f,'COMPONENT %s\nGEOMETRIES %s\nPHYSICS %s\nMATERIALS %s\n',comp,strjoin(cellstr(char(cp.geom.tags)),','),strjoin(cellstr(char(cp.physics.tags)),','),strjoin(cellstr(char(cp.material.tags)),',')); for mt=cellstr(char(cp.material.tags))', x=cp.material(mt{1});fprintf(f,'MATERIAL tag=%s selection=%s k=%s rho=%s cp=%s\n',mt{1},mat2str(double(x.selection.entities)'),safe(x,'thermalconductivity'),safe(x,'density'),safe(x,'heatcapacity'));end; for pt=cellstr(char(cp.physics.tags))',p=cp.physics(pt{1});for ft=cellstr(char(p.feature.tags))',x=p.feature(ft{1});fprintf(f,'FEATURE physics=%s tag=%s type=%s active=%d selection=%s Tinit=%s T0=%s HeatFluxType=%s h=%s Text=%s q0=%s Q0=%s\n',pt{1},ft{1},char(x.getType),x.isActive(),mat2str(double(x.selection.entities)'),safeF(x,'Tinit'),safeF(x,'T0'),safeF(x,'HeatFluxType'),safeF(x,'h'),safeF(x,'Text'),safeF(x,'q0'),safeF(x,'Q0'));end,end;end
fclose(f);clear c
end
function v=safe(x,key),try,v=char(x.propertyGroup('def').getString(key));catch,v='<unavailable>';end,end
function v=safeF(x,key),try,v=char(x.getString(key));catch,v='<unavailable>';end,end
