function smoke_test_comsol_partition_domains(outputFile)
% COMSOL 6.1 API smoke test: partition one block with a yz work plane.
% The model is in-memory only and is never saved.
import com.comsol.model.*
import com.comsol.model.util.*
m=ModelUtil.create(['partition_smoke_' regexprep(char(java.util.UUID.randomUUID),'-','')]);
cleanup=onCleanup(@()ModelUtil.remove(m.tag));
g=m.component.create('comp1',true).geom.create('geom1',3);
g.create('blk1','Block'); g.feature('blk1').set('size',{'1','1','1'}); g.run;
before=double(g.getNDomains);
g.create('wp1','WorkPlane'); g.feature('wp1').set('quickplane','yz'); g.feature('wp1').set('quickx','0.5');
g.create('par1','PartitionDomains'); par=g.feature('par1');
par.selection('domain').all; par.set('partitionwith','workplane'); par.set('workplane','wp1');
g.run;
after=double(g.getNDomains);
f=fopen(outputFile,'w'); if f<0,error('Cannot write %s',outputFile),end
fprintf(f,'before_domains=%d\nafter_domains=%d\n',before,after); fclose(f);
if before~=1 || after~=2,error('PartitionDomains smoke test expected 1 -> 2 Domains, got %d -> %d.',before,after),end
clear cleanup
end
