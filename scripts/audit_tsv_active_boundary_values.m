function audit_tsv_active_boundary_values(mphFile, outputFile)
% Read-only dump of active HT boundary/value expressions exactly as stored.
import com.comsol.model.util.*
m=mphload(mphFile,['ddm_bc_audit_' regexprep(char(java.util.UUID.randomUUID),'-','')]);
cleanup=onCleanup(@()ModelUtil.remove(m.tag)); h=m.component('comp1').physics('ht');
fid=fopen(outputFile,'w'); if fid<0,error('Cannot write %s',outputFile);end
closeFile=onCleanup(@()fclose(fid)); tags=cellstr(char(h.feature.tags));
for q=1:numel(tags)
 f=h.feature(tags{q}); if ~f.isActive(),continue,end
 typ=char(f.getType); sel=double(f.selection.entities)';
 if ~ismember(typ,{'TemperatureBoundary','HeatFluxBoundary','HeatSource'}),continue,end
 fprintf(fid,'TAG=%s TYPE=%s SELECTION_COUNT=%d\n',tags{q},typ,numel(sel));
 for key={'T0','HeatFluxType','h','Text','q0','Q0','P0'}
  try, fprintf(fid,'%s=%s\n',key{1},char(f.getString(key{1}))); catch, end
 end
 fprintf(fid,'SELECTION=%s\n',mat2str(sel));
end
clear closeFile cleanup
end
