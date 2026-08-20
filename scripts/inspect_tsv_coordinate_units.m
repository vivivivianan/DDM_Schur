function inspect_tsv_coordinate_units(mphFile, meshFile, outputFile)
% Read-only audit of COMSOL geometry units and raw MPHTXT coordinate extent.
% Does not save or alter the MPH or MPHTXT.
import com.comsol.model.util.*
m=mphload(mphFile,['ddm_unit_audit_' regexprep(char(java.util.UUID.randomUUID),'-','')]);
cleanup=onCleanup(@()ModelUtil.remove(m.tag));
g=m.component('comp1').geom('geom1');
raw=parse_mphtxt_entities(meshFile);
fid=fopen(outputFile,'w');
if fid<0, error('Cannot write %s',outputFile); end
closeFile=onCleanup(@()fclose(fid));
fprintf(fid,'SOURCE_MPH=%s\n',mphFile);
fprintf(fid,'GEOMETRY_LENGTH_UNIT=%s\n',char(g.lengthUnit()));
fprintf(fid,'GEOMETRY_NDOMAINS=%d\n',double(g.getNDomains));
fprintf(fid,'RAW_MPHTXT=%s\n',meshFile);
fprintf(fid,'RAW_COORD_MIN=%.17g,%.17g,%.17g\n',min(raw.V,[],1));
fprintf(fid,'RAW_COORD_MAX=%.17g,%.17g,%.17g\n',max(raw.V,[],1));
fprintf(fid,'RAW_COORD_EXTENT=%.17g,%.17g,%.17g\n',max(raw.V,[],1)-min(raw.V,[],1));
fprintf(fid,'RAW_TET_VOLUME_SUM=%.17g\n',totalTetVolume(raw));
clear closeFile cleanup
end

function v=totalTetVolume(m)
p=m.V; F=m.tet;
v=sum(abs(dot(p(F(:,2),:)-p(F(:,1),:),cross(p(F(:,3),:)-p(F(:,1),:),p(F(:,4),:)-p(F(:,1),:),2),2))/6);
end
