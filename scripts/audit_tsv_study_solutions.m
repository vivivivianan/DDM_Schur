function audit_tsv_study_solutions(mphFile,outputFile)
% Run std2 in memory and inventory solution/dataset bindings before sampling.
import com.comsol.model.util.*
m=mphload(mphFile,['study_audit_' regexprep(char(java.util.UUID.randomUUID),'-','')]); c=onCleanup(@()ModelUtil.remove(m.tag));
m.study('std2').run; f=fopen(outputFile,'w');
fprintf(f,'SOLUTIONS\n');
for tag=cellstr(char(m.sol.tags))', fprintf(f,'sol=%s\n',tag{1}); end
fprintf(f,'DATASETS\n');
for tag=cellstr(char(m.result.dataset.tags))'
 d=m.result.dataset(tag{1}); fprintf(f,'dataset=%s type=%s solution=%s\n',tag{1},char(d.getType()),safe(d,'solution'));
end
fprintf(f,'STUDY_TIME_PROPERTIES\n'); st=m.study('std2').feature('time'); fprintf(f,'solnum=%s tlist=%s\n',safe(st,'solnum'),safe(st,'tlist'));
fclose(f); clear c
end
function x=safe(node,key),try,x=char(node.getString(key));catch,x='<unavailable>';end,end
