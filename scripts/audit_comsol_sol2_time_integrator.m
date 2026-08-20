function audit_comsol_sol2_time_integrator(mphFile,outputFile)
% Read-only audit of the stored sol2 solver sequence.
import com.comsol.model.util.*
m=mphload(mphFile,['ddm_integrator_audit_' regexprep(char(java.util.UUID.randomUUID),'-','')]);
cleanup=onCleanup(@()ModelUtil.remove(m.tag)); sol=m.sol('sol2');
fid=fopen(outputFile,'w'); if fid<0,error('Cannot write %s',outputFile);end
fprintf(fid,'SOL=sol2\n'); dumpNode(fid,sol,'sol2');
for tag=cellstr(char(sol.feature.tags))', dumpNode(fid,sol.feature(tag{1}),tag{1}); end
fclose(fid); clear cleanup
end
function dumpNode(fid,node,tag)
fprintf(fid,'NODE=%s TYPE=%s\n',tag,char(node.getType));
keys={'tlist','tstepsbdf','timestepper','bdforder','maxorder','initialstepbdf','maxstepbdf','rtol','atolglobal','atol','estrat','timesel','linsolver','nonlin','maxiter'};
for k=1:numel(keys), try, fprintf(fid,'%s=%s\n',keys{k},char(node.getString(keys{k}))); catch, end, end
end
