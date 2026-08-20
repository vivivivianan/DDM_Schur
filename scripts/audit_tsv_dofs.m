function audit_tsv_dofs(mphFile,outputFile)
% Read-only DOF audit for the current COMSOL solution and raw C++ mesh.
import com.comsol.model.util.*
m=mphload(mphFile,['ddm_dof_audit_' regexprep(char(java.util.UUID.randomUUID),'-','')]);
cleanup=onCleanup(@()ModelUtil.remove(m.tag));
s=mphsolinfo(m,'soltag','sol2');
fid=fopen(outputFile,'w'); if fid<0,error('Cannot write %s',outputFile);end
closeFile=onCleanup(@()fclose(fid));
fprintf(fid,'COMSOL_SOLUTION=sol2\n');
dumpField(fid,'solvals',s,'solvals'); dumpField(fid,'sizes',s,'sizes'); dumpField(fid,'nummesh',s,'nummesh'); dumpField(fid,'sizesolvals',s,'sizesolvals'); dumpField(fid,'NU',s,'NU');
try
 x=mphxmeshinfo(m,'soltag','sol2');
 dumpField(fid,'XMesh_ndofs',x,'ndofs');
catch err
 fprintf(fid,'XMesh_error=%s\n',err.message);
end
fprintf(fid,'RAW_MESH_DOF_NOTE=C++ uses continuous P2 tetrahedral nodes; see cpp log Nodes(P2).\n');
clear closeFile cleanup
end
function dumpField(fid,label,s,key)
if isfield(s,key), v=s.(key); fprintf(fid,'%s=%s\n',label,mat2str(v)); else, fprintf(fid,'%s=<absent>\n',label); end
end
