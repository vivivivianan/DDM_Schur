function audit_tsv_saved_sol2_reference(mphFile,outputFile)
% Read the user-saved sol2 only; do not run a Study or mutate the source MPH.
import com.comsol.model.util.*
m=mphload(mphFile,['ddm_saved_sol2_audit_' regexprep(char(java.util.UUID.randomUUID),'-','')]);
cleanup=onCleanup(@()ModelUtil.remove(m.tag));
s=mphsolinfo(m,'soltag','sol2'); x=mphxmeshinfo(m,'soltag','sol2');
v=mpheval(m,'T','edim','domain','dataset','dset2','solnum','end');
fid=fopen(outputFile,'w'); if fid<0,error('Cannot write %s',outputFile);end
fprintf(fid,'SOURCE=%s\nSOLTAG=sol2\nNDOFS=%d\nTIME_STEPS=%d\nTIME_START_S=%.17g\nTIME_END_S=%.17g\nTMIN_K=%.17g\nTMAX_K=%.17g\nTAVG_K=%.17g\n',mphFile,double(x.ndofs),numel(s.solvals),double(s.solvals(1)),double(s.solvals(end)),min(v.d1),max(v.d1),mean(v.d1)); fclose(fid);
clear cleanup
end
