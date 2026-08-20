function extract_comsol_tmax_history(mphFile, outputCsv)
% Read-only extraction of Tmax at each saved sol2 time; no Study run/save.
import com.comsol.model.util.*
m=mphload(mphFile,['ddm_tmax_' regexprep(char(java.util.UUID.randomUUID),'-','')]);
cleanup=onCleanup(@()ModelUtil.remove(m.tag)); s=mphsolinfo(m,'soltag','sol2'); t=double(s.solvals(:));
fid=fopen(outputCsv,'w'); if fid<0,error('Cannot write %s',outputCsv);end
fprintf(fid,'time_s,tmax_K\n');
for k=1:numel(t)
 v=mpheval(m,'T','edim','domain','dataset','dset2','solnum',k);
 fprintf(fid,'%.17g,%.17g\n',t(k),max(v.d1));
end
fclose(fid); clear cleanup
end
