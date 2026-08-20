function run_tsv_comsol_reference_forced(mphFile,outputFile)
% Recompute std2 in memory after discarding stored solution data.
% The source MPH is never saved or overwritten.
import com.comsol.model.util.*
try
 m=mphload(mphFile,['tsv_reference_forced_' regexprep(char(java.util.UUID.randomUUID),'-','')]);
 cleanup=onCleanup(@()ModelUtil.remove(m.tag));
 sol=m.sol('sol2'); sol.clearSolutionData(); m.study('std2').run;
 x=mphxmeshinfo(m,'soltag','sol2'); si=mphsolinfo(m,'soltag','sol2'); timeEnd=double(si.solvals(end));
 v=mpheval(m,'T','edim','domain','dataset','dset2','solnum','end');
 q=mphint2(m,'ht.Q','volume','dataset','dset2','selection',double(m.component('comp1').physics('ht').feature('hs1').selection.entities),'solnum','end');
 fid=fopen(outputFile,'w'); if fid<0,error('Cannot write %s',outputFile);end
 fprintf(fid,'SOURCE=%s\nSTUDY=std2/time (forced recompute)\nTIME_END_S=%.17g\nCOMSOL_NDOFS=%d\nSAMPLES=%d\nTMIN_K=%.17g\nTMAX_K=%.17g\nTAVG_K=%.17g\nINTEGRAL_ht.Q_W=%.17g\n',mphFile,timeEnd,double(x.ndofs),numel(v.d1),min(v.d1),max(v.d1),mean(v.d1),q);fclose(fid);
 clear cleanup
catch err
 fid=fopen([outputFile '.error.txt'],'w'); fprintf(fid,'%s\n',getReport(err,'extended','hyperlinks','off')); fclose(fid); rethrow(err)
end
end
