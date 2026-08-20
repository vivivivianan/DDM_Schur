function run_tsv_comsol_reference_summary(mphFile,outputFile)
% Run the existing COMSOL transient study in memory only; never save the MPH.
import com.comsol.model.util.*
m=mphload(mphFile,['tsv_reference_' regexprep(char(java.util.UUID.randomUUID),'-','')]); c=onCleanup(@()ModelUtil.remove(m.tag));
m.study('std2').run;
values=mpheval(m,'T','edim','domain','solnum','end'); T=double(values.d1(:));
f=fopen(outputFile,'w'); fprintf(f,'SOURCE=%s\nSTUDY=std2/time\nTIME=200 s\nSAMPLES=%d\nTMIN_K=%.17g\nTMAX_K=%.17g\nTAVG_K=%.17g\n',mphFile,numel(T),min(T),max(T),mean(T)); fclose(f); clear c
end
