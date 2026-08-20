function integrate_tsv_heat_source(mphFile,outputFile)
% Read-only COMSOL reference integration of the actual thermal source field.
import com.comsol.model.util.*
m=mphload(mphFile,['q_integral_' regexprep(char(java.util.UUID.randomUUID),'-','')]); c=onCleanup(@()ModelUtil.remove(m.tag));
m.study('std2').run; hs=m.component('comp1').physics('ht').feature('hs1'); d=double(hs.selection.entities);
q=mphint2(m,'ht.Q','volume','dataset','dset2','selection',d,'solnum','end');
f=fopen(outputFile,'w'); fprintf(f,'FEATURE=hs1\nTYPE=%s\nACTIVE=%d\nP0=%s W\nQ0=%s W/m^3\nSELECTION_DOMAINS=%d\nINTEGRAL_ht.Q_W=%.17g\n',char(hs.getType()),hs.isActive(),char(hs.getString('P0')),char(hs.getString('Q0')),numel(d),q); fclose(f); clear c
end
