function inspect_comsol_x_partition_candidates(mphFile, outputCsv, outputReport)
% Read-only geometry audit for an x-direction DDM partition.  It deliberately
% uses COMSOL Domain IDs and geometry vertices, never mesh element positions.
import com.comsol.model.util.*
m=mphload(mphFile,['ddm_x_audit_' regexprep(char(java.util.UUID.randomUUID),'-','')]);
cleanup=onCleanup(@()ModelUtil.remove(m.tag));
c=m.component('comp1'); g=c.geom('geom1'); h=c.physics('ht');
solid=h.feature('solid1'); domains=sort(double(solid.selection.entities)');
if isempty(domains), error('Heat Transfer solid1 selection is empty.'); end
xmin=zeros(numel(domains),1); xmax=xmin; xcenter=xmin;
for q=1:numel(domains)
  p=mphgetcoords(m,'geom1','domain',domains(q));
  if isempty(p), error('No geometry coordinates returned for Domain %d.',domains(q)); end
  xmin(q)=min(p(1,:)); xmax(q)=max(p(1,:)); xcenter(q)=mean(p(1,:));
end
T=table(domains(:),xmin,xcenter,xmax,xmax-xmin,'VariableNames',{'domain_id','xmin_m','xcenter_m','xmax_m','width_m'});
T=sortrows(T,'xcenter_m'); writetable(T,outputCsv);
u=unique(round(xcenter,15)); gaps=diff(u); [gapSorted,ii]=sort(gaps,'descend');
f=fopen(outputReport,'w'); if f<0,error('Cannot write %s',outputReport),end
closeFile=onCleanup(@()fclose(f));
fprintf(f,'SOURCE_MPH=%s\nGEOMETRY_LENGTH_UNIT=%s\nHEAT_DOMAIN_COUNT=%d\n',mphFile,char(g.lengthUnit()),numel(domains));
fprintf(f,'X_RANGE_M=%.17g,%.17g\nUNIQUE_DOMAIN_X_CENTERS=%d\n',min(xmin),max(xmax),numel(u));
fprintf(f,'\nLargest gaps between unique Domain x-centers (candidate partition cuts):\nrank,gap_m,left_center_m,right_center_m\n');
for q=1:min(40,numel(ii)), i=ii(q); fprintf(f,'%d,%.17g,%.17g,%.17g\n',q,gapSorted(q),u(i),u(i+1)); end
fprintf(f,'\nDomain x-center distribution (m):\n');
for q=1:height(T),fprintf(f,'%d,%.17g,%.17g,%.17g\n',T.domain_id(q),T.xmin_m(q),T.xcenter_m(q),T.xmax_m(q));end
clear closeFile cleanup
end
