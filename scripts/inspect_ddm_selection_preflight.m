function inspect_ddm_selection_preflight(mphFile,outputFile)
% Read-only preflight: selection inventory and Time Dependent tlist.
import com.comsol.model.util.*
m=mphload(mphFile,['ddm_preflight_' regexprep(char(java.util.UUID.randomUUID),'-','')]); c=onCleanup(@()ModelUtil.remove(m.tag));
cp=m.component('comp1'); f=fopen(outputFile,'w');
fprintf(f,'COMPONENT=comp1\nGEOMETRY=%s\nPHYSICS=%s\n',strjoin(cellstr(char(cp.geom.tags)),','),strjoin(cellstr(char(cp.physics.tags)),','));
fprintf(f,'DDM SELECTIONS\n');
for tag=cellstr(char(cp.selection.tags))'
  s=cp.selection(tag{1}); label=readLabel(s); isDdm=~isempty(regexp(tag{1},'^ddm_sd_\d+$','once','ignorecase'))||~isempty(regexp(label,'^ddm_sd_\d+$','once','ignorecase'));
  if isDdm
    x=mphgetselection(s); fprintf(f,'tag=%s label=%s type=%s dimension=%d entities=%s count=%d\n',tag{1},label,char(s.getType()),double(x.dimension),mat2str(double(x.entities(:))'),numel(x.entities));
  end
end
fprintf(f,'STUDIES\n');
for st=cellstr(char(m.study.tags))'
  features=cellstr(char(m.study(st{1}).feature.tags));
  for ft=features'
    feature=m.study(st{1}).feature(ft{1}); if strcmpi(char(feature.getType()),'Transient')||strcmpi(ft{1},'time')
      try,tlist=char(feature.getString('tlist'));catch,tlist='<unavailable>';end
      fprintf(f,'study=%s feature=%s type=%s tlist=%s\n',st{1},ft{1},char(feature.getType()),tlist);
    end
  end
end
fclose(f); clear c
end
function label=readLabel(s),try,label=char(s.label());catch,label=char(s.label);end,end
