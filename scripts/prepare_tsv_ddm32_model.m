function prepare_tsv_ddm32_model(sourceMph, targetMph, reportFile, dryRun)
% Split the TSV model into 32 x-direction domain groups in a new MPH only.
% `dryRun=true` mutates only the in-memory copy and never saves it.
if nargin<4, dryRun=true; end
import com.comsol.model.util.*
m=mphload(sourceMph,['tsv_ddm32_' regexprep(char(java.util.UUID.randomUUID),'-','')]);
cleanup=onCleanup(@()ModelUtil.remove(m.tag)); c=m.component('comp1'); g=c.geom('geom1'); h=c.physics('ht');
if ~strcmp(char(g.lengthUnit()),'m'), error('Expected geometry length unit m.'); end
if ~dryRun && isfile(targetMph), error('Refusing to overwrite existing target MPH: %s',targetMph); end
for cut=1:31
  x=(3+6*(cut-1))*1e-6; wp=sprintf('ddm32_wp_%02d',cut); par=sprintf('ddm32_par_%02d',cut);
  g.create(wp,'WorkPlane'); g.feature(wp).set('quickplane','yz'); g.feature(wp).set('quickx',sprintf('%.17g',x));
  g.create(par,'PartitionDomains'); g.feature(par).selection('domain').all; g.feature(par).set('partitionwith','workplane'); g.feature(par).set('workplane',wp);
end
g.run;
solid=h.feature('solid1'); heatDomains=sort(double(solid.selection.entities)'); n=double(g.getNDomains);
if numel(heatDomains)~=n || ~isequal(heatDomains,1:n), error('Post-partition Heat Transfer selection does not cover every Domain exactly once. geometry=%d heat=%d',n,numel(heatDomains)); end
% Every post-partition Domain must lie wholly inside one 6 um x cell.
groups=cell(32,1); crossings=[];
for d=heatDomains
 p=mphgetcoords(m,'geom1','domain',d); lo=min(p(1,:)); hi=max(p(1,:)); xc=mean(p(1,:)); idx=floor((xc+3e-6)/6e-6)+1; idx=max(1,min(32,idx));
 if idx==1, left=-4e-6; right=3e-6; else, left=3e-6+(idx-2)*6e-6; right=left+6e-6; end
 if idx==32,right=190e-6;end
 if lo<left-1e-12 || hi>right+1e-12, crossings(end+1)=d; continue,end %#ok<AGROW>
 groups{idx}(end+1)=d; %#ok<AGROW>
end
if ~isempty(crossings), error('Domains still cross a 32-way x cut: %s',mat2str(crossings)); end
if any(cellfun(@isempty,groups)), error('At least one ddm_sd group is empty.'); end
% A source model may already carry a previous DDM partition.  Remove only
% selections identified by the public ddm_sd_### naming convention (tag or
% label) before installing this derived 32-way partition.
tags=cellstr(char(c.selection.tags));
for q=1:numel(tags)
 tagMatch=~isempty(regexp(tags{q},'^ddm_sd_\d+$','once'));
 try, oldLabel=char(c.selection(tags{q}).label()); catch, oldLabel=''; end
 labelMatch=~isempty(regexp(oldLabel,'^ddm_sd_\d+$','once'));
 if tagMatch||labelMatch, c.selection.remove(tags{q}); end
end
for q=1:32
 tag=sprintf('ddm_sd_%03d',q); c.selection.create(tag,'Explicit'); c.selection(tag).label(sprintf('DDM x group %03d',q)); c.selection(tag).geom('geom1',3); c.selection(tag).set(groups{q});
end
f=fopen(reportFile,'w'); if f<0,error('Cannot write %s',reportFile),end
fprintf(f,'SOURCE=%s\nTARGET=%s\nDRY_RUN=%d\nPOST_PARTITION_GEOMETRY_DOMAINS=%d\nPOST_PARTITION_HEAT_DOMAINS=%d\n',sourceMph,targetMph,dryRun,n,numel(heatDomains));
for q=1:32, fprintf(f,'ddm_sd_%03d,count=%d,domains=%s\n',q,numel(groups{q}),mat2str(groups{q}));end; fclose(f);
if ~dryRun, mphsave(m,targetMph); end
clear cleanup
end
