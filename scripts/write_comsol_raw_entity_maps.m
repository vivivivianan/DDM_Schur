function write_comsol_raw_entity_maps(mphFile,out)
% Read-only fallback: entity maps remain available even when physics is unsupported.
import com.comsol.model.util.*
m=mphload(mphFile,['maps_' regexprep(char(java.util.UUID.randomUUID),'-','')]); c=onCleanup(@()ModelUtil.remove(m.tag)); ct=cellstr(char(m.component.tags)); if numel(ct)~=1,error('Ambiguous component'),end;cp=m.component(ct{1});gt=cellstr(char(cp.geom.tags));if numel(gt)~=1,error('Ambiguous geometry'),end;g=cp.geom(gt{1});nd=double(g.getNDomains);nb=double(g.getNBoundaries);mats=cellstr(char(cp.material.tags));owner=cell(1,nd);props=containers.Map;
for q=1:numel(mats),x=cp.material(mats{q});props(mats{q})={safe(x,'thermalconductivity'),safe(x,'density'),safe(x,'heatcapacity')};for d=double(x.selection.entities)',if d>=1&&d<=nd,owner{d}=mats{q};end,end,end
f=fopen(fullfile(out,'comsol_domain_map.csv'),'w');fprintf(f,'domain_id,material_tag,k,rho,Cp\n');for d=1:nd,if isempty(owner{d}),owner{d}='<unassigned>';fprintf(f,'%d,%s,,,\n',d,owner{d});else,p=props(owner{d});fprintf(f,'%d,%s,%s,%s,%s\n',d,owner{d},p{1},p{2},p{3});end,end;fclose(f);
f=fopen(fullfile(out,'comsol_boundary_map.csv'),'w');fprintf(f,'boundary_id,adjacent_domain_ids\n');for b=1:nb,a=double(mphgetadj(m,gt{1},'domain','boundary',b))';fprintf(f,'%d,"%s"\n',b,strjoin(arrayfun(@num2str,a,'UniformOutput',false),';'));end;fclose(f);
f=fopen(fullfile(out,'model_selection_inventory.txt'),'w');fprintf(f,'component=%s\ngeometry=%s\ncomponent selections=%s\nstudies:\n',ct{1},gt{1},strjoin(cellstr(char(cp.selection.tags)),','));for st=cellstr(char(m.study.tags))',ft=cellstr(char(m.study(st{1}).feature.tags));fprintf(f,'  %s features=%s\n',st{1},strjoin(ft,','));for q=1:numel(ft),if strcmpi(ft{q},'time'),fprintf(f,'    tlist=%s\n',char(m.study(st{1}).feature(ft{q}).getString('tlist')));end,end,end;fclose(f);clear c
end
function v=safe(x,key),try,v=char(x.propertyGroup('def').getString(key));catch,v='<unavailable>';end,end
