function inspect_comsol_case(mphFile,outputDir)
% Read-only inspection front-end for an existing COMSOL model.
if ~isfile(mphFile), error('MPH not found: %s',mphFile); end
if ~isfolder(outputDir), mkdir(outputDir); end
import com.comsol.model.util.*
model=mphload(mphFile,['inspect_' regexprep(char(java.util.UUID.randomUUID),'-','')]); cleanup=onCleanup(@()ModelUtil.remove(model.tag));
cfg=struct('component_tag','','geometry_tag','','heat_physics_tag','');
try
 info=inspect_comsol_heat_model(model,cfg);
catch ME
 fid=fopen(fullfile(outputDir,'model_inspection.txt'),'w'); fprintf(fid,'SOURCE MPH: %s\n\nINSPECTION FAILED\n%s\n',mphFile,getReport(ME,'extended','hyperlinks','off')); fclose(fid); rethrow(ME)
end
fid=fopen(fullfile(outputDir,'model_inspection.txt'),'w');
fprintf(fid,'SOURCE MPH: %s\nCOMSOL runtime: %s\nMATLAB: %s\ncomponent_tag: %s\ngeometry_tag: %s\nheat_physics_tag: %s\nDomains: %s\nBoundaries: %s\n\nMATERIALS\n',mphFile,info.comsolVersion,version,info.component,info.geometry,info.heatPhysics,mat2str(info.domains),mat2str(info.boundaries));
for q=1:numel(info.materials),m=info.materials(q);fprintf(fid,'tag=%s domains=%s k=%.17g rho=%.17g Cp=%.17g\n',m.tag,mat2str(m.domains),m.k,m.rho,m.cp);end
p=info.physics;fprintf(fid,'\nINITIAL TEMPERATURE: %.17g K\nDISABLED HEAT FEATURES: %s\nDIRICHLET [boundary,T]:\n%s\nCONVECTION [boundary,h,Tinf]:\n%s\nHEAT FLUX [boundary,q_in]:\n%s\nHEAT SOURCE DENSITY [domain,Q]:\n%s\nHEAT SOURCE TOTAL POWER:\n',p.initialTemperature,strjoin(p.disabledFeatures,', '),mat2str(p.dirichlet),mat2str(p.convection),mat2str(p.heatflux),mat2str(p.source));for q=1:numel(p.sourceTotalPower),x=p.sourceTotalPower(q);fprintf(fid,'tag=%s domains=%s P=%.17g W\n',x.tag,mat2str(x.domains),x.powerW);end;fprintf(fid,'THERMAL INSULATION BOUNDARIES:\n%s\n',mat2str(p.insulation));
st=cellstr(char(model.study.tags));fprintf(fid,'\nSTUDIES\n');for q=1:numel(st),ft=cellstr(char(model.study(st{q}).feature.tags));fprintf(fid,'%s features=%s\n',st{q},strjoin(ft,','));for k=1:numel(ft),if strcmpi(ft{k},'time'),fprintf(fid,'  tlist=%s\n',char(model.study(st{q}).feature(ft{k}).getString('tlist')));end,end,end
fprintf(fid,'\nCOMPONENT NAMED SELECTIONS\n%s\n',strjoin(cellstr(char(model.component(info.component).selection.tags)),','));fclose(fid);
fid=fopen(fullfile(outputDir,'comsol_domain_map.csv'),'w');fprintf(fid,'domain_id,material_tag,k,rho,Cp\n');for d=info.domains,m=info.materials(arrayfun(@(x)ismember(d,x.domains),info.materials));fprintf(fid,'%d,%s,%.17g,%.17g,%.17g\n',d,m.tag,m.k,m.rho,m.cp);end;fclose(fid);
fid=fopen(fullfile(outputDir,'comsol_boundary_map.csv'),'w');fprintf(fid,'boundary_id,adjacent_domain_ids\n');for b=info.boundaries,fprintf(fid,'%d,"%s"\n',b,strjoin(arrayfun(@num2str,info.boundaryDomains{b},'UniformOutput',false),';'));end;fclose(fid);
clear cleanup
end
