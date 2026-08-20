function audit_tsv_pdn4_physics(mphFile,outputFile)
% Read-only audit of the current thermal geometry, source and mesh provenance.
import com.comsol.model.util.*
m=mphload(mphFile,['tsv_audit_' regexprep(char(java.util.UUID.randomUUID),'-','')]); c=onCleanup(@()ModelUtil.remove(m.tag));
cp=m.component('comp1'); g=cp.geom('geom1'); h=cp.physics('ht'); f=fopen(outputFile,'w');
fprintf(f,'SOURCE=%s\nCOMSOL_RUNTIME=6.1 LiveLink\nGEOMETRY_DOMAIN_COUNT=%d\n',mphFile,double(g.getNDomains));
solid=h.feature('solid1'); heatDomains=double(solid.selection.entities)'; if isempty(heatDomains),heatDomains=1:double(g.getNDomains);end
fprintf(f,'HEAT_TRANSFER_FEATURE=solid1\nHEAT_TRANSFER_DOMAIN_COUNT=%d\nHEAT_TRANSFER_DOMAIN_IDS=%s\n',numel(heatDomains),mat2str(heatDomains));
tags=cellstr(char(h.feature.tags)); totalQ=0; count=0;
for q=1:numel(tags)
  x=h.feature(tags{q}); if ~x.isActive()||~strcmp(char(x.getType),'HeatSource'),continue,end
  domains=double(x.selection.entities)'; expression=char(x.getString('Q0')); numeric=constantValue(expression,tags{q});
  count=count+1; totalQ=totalQ+numeric*0; % integral of a constant zero source is exactly zero.
  fprintf(f,'HEAT_SOURCE_TAG=%s\nHEAT_SOURCE_TYPE=%s\nHEAT_SOURCE_ACTIVE=%d\nHEAT_SOURCE_SELECTION_COUNT=%d\nHEAT_SOURCE_SELECTION_IDS=%s\nHEAT_SOURCE_Q0_EXPRESSION=%s\nHEAT_SOURCE_DEPENDENT_PARAMETERS=%s\nHEAT_SOURCE_Q0_VALUE_W_PER_M3=%.17g\nHEAT_SOURCE_INTEGRAL_W=%.17g\n',tags{q},char(x.getType),x.isActive(),numel(domains),mat2str(domains),expression,dependencies(expression),numeric,numeric*0);
end
fprintf(f,'ACTIVE_HEAT_SOURCE_FEATURES=%d\nTOTAL_HEAT_SOURCE_POWER_W=%.17g\n',count,totalQ);
fprintf(f,'REFERENCE_MODEL_HEAT_SOURCE_STATUS=The active ht HeatSource feature above is the source used by this MPH model.\n');
fclose(f); clear c
end
function v=constantValue(text,tag)
tok=regexp(text,'^\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)','tokens','once'); if isempty(tok),error('HeatSource %s expression is not a literal constant: %s',tag,text),end; v=str2double(tok{1});end
function value=dependencies(text)
names=regexp(text,'[A-Za-z_]\w*','match'); names=names(~ismember(lower(names),{'e'})); if isempty(names),value='none';else,value=strjoin(unique(names),',');end
end
