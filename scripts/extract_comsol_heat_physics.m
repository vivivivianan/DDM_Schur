function p=extract_comsol_heat_physics(model,info)
h=model.component(info.component).physics(info.heatPhysics); tags=cellstr(char(h.feature.tags)); p=struct('initialTemperature',[],'dirichlet',[],'convection',[],'heatflux',[],'source',[],'sourceTotalPower',struct('tag',{},'domains',{},'powerW',{}),'insulation',[],'disabledFeatures',{{}});
for q=1:numel(tags)
 f=h.feature(tags{q});
 % COMSOL retains disabled nodes in the model tree. They must be recorded but
 % must not participate in the exported physical problem.
 if ~f.isActive(), p.disabledFeatures{end+1}=sprintf('%s (%s)',tags{q},char(f.getType)); continue, end
 typ=char(f.getType); sel=double(f.selection.entities)'; if isempty(sel),continue,end
 switch typ
 case 'TemperatureBoundary', p.dirichlet=[p.dirichlet; rows(sel,const(f.getString('T0'),tags{q}),tags{q})];
 case 'HeatFluxBoundary'
  mode=char(f.getString('HeatFluxType')); if contains(lower(mode),'convective'), p.convection=[p.convection; rows2(sel,const(f.getString('h'),tags{q}),const(f.getString('Text'),tags{q}),tags{q})]; else, p.heatflux=[p.heatflux; rows(sel,const(f.getString('q0'),tags{q}),tags{q})]; end
 case 'HeatSource'
  % COMSOL 6.1 Heat Source supports either Q0 [W/m^3] or Heat rate
  % P0 [W]. The latter is a selection-total power (verified by integral
  % of ht.Q on TSV PDN4); it is distributed by Domain volume downstream.
  q0=const(f.getString('Q0'),tags{q}); [hasP0,p0]=optionalConst(f,'P0',tags{q});
  if hasP0&&abs(p0)>0
    if abs(q0)>0, error('HeatSource %s has both nonzero Q0 and P0; explicit mode support is required.',tags{q}),end
    p.sourceTotalPower(end+1)=struct('tag',tags{q},'domains',sel,'powerW',p0); %#ok<AGROW>
  else
    p.source=[p.source; rows(sel,q0,tags{q})];
  end
 case 'ThermalInsulation', p.insulation=[p.insulation; string(sel(:))];
 case {'InitialValues','init'}, p.initialTemperature=const(f.getString('Tinit'),tags{q});
 otherwise
  if ~ismember(tags{q},{'solid1','ins1','idi1','ltneb1','os1','cib1','dcont1'}), error('Unsupported material/physics feature tag=%s type=%s selection=%s property=feature',tags{q},typ,mat2str(sel)); end
 end
end
if isempty(p.initialTemperature), error('No supported InitialValues Tinit in heat physics %s.',info.heatPhysics);end
end
function a=rows(ids,v,tag),a=[double(ids(:)),repmat(v,numel(ids),1)];end
function a=rows2(ids,x,y,tag),a=[double(ids(:)),repmat(x,numel(ids),1),repmat(y,numel(ids),1)];end
function x=const(v,tag),tok=regexp(char(v),'^\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)','tokens','once');if isempty(tok),error('unsupported expression in physics feature %s',tag),end,x=str2double(tok{1});end
function [present,x]=optionalConst(f,key,tag),try,x=const(f.getString(key),tag);present=true;catch,present=false;x=0;end,end
