function info=inspect_comsol_heat_model(model,cfg)
% Read model entities/features; never infer physics from coordinates.
ct=cellstr(char(model.component.tags)); if isempty(cfg.component_tag), if numel(ct)~=1,error('Ambiguous component; set component_tag.'),end, comp=ct{1}; else,comp=char(cfg.component_tag);end
c=model.component(comp); gt=cellstr(char(c.geom.tags)); if isempty(cfg.geometry_tag), if numel(gt)~=1,error('Ambiguous geometry; set geometry_tag.'),end,geom=gt{1};else,geom=char(cfg.geometry_tag);end
pt=cellstr(char(c.physics.tags)); if isempty(cfg.heat_physics_tag), candidates={}; for q=1:numel(pt), if contains(lower(char(c.physics(pt{q}).getType)),'heat'),candidates{end+1}=pt{q};end,end; if numel(candidates)~=1,error('Ambiguous/missing heat physics; set heat_physics_tag.'),end,ht=candidates{1};else,ht=char(cfg.heat_physics_tag);end
g=c.geom(geom); nd=double(g.getNDomains); nb=double(g.getNBoundaries);
% COMSOL MPHTXT coordinates are expressed in the geometry length unit.
% Convert that unit once, here, for the C++ solver's SI-coordinate input.
% Unknown units deliberately fail instead of silently producing bad physics.
lengthUnit=char(g.lengthUnit()); coordinateScale=geometryLengthUnitToSI(lengthUnit);
info=struct('component',comp,'geometry',geom,'heatPhysics',ht,'domains',1:nd,'boundaries',1:nb,'comsolVersion','6.1 (reported by LiveLink runtime)', 'geometryLengthUnit',lengthUnit,'coordinateScale',coordinateScale);
info.boundaryDomains=cell(1,nb); for b=1:nb, info.boundaryDomains{b}=double(mphgetadj(model,geom,'domain','boundary',b))'; end
% Heat Transfer domains come from the Solid Heat Transfer feature.  An empty
% feature selection means COMSOL's default: every geometry Domain.
h=c.physics(ht); heatDomains=[]; ft=cellstr(char(h.feature.tags));
for q=1:numel(ft)
 f=h.feature(ft{q}); if contains(lower(char(f.getType)),'solidheattransfermodel')
   d=double(f.selection.entities)'; if isempty(d), d=info.domains; end; heatDomains=[heatDomains,d]; %#ok<AGROW>
 end
end
if isempty(heatDomains), error('Could not determine Heat Transfer Domains from physics %s.',ht); end
info.heatDomains=unique(heatDomains);
info.materials=extract_comsol_materials(model,info); info.physics=extract_comsol_heat_physics(model,info);
end

function scale=geometryLengthUnitToSI(unit)
% Stable COMSOL 6.1/6.4 geometry length-unit mapping.  MPHTXT itself is
% retained byte-for-byte; this scale is only written into the DDM input.
switch lower(strtrim(unit))
 case {'m','meter','metre'}, scale=1;
 case {'dm'}, scale=1e-1;
 case {'cm'}, scale=1e-2;
 case {'mm'}, scale=1e-3;
 case {'um','µm','μm'}, scale=1e-6;
 case {'nm'}, scale=1e-9;
 otherwise, error('Unsupported COMSOL geometry length unit "%s". Add an explicit validated SI conversion before export.',unit)
end
end
