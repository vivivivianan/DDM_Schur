function materials=extract_comsol_materials(model,info)
c=model.component(info.component); tags=cellstr(char(c.material.tags)); materials=struct('tag',{},'domains',{},'k',{},'rho',{},'cp',{});
for q=1:numel(tags)
 m=c.material(tags{q}); d=double(m.selection.entities)'; if isempty(d), continue,end
 try, pg=m.propertyGroup('def'); k=parseConst(pg.getString('thermalconductivity'),'thermalconductivity',tags{q}); rho=parseConst(pg.getString('density'),'density',tags{q}); cp=parseConst(pg.getString('heatcapacity'),'heatcapacity',tags{q}); catch ME, error('Unsupported material/physics feature tag=%s type=material selection=%s property=%s',tags{q},mat2str(d),ME.message); end
 if numel(k)~=1 && ~(numel(k)==3&&max(abs(k-k(1)))<1e-12), error('Unsupported material/physics feature tag=%s type=anisotropic selection=%s property=thermalconductivity',tags{q},mat2str(d)); end
 materials(end+1)=struct('tag',tags{q},'domains',d,'k',k(1),'rho',rho(1),'cp',cp(1)); %#ok<AGROW>
end
for d=info.domains
 n=sum(arrayfun(@(x)ismember(d,x.domains),materials)); if n~=1,error('Domain %d must have exactly one supported material; found %d.',d,n);end
end
end
function x=parseConst(v,property,tag)
if iscell(v), v=char(v{1}); else,v=char(v);end
tok=regexp(v,'^\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)','tokens','once'); if isempty(tok),error('unsupported expression for %s in %s',property,tag),end; x=str2double(tok{1});
end
