function cfg=resolve_ddm_subdomains_from_selections(model,cfg,info)
% Resolve DDM domains from component Explicit selections named ddm_sd_###.
% Named selections take priority over legacy JSON domain_ids.  The selection
% entities are COMSOL Domain IDs; no mesh entity numbering is involved here.
c=model.component(info.component); tags=cellstr(char(c.selection.tags)); found=struct('id',{},'domains',{},'tag',{},'label',{});
for q=1:numel(tags)
  tag=tags{q}; sel=c.selection(tag); label=selectionLabel(sel,tag);
  tagId=parseId(tag); labelId=parseId(label);
  if isempty(tagId)&&isempty(labelId), continue, end
  if ~isempty(tagId)&&~isempty(labelId)&&tagId~=labelId
    error('DDM selection tag "%s" and label "%s" specify different IDs.',tag,label);
  end
  id=tagId; if isempty(id), id=labelId; end
  typ=selectionType(sel,tag);
  if ~strcmpi(typ,'Explicit'), error('DDM selection %s must be Explicit, got %s.',tag,typ); end
  S=mphgetselection(sel);
  if ~isfield(S,'dimension')||double(S.dimension)~=3
    error('DDM selection %s must be a Domain-level (dimension 3) selection.',tag);
  end
  domains=double(S.entities(:))';
  if isempty(domains), error('DDM selection %s is empty.',tag); end
  if any(~ismember(domains,info.domains)), error('DDM selection %s contains invalid Domain IDs %s.',tag,mat2str(domains)); end
  found(end+1)=struct('id',id,'domains',sort(domains),'tag',tag,'label',label); %#ok<AGROW>
end
if isempty(found)
  if isempty(cfg.subdomains), error('No ddm_sd_### Explicit Domain selections found and JSON has no legacy subdomains/domain_ids.'); end
  for q=1:numel(cfg.subdomains)
    if isempty(cfg.subdomains(q).domain_ids), error('No named DDM selections: legacy JSON Subdomain %d needs domain_ids.',cfg.subdomains(q).id); end
  end
  cfg.partition_source='json_domain_ids'; return
end
ids=[found.id]; if numel(unique(ids))~=numel(ids), error('Two ddm_sd_* selections map to the same Subdomain ID.'); end
found=found(orderById(found));
% If the legacy mapping is present too, it is a consistency assertion only.
for q=1:numel(cfg.subdomains)
  js=cfg.subdomains(q);
  ix=find([found.id]==js.id,1);
  if isempty(ix), error('JSON subdomain %d has no matching ddm_sd_* selection.',js.id); end
  if ~isempty(js.domain_ids)&&~isequal(sort(double(js.domain_ids(:))'),found(ix).domains)
    error('JSON domain_ids for Subdomain %d disagree with selection %s.',js.id,found(ix).tag);
  end
end
new=struct('id',{},'domain_ids',{},'mesh',{});
for q=1:numel(found)
  ix=find(arrayfun(@(x)x.id==found(q).id,cfg.subdomains),1);
  h=cfg.default_mesh_hmax;
  if ~isempty(ix), h=cfg.subdomains(ix).mesh.hmax; end
  new(end+1)=struct('id',found(q).id,'domain_ids',found(q).domains,'mesh',struct('hmax',h)); %#ok<AGROW>
end
allDomains=[found.domains];
if numel(unique(allDomains))~=numel(allDomains), error('A COMSOL Domain occurs in more than one ddm_sd_* selection.'); end
if ~isequal(sort(allDomains),sort(info.heatDomains))
  error('ddm_sd_* selections must cover each Heat Transfer Domain exactly once. selections=%s heat_domains=%s',mat2str(sort(allDomains)),mat2str(sort(info.heatDomains)));
end
cfg.subdomains=new; cfg.partition_source='comsol_named_selections';
end
function out=orderById(x),[~,out]=sort([x.id]);end
function id=parseId(value)
m=regexp(char(value),'^ddm_sd_(\d+)$','tokens','once','ignorecase'); if isempty(m),id=[];else,id=str2double(m{1});end
end
function value=selectionLabel(sel,fallback)
try, value=char(sel.label()); catch, try, value=char(sel.label); catch, value=fallback; end, end
end
function value=selectionType(sel,tag)
try, value=char(sel.getType()); catch err, error('Could not read selection type for %s: %s',tag,err.message); end
end
