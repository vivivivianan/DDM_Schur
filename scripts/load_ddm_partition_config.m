function cfg=load_ddm_partition_config(path)
% Load DDM run options and the optional legacy Domain-ID partition.
if ~isfile(path), error('Partition configuration does not exist: %s',path); end
cfg=jsondecode(fileread(path));
if ~isfield(cfg,'subdomains'), cfg.subdomains=struct([]); end
if ~isfield(cfg,'solver'), cfg.solver=struct(); end
if ~isfield(cfg,'component_tag'), cfg.component_tag=''; end
if ~isfield(cfg,'geometry_tag'), cfg.geometry_tag=''; end
if ~isfield(cfg,'heat_physics_tag'), cfg.heat_physics_tag=''; end
if ~isfield(cfg,'default_mesh_hmax'), cfg.default_mesh_hmax=0.002; end
if ~isscalar(cfg.default_mesh_hmax)||cfg.default_mesh_hmax<=0, error('default_mesh_hmax must be positive metres.'); end
if ~isfield(cfg.solver,'dt'), cfg.solver.dt=1.0; end
if ~isfield(cfg.solver,'t_end'), cfg.solver.t_end=100.0; end
if ~isempty(cfg.subdomains), ids=[cfg.subdomains.id]; if numel(unique(ids))~=numel(ids), error('Subdomain IDs must be unique.'); end, end
for k=1:numel(cfg.subdomains)
 s=cfg.subdomains(k);
 if ~isfield(s,'domain_ids'), cfg.subdomains(k).domain_ids=[]; end
 if ~isfield(s,'mesh')||~isfield(s.mesh,'hmax'), cfg.subdomains(k).mesh=struct('hmax',cfg.default_mesh_hmax); end
 if cfg.subdomains(k).mesh.hmax<=0, error('Subdomain %g mesh.hmax must be positive metres.',s.id); end
end
if ~isfield(cfg.solver,'dt')||~isfield(cfg.solver,'t_end')||cfg.solver.dt<=0||cfg.solver.t_end<=0, error('solver.dt and solver.t_end must be positive.'); end
end
