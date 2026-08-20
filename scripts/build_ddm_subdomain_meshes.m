function paths=build_ddm_subdomain_meshes(model,cfg,info,out)
% Create independent in-memory mesh sequences and export COMSOL-native MPHTXT.
% local_refinement is applied only to these new in-memory sequences; the
% source MPH, physics, selections and existing mesh sequences are untouched.
raw=fullfile(out,'raw_mesh'); if ~isfolder(raw),mkdir(raw),end; c=model.component(info.component); paths=cell(1,numel(cfg.subdomains));
for s=1:numel(cfg.subdomains)
 tag=sprintf('ddm_mesh_%03d',s); if any(strcmp(cellstr(char(c.mesh.tags)),tag)), c.mesh.remove(tag); end
 ms=c.mesh.create(tag); ms.create('size1','Size'); ms.feature('size1').selection.geom(info.geometry,3); ms.feature('size1').selection.set(cfg.subdomains(s).domain_ids); ms.feature('size1').set('custom','on'); ms.feature('size1').set('hmax',sprintf('%.17g[m]',cfg.subdomains(s).mesh.hmax)); ms.feature('size1').set('hmin',sprintf('%.17g[m]',cfg.subdomains(s).mesh.hmax/4));
 if isfield(cfg,'local_refinement') && isfield(cfg.local_refinement,'subdomain') && cfg.local_refinement.subdomain==s
  lr=cfg.local_refinement;
  if ~isfield(lr,'domain_ids') || isempty(lr.domain_ids) || ~isfield(lr,'hmax') || ~isscalar(lr.hmax) || lr.hmax<=0, error('local_refinement requires nonempty domain_ids and positive hmax.'); end
  if ~all(ismember(lr.domain_ids,cfg.subdomains(s).domain_ids)), error('local_refinement Domain IDs must belong to its target subdomain.'); end
  ms.create('size_local','Size'); ms.feature('size_local').selection.geom(info.geometry,3); ms.feature('size_local').selection.set(lr.domain_ids); ms.feature('size_local').set('custom','on'); ms.feature('size_local').set('hmax',sprintf('%.17g[m]',lr.hmax)); ms.feature('size_local').set('hmin',sprintf('%.17g[m]',lr.hmax/4));
  fprintf('SD%d local COMSOL Size refinement: Domains=%s, hmax=%.17g m\n',s,mat2str(lr.domain_ids),lr.hmax);
 end
 ms.create('ftet1','FreeTet'); ms.feature('ftet1').selection.geom(info.geometry,3); ms.feature('ftet1').selection.set(cfg.subdomains(s).domain_ids); ms.run;
 paths{s}=fullfile(raw,sprintf('subdomain_%03d.mphtxt',s)); ms.export(paths{s});
end
end
