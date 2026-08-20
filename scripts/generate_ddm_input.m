function inputPath=generate_ddm_input(out,cfg,info,meshes,interfaces,validation)
% Emit only the existing C++ input grammar, after raw-MPHTXT validation.
integrator='backward_euler';if isfield(cfg.solver,'time_integrator'),integrator=cfg.solver.time_integrator;end
inputPath=fullfile(out,'ddm_input.txt'); f=fopen(inputPath,'w');if f<0,error('Cannot write %s',inputPath),end;c=onCleanup(@()fclose(f));fprintf(f,'# Generated from COMSOL selections and raw MPHTXT entity validation.\n# COMSOL geometry length unit: %s; MPHTXT coordinates -> SI scale.\ncoordinate_scale = %.17g\ntime_integrator = %s\ntime_step = %.17g\ntime_steps = %d\ninitial_temperature = %.17g\noutput_dir = ddm_output\n',info.geometryLengthUnit,info.coordinateScale,integrator,cfg.solver.dt,round(cfg.solver.t_end/cfg.solver.dt),info.physics.initialTemperature);
for s=1:numel(cfg.subdomains),fprintf(f,'domain = raw_mesh/subdomain_%03d.mphtxt, SD%d, 1, 1000, 1000\n',s,cfg.subdomains(s).id);end
for d=info.heatDomains, mat=info.materials(arrayfun(@(x)ismember(d,x.domains),info.materials));s=find(arrayfun(@(x)ismember(d,x.domain_ids),cfg.subdomains),1);fprintf(f,'domain_material = %d,%d,%s,%.17g,%.17g,%.17g,%.17g,%.17g\n',s-1,d,mat.tag,mat.k,mat.k,mat.k,mat.rho,mat.cp);end
writeBoundary(f,'dirichlet',info.physics.dirichlet,info,cfg,meshes);writeBoundary(f,'convection',info.physics.convection,info,cfg,meshes);writeBoundary(f,'heat_flux',info.physics.heatflux,info,cfg,meshes);
for q=1:size(info.physics.source,1),d=info.physics.source(q,1);s=find(arrayfun(@(x)ismember(d,x.domain_ids),cfg.subdomains),1);if isempty(s),error('Heat source domain %d is not partitioned.',d),end;fprintf(f,'heat_source = %d,%d,%.17g\n',s-1,d,info.physics.source(q,2)*domainVolume(meshes{s},d,info.coordinateScale));end
for q=1:numel(info.physics.sourceTotalPower)
 src=info.physics.sourceTotalPower(q); volumes=zeros(size(src.domains)); owners=zeros(size(src.domains));
 for j=1:numel(src.domains),d=src.domains(j);owner=find(arrayfun(@(x)ismember(d,x.domain_ids),cfg.subdomains),1);if isempty(owner),error('Heat-rate source domain %d is not partitioned.',d),end;owners(j)=owner;volumes(j)=domainVolume(meshes{owner},d,info.coordinateScale);end
 totalVolume=sum(volumes);if totalVolume<=0,error('Heat-rate source %s has zero selected volume.',src.tag),end
 for j=1:numel(src.domains),fprintf(f,'heat_source = %d,%d,%.17g\n',owners(j)-1,src.domains(j),src.powerW*volumes(j)/totalVolume);end
end
for q=1:numel(validation.interfaces),x=validation.interfaces(q);fprintf(f,'interface = %d,%d,%s,%s\n',x.left-1,x.right-1,list(x.tri),list(x.tri));end
end
function writeBoundary(f,key,a,info,cfg,meshes)
for q=1:size(a,1),b=a(q,1);ad=info.boundaryDomains{b};if numel(ad)~=1,error('%s boundary %d is not an external single-domain COMSOL selection.',key,b),end;s=find(arrayfun(@(x)ismember(ad,x.domain_ids),cfg.subdomains),1);if isempty(s),error('%s boundary %d belongs to an unpartitioned domain.',key,b),end;e=comsolBoundaryToMphtxtTriEntity(b);if ~any(meshes{s}.triEntity==e),error('%s COMSOL boundary %d converts to absent MPHTXT entity %d.',key,b,e),end;fprintf(f,[key ' = %d,%d'],s-1,e);for k=2:size(a,2),fprintf(f,',%.17g',a(q,k));end;fprintf(f,'\n');end
end
function id=comsolBoundaryToMphtxtTriEntity(id),id=id-1;end
function v=domainVolume(m,d,coordinateScale),F=m.tet(m.tetEntity==d,:);p=m.V*coordinateScale;v=sum(abs(dot(p(F(:,2),:)-p(F(:,1),:),cross(p(F(:,3),:)-p(F(:,1),:),p(F(:,4),:)-p(F(:,1),:),2),2))/6);end
function s=list(x),s=strjoin(arrayfun(@num2str,x,'UniformOutput',false),';');end
