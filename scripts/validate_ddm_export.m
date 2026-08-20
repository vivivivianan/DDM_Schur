function validation=validate_ddm_export(meshes,interfaces,info,cfg)
% Validate conversion existence, per-entity geometry, and conformality.
for s=1:numel(meshes), want=sort(cfg.subdomains(s).domain_ids(:)); got=sort(unique(meshes{s}.tetEntity));if ~isequal(want,got),error('Subdomain %d MPHTXT domains mismatch: expected %s got %s',s,mat2str(want),mat2str(got));end,end
validation=struct; validation.interfaces=struct('left',{},'right',{},'comsol',{},'tri',{},'A',{},'B',{},'sameNodes',{},'sameTriangulation',{}); for q=1:numel(interfaces)
 i=interfaces(q); ids=comsolBoundaryToMphtxtTriEntity(i.comsolBoundaryIds); A=stats(meshes{i.left},ids);B=stats(meshes{i.right},ids);if A.triangles==0||B.triangles==0||A.area<=0||B.area<=0,error('Interface SD%d-SD%d does not map to raw MPHTXT triangles.',i.left,i.right),end; rel=abs(A.area-B.area)/max([A.area B.area eps]);if rel>1e-9,error('Interface SD%d-SD%d area mismatch %g.',i.left,i.right,rel),end; validation.interfaces(end+1)=struct('left',i.left,'right',i.right,'comsol',i.comsolBoundaryIds,'tri',ids,'A',A,'B',B,'sameNodes',sameNodes(meshes{i.left},ids,meshes{i.right},ids),'sameTriangulation',sameTriangles(meshes{i.left},ids,meshes{i.right},ids));end
end
function id=comsolBoundaryToMphtxtTriEntity(id),id=id-1;end
function s=stats(m,id),F=m.tri(ismember(m.triEntity,id),:);p=m.V;n=unique(F(:));a=.5*vecnorm(cross(p(F(:,2),:)-p(F(:,1),:),p(F(:,3),:)-p(F(:,1),:),2),2,2);s=struct('triangles',size(F,1),'nodes',numel(n),'area',sum(a),'bbox',[min(p(n,:),[],1);max(p(n,:),[],1)]);end
function x=sameNodes(a,ia,b,ib),x=isequal(sortrows(unique(round(a.V(unique(a.tri(ismember(a.triEntity,ia),:)),:)*1e12)/1e12,'rows')),sortrows(unique(round(b.V(unique(b.tri(ismember(b.triEntity,ib),:)),:)*1e12)/1e12,'rows')));end
function x=sameTriangles(a,ia,b,ib)
% Compare geometric triangles, not local node numbering. This is diagnostic
% only; nonconforming interfaces are valid and expected by the DDM solver.
if ~sameNodes(a,ia,b,ib), x=false; return, end
x=isequal(canonicalTriangles(a,ia),canonicalTriangles(b,ib));
end
function C=canonicalTriangles(m,ids)
F=m.tri(ismember(m.triEntity,ids),:); P=round(m.V*1e12)/1e12;
C=zeros(size(F,1),9);
for q=1:size(F,1), v=sortrows(P(F(q,:),:)); C(q,:)=reshape(v.',1,[]); end
C=sortrows(C);
end
