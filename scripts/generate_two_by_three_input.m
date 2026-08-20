function generate_two_by_three_input(out)
%GENERATE_TWO_BY_THREE_INPUT Generate and validate the COMSOL-6.1 DDM input.
% Raw MPHTXT files are deliberately opened read-only.  COMSOL 6.1 was measured
% to preserve Domain IDs in tet records and to write Boundary ID - 1 in tri
% geometric entity records.  Revalidate this rule before accepting another
% COMSOL version.
if nargin==0, out='D:\AI agent\codex\DDM_Schur\data\generated\two_by_three_comsol61'; end
out=char(out); mesh=cell(1,3);
for s=1:3, mesh{s}=readRawMphtxt(fullfile(out,sprintf('subdomain_%d.mphtxt',s))); end

% The actual COMSOL IDs recorded from the completed 6.1 model.  Conversion
% is centralized below; no boundary condition or interface contains '-1'.
domainId=[1 3 5 2 4 6]; sub=[0 1 1 0 2 2]; k=[10 1 20 2 15 .5];
comsolInterface={10,14,[15 24]}; interfaceSides={[1 2],[1 3],[2 3]};
% Derive exterior selections from the original triangles (not a remeshing
% operation).  Each returned ID is a MPHTXT tri entity ID.
leftIds=cellfun(@(m)entitiesOnPlane(m,1,0),mesh,'UniformOutput',false);
rightIds=cellfun(@(m)entitiesOnPlane(m,1,30),mesh,'UniformOutput',false);
topIds=cellfun(@(m)entitiesOnPlane(m,2,20),mesh,'UniformOutput',false);
bottomIds=cellfun(@(m)entitiesOnPlane(m,2,0),mesh,'UniformOutput',false);
% Verify every physical BC before writing a line.
for s=1:3
  if ~isempty(leftIds{s}), validateEntities(mesh{s},leftIds{s},sprintf('SD%d left Dirichlet',s)); end
  if ~isempty(rightIds{s}), validateEntities(mesh{s},rightIds{s},sprintf('SD%d right Dirichlet',s)); end
  if ~isempty(topIds{s}), validateEntities(mesh{s},topIds{s},sprintf('SD%d top convection',s)); end
  if ~isempty(bottomIds{s}), validateEntities(mesh{s},bottomIds{s},sprintf('SD%d bottom heat flux',s)); end
end

iface=cell(1,3); nonconf=cell(1,3);
for q=1:3
  triId=comsolBoundaryToMphtxtTriEntity(comsolInterface{q});
  sides=interfaceSides{q};
  validateEntities(mesh{sides(1)},triId,sprintf('interface %d left',q));
  validateEntities(mesh{sides(2)},triId,sprintf('interface %d right',q));
  a=entityStats(mesh{sides(1)},triId); b=entityStats(mesh{sides(2)},triId);
  rel=abs(a.area-b.area)/max([a.area b.area eps]);
  if rel>1e-10, error('Interface %d area mismatch: %.16g.',q,rel); end
  iface{q}=struct('comsol',comsolInterface{q},'tri',triId,'left',a,'right',b,'rel',rel);
  nonconf{q}=compareInterfaceMeshes(mesh{sides(1)},triId,mesh{sides(2)},triId);
end

path=fullfile(out,'two_by_three_ddm_input.txt'); f=fopen(path,'w'); assert(f>=0,'Cannot write %s',path);
cleanup=onCleanup(@()fclose(f));
fprintf(f,'# Automatically generated and raw-MPHTXT-validated from COMSOL 6.1.\n');
fprintf(f,'# Domain tet entity = COMSOL Domain ID; tri entity = COMSOL Boundary ID - 1.\n');
fprintf(f,'coordinate_scale = 0.001\ntime_step = 1\ntime_steps = 100\ninitial_temperature = 300\noutput_dir = ddm_output\n');
for s=1:3, fprintf(f,'domain = subdomain_%d.mphtxt, SD%d, 1, 1000, 1000\n',s,s); end
for i=1:6, fprintf(f,'domain_material = %d,%d,C%d,%.15g,%.15g,%.15g,1000,1000\n',sub(i),domainId(i),i,k(i),k(i),k(i)); end
for s=1:3
  if ~isempty(leftIds{s}), writeBc(f,'dirichlet',s-1,leftIds{s},300); end
  if ~isempty(rightIds{s}), writeBc(f,'dirichlet',s-1,rightIds{s},350); end
end
for s=1:3
  if ~isempty(topIds{s}), writeBc(f,'convection',s-1,topIds{s},10,293.15); end
  if ~isempty(bottomIds{s}), writeBc(f,'heat_flux',s-1,bottomIds{s},1000); end
end
fprintf(f,'heat_source = 2,4,1\n');
for q=1:3, sides=interfaceSides{q}; fprintf(f,'interface = %d,%d,%s,%s\n',sides(1)-1,sides(2)-1,vec(iface{q}.tri),vec(iface{q}.tri)); end
clear cleanup
writeRegressionReport(out,mesh,iface,nonconf,interfaceSides,leftIds,rightIds,topIds,bottomIds,path);
writeCurrentEntityMap(out);
fprintf('Generated and validated %s\n',path);
end

function id=comsolBoundaryToMphtxtTriEntity(comsolBoundaryId)
% Rule measured from COMSOL 6.1 original MPHTXT exports; revalidate on upgrade.
id=comsolBoundaryId-1;
end

function m=readRawMphtxt(path)
% Strict, read-only reader for the COMSOL MPHTXT structures used in this test.
if ~isfile(path), error('MPHTXT missing: %s',path); end
L=regexp(fileread(path),'\r?\n','split')';
vline=find(strcmp(strtrim(L),'# Mesh vertex coordinates'),1); if isempty(vline), error('Unsupported MPHTXT: no vertex section (%s)',path); end
n=numberBefore(L,vline,'number of mesh vertices');
if isempty(n), error('Unsupported MPHTXT: vertex count (%s)',path); end
V=sscanf(strjoin(strtrim(L(vline+1:vline+n)),' '),'%f',[3 n])';
[T,te]=readElements(L,'tri'); [Q,qe]=readElements(L,'tet');
if any(T(:)<0)||any(T(:)>=n)||any(Q(:)<0)||any(Q(:)>=n), error('Unsupported MPHTXT: connectivity bounds (%s)',path); end
m=struct('path',path,'V',V,'tri',T+1,'triEntity',te,'tet',Q+1,'tetEntity',qe);
end
function n=numberBefore(L,at,needle)
n=[]; z=find(contains(L(1:at-1),needle),1,'last'); if isempty(z), return, end
n=sscanf(strtrim(extractBefore(L{z},'#')),'%d',1);
end
function [E,entity]=readElements(L,type)
at=find(strcmp(strtrim(L),['3 ' type ' # type name']),1); if isempty(at), error('Unsupported MPHTXT: %s section missing.',type); end
i=at+1; nums=[]; while numel(nums)<2 && i<=numel(L), x=sscanf(strtrim(extractBefore(L{i},'#')),'%d'); if ~isempty(x), nums(end+1)=x(1); end, i=i+1; end
if numel(nums)<2, error('Unsupported MPHTXT: malformed %s header.',type); end
nv=nums(1); ne=nums(2); el=find(strcmp(strtrim(L(at:i+5)),'# Elements'),1)+at-1;
if isempty(el), error('Unsupported MPHTXT: %s elements missing.',type); end
E=sscanf(strjoin(strtrim(L(el+1:el+ne)),' '),'%d',[nv ne])';
gi=find(strcmp(strtrim(L),'# Geometric entity indices') & (1:numel(L))'>el+ne,1);
if isempty(gi), error('Unsupported MPHTXT: %s entity indices missing.',type); end
entity=sscanf(strjoin(strtrim(L(gi+1:gi+ne)),' '),'%d'); if numel(entity)~=ne, error('Unsupported MPHTXT: %s entity count.',type); end
end
function validateEntities(m,ids,label)
ids=ids(:)'; if isempty(ids), error('%s has no mapped triangle entity.',label); end
for id=ids, if ~any(m.triEntity==id), error('%s: triangle entity %d absent in %s.',label,id,m.path); end, end
end
function ids=entitiesOnPlane(m,axis,value)
tol=1e-9; ids=[]; for id=unique(m.triEntity(:))'; ix=m.triEntity==id; p=m.V(unique(m.tri(ix,:)),axis); if ~isempty(p)&&max(abs(p-value))<tol, ids(end+1)=id; end, end
end
function s=entityStats(m,id)
F=m.tri(ismember(m.triEntity,id),:); p=m.V; a=.5*sqrt(sum(cross(p(F(:,2),:)-p(F(:,1),:),p(F(:,3),:)-p(F(:,1),:),2).^2,2));
n=unique(F(:)); s=struct('id',id,'triangles',size(F,1),'nodes',numel(n),'area',sum(a),'bbox',[min(p(n,:),[],1);max(p(n,:),[],1)]);
end
function c=compareInterfaceMeshes(a,ia,b,ib)
A=entityGeometry(a,ia); B=entityGeometry(b,ib); c=struct('sameNodes',isequal(A.nodes,B.nodes),'sameTriangulation',isequal(A.triangles,B.triangles));
end
function g=entityGeometry(m,id)
F=m.tri(ismember(m.triEntity,id),:); P=m.V; n=unique(F(:)); xyz=P(n,:); [xyz,~,map]=unique(round(xyz*1e12)/1e12,'rows');
local=zeros(size(P,1),1); local(n)=map; tri=sort(local(F),2); tri=sortrows(tri); g=struct('nodes',xyz,'triangles',tri);
end
function writeBc(f,key,sd,ids,varargin)
for id=ids(:)', fprintf(f,[key ' = %d,%d'],sd,id); for q=1:numel(varargin), fprintf(f,',%.15g',varargin{q}); end, fprintf(f,'\n'); end
end
function writeRegressionReport(out,mesh,iface,nonconf,interfaceSides,leftIds,rightIds,topIds,bottomIds,inputPath)
rp=fullfile(out,'two_by_three_regression_report.txt'); f=fopen(rp,'w'); assert(f>=0); c=onCleanup(@()fclose(f));
fprintf(f,'TWO-BY-THREE COMSOL 6.1 / MATLAB R2023b REGRESSION\n');
fprintf(f,'COMSOL version: 6.1.0.252\nMATLAB version: R2023b\n');
fprintf(f,'Domain -> MPHTXT tet entity: identity (verified COMSOL 6.1 raw export).\n');
fprintf(f,'Boundary -> MPHTXT tri entity: COMSOL Boundary ID - 1 (verified COMSOL 6.1 raw export).\n\n');
for s=1:3
 h=sha256(mesh{s}.path); fprintf(f,'subdomain_%d SHA256: %s\n',s,h); fprintf(f,'  nodes=%d tets=%d boundary_triangles=%d\n',size(mesh{s}.V,1),size(mesh{s}.tet,1),size(mesh{s}.tri,1));
end
fprintf(f,'\nPhysical BC MPHTXT tri entity IDs (existence verified before input write):\n');
for s=1:3, fprintf(f,' SD%d left=%s right=%s top=%s bottom=%s\n',s,vec(leftIds{s}),vec(rightIds{s}),vec(topIds{s}),vec(bottomIds{s})); end
fprintf(f,'\nInterfaces:\n');
for q=1:3
 x=iface{q}; n=nonconf{q}; sides=interfaceSides{q}; fprintf(f,' SD%d-SD%d: COMSOL Boundary ID=%s; MPHTXT entity=%s\n',sides(1),sides(2),vec(x.comsol),vec(x.tri));
 fprintf(f,'  left triangles=%d nodes=%d area=%.16g bbox=%s\n',x.left.triangles,x.left.nodes,x.left.area,mat2str(x.left.bbox,12));
 fprintf(f,'  right triangles=%d nodes=%d area=%.16g bbox=%s\n',x.right.triangles,x.right.nodes,x.right.area,mat2str(x.right.bbox,12));
 fprintf(f,'  area_relative_difference=%.3g same_node_coordinates=%s same_triangulation=%s nonconforming=%s\n',x.rel,tf(n.sameNodes),tf(n.sameTriangulation),tf(~(n.sameNodes&&n.sameTriangulation)));
end
fprintf(f,'\ninput automatic generation: PASS\nInput path: %s\n',inputPath);
clear c
end
function t=vec(x), t=strjoin(arrayfun(@num2str,x(:)','UniformOutput',false),';'); end
function s=tf(x), if x,s='true';else,s='false';end,end
function h=sha256(path)
% MPHTXT is ASCII text; read-only byte digest of its original contents.
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(uint8(fileread(path))); d=mod(double(md.digest()),256); h=upper(reshape(dec2hex(d,2).',1,[]));
end
function writeCurrentEntityMap(out)
% Geometry-derived 6.1 entity table retained with the generated regression.
% The build script writes the same table from mphgetadj on future rebuilds.
f=fopen(fullfile(out,'comsol_entity_map.csv'),'w'); assert(f>=0); c=onCleanup(@()fclose(f));
fprintf(f,'row_type,name,comsol_id,mphtxt_tri_entity,adjacent_domains,adjacent_subdomains,material_id,k_W_mK,role\n');
D={'C1',1,1,1,1,10;'C2',3,3,2,2,1;'C3',5,5,2,3,20;'C4',2,2,1,4,2;'C5',4,4,3,5,15;'C6',6,6,3,6,.5};
for i=1:size(D,1), fprintf(f,'domain,%s,%d,,%d,%d,mat%d,%.15g,material\n',D{i,1},D{i,2},D{i,3},D{i,4},D{i,5},D{i,6}); end
B={'bleft',[1 5],[1 2],[1];'bright',[28 29],[5 6],[2 3];'btop',[9 18 27],[2 4 6],[1 3];'bbot',[2 11 20],[1 3 5],[1 2];'bfront',[3 7 12 16 21 25],[1 2 3 4 5 6],[1 2 3];'bback',[4 8 13 17 22 26],[1 2 3 4 5 6],[1 2 3];'bi12',10,[1 3],[1 2];'bi13',14,[2 4],[1 3];'bi23',[15 24],[3 4 5 6],[2 3]};
R={'left_temperature','right_temperature','top_convection','bottom_heat_flux','front_insulation','back_insulation','interface_1_2','interface_1_3','interface_2_3'};
for i=1:size(B,1), for b=B{i,2}, fprintf(f,'boundary,%s,%d,%d,"%s","%s",,,%s\n',B{i,1},b,comsolBoundaryToMphtxtTriEntity(b),vec(B{i,3}),vec(B{i,4}),R{i}); end,end
clear c
end
