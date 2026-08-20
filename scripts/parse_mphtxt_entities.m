function m=parse_mphtxt_entities(path)
% Strict read-only parser. Unknown/malformed COMSOL MPHTXT fails explicitly.
if ~isfile(path),error('MPHTXT missing: %s',path),end; L=regexp(fileread(path),'\r?\n','split')'; v=find(strcmp(strtrim(L),'# Mesh vertex coordinates'),1); if isempty(v),error('Unsupported MPHTXT: vertex section missing'),end
n=sscanf(strtrim(extractBefore(L{find(contains(L(1:v),'number of mesh vertices'),1,'last')},'#')),'%d',1); V=sscanf(strjoin(strtrim(L(v+1:v+n)),' '),'%f',[3 n])'; [T,te]=elements(L,'tri');[Q,qe]=elements(L,'tet'); if any(T(:)<0)||any(T(:)>=n)||any(Q(:)<0)||any(Q(:)>=n),error('Unsupported MPHTXT: invalid connectivity'),end
m=struct('path',path,'V',V,'tri',T+1,'triEntity',te,'tet',Q+1,'tetEntity',qe);
end
function [E,G]=elements(L,t)
a=find(strcmp(strtrim(L),['3 ' t ' # type name']),1);if isempty(a),error('Unsupported MPHTXT: %s missing',t),end;i=a+1;x=[];while numel(x)<2,z=sscanf(strtrim(extractBefore(L{i},'#')),'%d');if ~isempty(z),x(end+1)=z(1);end;i=i+1;end
nv=x(1);ne=x(2);e=find(strcmp(strtrim(L),'# Elements')&(1:numel(L))'>a,1);E=sscanf(strjoin(strtrim(L(e+1:e+ne)),' '),'%d',[nv ne])';g=find(strcmp(strtrim(L),'# Geometric entity indices')&(1:numel(L))'>e+ne,1);if isempty(g),error('Unsupported MPHTXT: geometric entity indices missing'),end;G=sscanf(strjoin(strtrim(L(g+1:g+ne)),' '),'%d');if numel(G)~=ne,error('Unsupported MPHTXT: entity index count'),end
end
