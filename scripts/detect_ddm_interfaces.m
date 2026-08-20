function interfaces=detect_ddm_interfaces(info,cfg)
% Domain-boundary adjacency + config mapping only; never coordinate heuristics.
owner=containers.Map('KeyType','double','ValueType','double'); for s=1:numel(cfg.subdomains),for d=cfg.subdomains(s).domain_ids(:)',owner(double(d))=s;end,end
interfaces=struct('left',{},'right',{},'comsolBoundaryIds',{}); pairs=containers.Map;
for b=info.boundaries
 d=info.boundaryDomains{b}; if numel(d)~=2||~isKey(owner,d(1))||~isKey(owner,d(2)),continue,end; a=owner(d(1));z=owner(d(2));if a==z,continue,end; key=sprintf('%d_%d',min(a,z),max(a,z));if ~isKey(pairs,key),pairs(key)=b;else,pairs(key)=[pairs(key),b];end
end
K=keys(pairs);for q=1:numel(K),x=sscanf(K{q},'%d_%d');interfaces(end+1)=struct('left',x(1),'right',x(2),'comsolBoundaryIds',pairs(K{q}));end
end
