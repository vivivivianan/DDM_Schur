function audit_mphtxt_internal_connectivity(meshFiles, outputFile)
% Read-only audit: count cross-COMSOL-domain tetrahedral shared faces.
% A zero count for expected internal material interfaces would reveal that
% the raw mesh topology cannot support continuous intra-subdomain conduction.
fid=fopen(outputFile,'w'); if fid<0, error('Cannot write %s',outputFile); end
closeFile=onCleanup(@()fclose(fid));
for s=1:numel(meshFiles)
 m=parse_mphtxt_entities(meshFiles{s}); F=tetFaces(m.tet); owner=repmat(m.tetEntity,4,1);
 [Fs,ix]=sort(F,2); [~,ord]=sortrows(Fs); Fs=Fs(ord,:); owner=owner(ord);
 sameFace=all(Fs(2:end,:)==Fs(1:end-1,:),2);
 cross=find(sameFace & owner(2:end)~=owner(1:end-1));
 fprintf(fid,'SUBDOMAIN=%d\n',s);
 fprintf(fid,'TETS=%d\n',size(m.tet,1));
 fprintf(fid,'CROSS_DOMAIN_SHARED_FACES=%d\n',numel(cross));
 if isempty(cross), fprintf(fid,'CROSS_DOMAIN_PAIRS=\n');
 else
  pairs=sort([owner(cross),owner(cross+1)],2);
  [u,~,g]=unique(pairs,'rows'); counts=accumarray(g,1);
  for q=1:size(u,1), fprintf(fid,'PAIR=%d,%d SHARED_FACES=%d\n',u(q,1),u(q,2),counts(q)); end
 end
end
clear closeFile
end

function F=tetFaces(T)
F=[T(:,[1 2 3]);T(:,[1 2 4]);T(:,[1 3 4]);T(:,[2 3 4])];
end
