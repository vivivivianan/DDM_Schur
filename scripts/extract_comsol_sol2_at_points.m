function extract_comsol_sol2_at_points(mphFile,pointTxt,outputCsv,summaryFile)
% Read-only comparison of a saved COMSOL sol2 field at supplied SI points.
% pointTxt must contain x_m y_m z_m temperature_K with no header.
import com.comsol.model.util.*

m=mphload(mphFile,['ddm_point_compare_' regexprep(char(java.util.UUID.randomUUID),'-','')]);
cleanup=onCleanup(@()ModelUtil.remove(m.tag));
sol=mphsolinfo(m,'soltag','sol2');
if isempty(sol.solvals), error('sol2 has no saved solution times.'); end
native=mpheval(m,'T','edim','domain','dataset','dset2','solnum','end');
nativeT=double(native.d1(:));
points=readmatrix(pointTxt,'FileType','text');
if size(points,2)~=4 || any(~isfinite(points(:)))
    error('Expected finite four-column x_m y_m z_m temperature_K data.');
end

n=size(points,1); chunk=25000; comsolT=nan(n,1); validPoint=false(n,1);
fid=fopen(outputCsv,'w');
if fid<0, error('Cannot write %s',outputCsv); end
fprintf(fid,'x_m,y_m,z_m,ddm_m3_temperature_k,comsol_temperature_k,error_ddm_minus_comsol_k,valid_comsol_point\n');
for first=1:chunk:n
    last=min(n,first+chunk-1); idx=first:last;
    % ext=0 disables extrapolation beyond the actual COMSOL mesh.
    values=double(mphinterp(m,'T','coord',points(idx,1:3).','dataset','dset2','solnum','end','ext',0));
    values=values(:);
    if numel(values)~=numel(idx), fclose(fid); error('COMSOL interpolation size mismatch for rows %d through %d.',first,last); end
    comsolT(idx)=values;
    validPoint(idx)=isfinite(values);
    rows=[points(idx,1:3),points(idx,4),values,points(idx,4)-values,double(validPoint(idx))];
    fprintf(fid,'%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d\n',rows.');
end
fclose(fid);
if ~any(validPoint), error('No supplied DDM point lies in the COMSOL mesh with extrapolation disabled.'); end

valid=find(validPoint); err=points(valid,4)-comsolT(valid);
[ddmMin,ddmMinIndex]=min(points(:,4)); [ddmMax,ddmMaxIndex]=max(points(:,4));
[sampleMin,localSampleMinIndex]=min(comsolT(valid)); sampleMinIndex=valid(localSampleMinIndex);
[sampleMax,localSampleMaxIndex]=max(comsolT(valid)); sampleMaxIndex=valid(localSampleMaxIndex);
[maxAbs,localMaxAbsIndex]=max(abs(err)); maxAbsIndex=valid(localMaxAbsIndex);
fid=fopen(summaryFile,'w');
if fid<0, error('Cannot write %s',summaryFile); end
fprintf(fid,['SOURCE=%s\nSOLTAG=sol2\nTIME_END_S=%.17g\nPOINT_COUNT=%d\nVALID_POINT_COUNT=%d\nINVALID_POINT_COUNT=%d\n' ...
    'COMSOL_NATIVE_TMIN_K=%.17g\nCOMSOL_NATIVE_TMAX_K=%.17g\n' ...
    'DDM_TMIN_K=%.17g\nDDM_TMIN_XYZ_M=%.17g %.17g %.17g\n' ...
    'DDM_TMAX_K=%.17g\nDDM_TMAX_XYZ_M=%.17g %.17g %.17g\n' ...
    'COMSOL_AT_DDM_POINTS_TMIN_K=%.17g\nCOMSOL_AT_DDM_POINTS_TMIN_XYZ_M=%.17g %.17g %.17g\n' ...
    'COMSOL_AT_DDM_POINTS_TMAX_K=%.17g\nCOMSOL_AT_DDM_POINTS_TMAX_XYZ_M=%.17g %.17g %.17g\n' ...
    'TMIN_DIFFERENCE_DDM_MINUS_COMSOL_SAMPLED_K=%.17g\n' ...
    'ERROR_RMS_K=%.17g\nERROR_MEAN_K=%.17g\nERROR_MAX_ABSOLUTE_K=%.17g\n' ...
    'ERROR_MAX_ABSOLUTE_XYZ_M=%.17g %.17g %.17g\nERROR_AT_MAX_ABSOLUTE_K=%.17g\n'], ...
    mphFile,double(sol.solvals(end)),n,numel(valid),n-numel(valid),min(nativeT),max(nativeT), ...
    ddmMin,points(ddmMinIndex,1),points(ddmMinIndex,2),points(ddmMinIndex,3), ...
    ddmMax,points(ddmMaxIndex,1),points(ddmMaxIndex,2),points(ddmMaxIndex,3), ...
    sampleMin,points(sampleMinIndex,1),points(sampleMinIndex,2),points(sampleMinIndex,3), ...
    sampleMax,points(sampleMaxIndex,1),points(sampleMaxIndex,2),points(sampleMaxIndex,3), ...
    ddmMin-sampleMin,sqrt(mean(err.^2)),mean(err),maxAbs, ...
    points(maxAbsIndex,1),points(maxAbsIndex,2),points(maxAbsIndex,3),err(maxAbsIndex));
fclose(fid);
clear cleanup
end
