function build_two_by_three_comsol(outputDir,mliDir,port)
if nargin<1, outputDir='D:\AI agent\codex\DDM_Schur\data\generated\two_by_three_comsol61'; end
if nargin<2, mliDir='D:\Program Files\COMSOL\COMSOL61\Multiphysics\mli'; end
if nargin<3, port=20361; end
if ~isfolder(outputDir), mkdir(outputDir); end
addpath(mliDir); mphstart('localhost',port); import com.comsol.model.*; import com.comsol.model.util.*
model=ModelUtil.create('two_by_three_cubes'); model.component.create('comp1',true); model.component('comp1').geom.create('geom1',3); model.component('comp1').geom('geom1').lengthUnit('mm');
tags={'c1','c2','c3','c4','c5','c6'}; names={'C1','C2','C3','C4','C5','C6'}; pos={'0','0','0';'10','0','0';'20','0','0';'0','10','0';'10','10','0';'20','10','0'};
for i=1:6, model.component('comp1').geom('geom1').create(tags{i},'Block'); model.component('comp1').geom('geom1').feature(tags{i}).set('size',{'10','10','10'}); model.component('comp1').geom('geom1').feature(tags{i}).set('pos',pos(i,:)); model.component('comp1').geom('geom1').feature(tags{i}).set('selresult','on'); end
model.component('comp1').geom('geom1').create('uni1','Union'); model.component('comp1').geom('geom1').feature('uni1').selection('input').set(tags); model.component('comp1').geom('geom1').feature('uni1').set('intbnd',true);
btag={'bleft','bright','btop','bbot','bfront','bback','bi12','bi13','bi23'}; box=[-.001 .001 -.001 20.001 -.001 10.001;29.999 30.001 -.001 20.001 -.001 10.001;-.001 30.001 19.999 20.001 -.001 10.001;-.001 30.001 -.001 .001 -.001 10.001;-.001 30.001 -.001 20.001 -.001 .001;-.001 30.001 -.001 20.001 9.999 10.001;9.999 10.001 -.001 10.001 -.001 10.001;9.999 10.001 9.999 20.001 -.001 10.001;9.999 30.001 9.999 10.001 -.001 10.001];
for i=1:numel(btag), f=model.component('comp1').geom('geom1').create(btag{i},'BoxSelection'); f.set('entitydim',2); f.set('xmin',num2str(box(i,1))); f.set('xmax',num2str(box(i,2))); f.set('ymin',num2str(box(i,3))); f.set('ymax',num2str(box(i,4))); f.set('zmin',num2str(box(i,5))); f.set('zmax',num2str(box(i,6))); f.set('condition','inside'); end
model.component('comp1').geom('geom1').run;
did=zeros(1,6); for i=1:6, did(i)=mphgetselection(model.component('comp1').selection(['geom1_' tags{i} '_dom'])).entities(1); end
% DDM partition is stored with the model, rather than duplicated in JSON.
% The numeric suffix is the DDM Subdomain ID consumed by the generic exporter.
ddmDomains={did([1 4]),did([2 3]),did([5 6])};
for s=1:3
  st=sprintf('ddm_sd_%03d',s); model.component('comp1').selection.create(st,'Explicit');
  model.component('comp1').selection(st).label(st); model.component('comp1').selection(st).geom('geom1',3);
  model.component('comp1').selection(st).set(ddmDomains{s});
end
% Geometry BoxSelection nodes identify the faces, but do not themselves create
% component selection tags in 6.1. Query the resulting COMSOL boundary IDs.
bid=cell(1,9); for i=1:9, bid{i}=mphselectbox(model,'geom1',[box(i,1:2);box(i,3:4);box(i,5:6)],'boundary')'; end
k=[10 1 20 2 15 .5]; for i=1:6, mt=['mat' num2str(i)]; model.component('comp1').material.create(mt,'Common'); model.component('comp1').material(mt).selection.set(did(i)); ks=[num2str(k(i)) '[W/(m*K)]']; model.component('comp1').material(mt).propertyGroup('def').set('thermalconductivity',{ks '0' '0' '0' ks '0' '0' '0' ks}); model.component('comp1').material(mt).propertyGroup('def').set('density','1000[kg/m^3]'); model.component('comp1').material(mt).propertyGroup('def').set('heatcapacity','1000[J/(kg*K)]'); end
% COMSOL 6.1 application ID is HeatTransfer; its displayed physics name is
% Heat Transfer in Solids (confirmed from this installation's model files).
ht=model.component('comp1').physics.create('ht','HeatTransfer','geom1'); ht.feature('init1').set('Tinit','300[K]');
ht.create('temp_left','TemperatureBoundary',2); ht.feature('temp_left').selection.set(bid{1}); ht.feature('temp_left').set('T0','300[K]'); ht.create('temp_right','TemperatureBoundary',2); ht.feature('temp_right').selection.set(bid{2}); ht.feature('temp_right').set('T0','350[K]');
ht.create('conv_top','HeatFluxBoundary',2); ht.feature('conv_top').selection.set(bid{3}); ht.feature('conv_top').set('HeatFluxType','ConvectiveHeatFlux'); ht.feature('conv_top').set('h','10[W/(m^2*K)]'); ht.feature('conv_top').set('Text','293.15[K]');
ht.create('flux_bottom','HeatFluxBoundary',2); ht.feature('flux_bottom').selection.set(bid{4}); ht.feature('flux_bottom').set('HeatFluxType','GeneralInwardHeatFlux'); ht.feature('flux_bottom').set('q0','1000[W/m^2]');
ht.create('ins_front','ThermalInsulation',2); ht.feature('ins_front').selection.set(bid{5}); ht.create('ins_back','ThermalInsulation',2); ht.feature('ins_back').selection.set(bid{6}); ht.create('src_c5','HeatSource',3); ht.feature('src_c5').selection.set(did(5)); ht.feature('src_c5').set('Q0','1e6[W/m^3]');
model.component('comp1').mesh.create('mesh_ref'); model.component('comp1').mesh('mesh_ref').feature('size').set('hauto',2); model.component('comp1').mesh('mesh_ref').run;
model.study.create('std_stationary'); model.study('std_stationary').create('stat','Stationary'); model.study('std_stationary').run;
model.study.create('std_time'); model.study('std_time').create('time','Transient'); model.study('std_time').feature('time').set('tlist','range(0,1,100)'); model.study('std_time').run;
model.component('comp1').mesh('mesh_ref').export(fullfile(outputDir,'reference_full.mphtxt'));
sd={did([1 4]),did([2 3]),did([5 6])}; h=[2.5 1.8 2.2];
for s=1:3, mt=['mesh_sd' num2str(s)]; model.component('comp1').mesh.create(mt); % Size must precede FreeTet in the mesh sequence; otherwise COMSOL uses the default size.
model.component('comp1').mesh(mt).create('size1','Size'); model.component('comp1').mesh(mt).feature('size1').selection.geom('geom1',3); model.component('comp1').mesh(mt).feature('size1').selection.set(sd{s}); model.component('comp1').mesh(mt).feature('size1').set('custom','on'); model.component('comp1').mesh(mt).feature('size1').set('hmax',[num2str(h(s)) '[mm]']); model.component('comp1').mesh(mt).feature('size1').set('hmin',[num2str(h(s)/4) '[mm]']);
model.component('comp1').mesh(mt).create('ftet1','FreeTet'); model.component('comp1').mesh(mt).feature('ftet1').selection.geom('geom1',3); model.component('comp1').mesh(mt).feature('ftet1').selection.set(sd{s}); model.component('comp1').mesh(mt).run; model.component('comp1').mesh(mt).export(fullfile(outputDir,['subdomain_' num2str(s) '.mphtxt'])); end
write_map(fullfile(outputDir,'comsol_entity_map.csv'),model,did,bid,names,k,btag); mphsave(model,fullfile(outputDir,'two_by_three_cubes.mph')); disp(did);
end
function write_map(path,model,did,bid,names,k,btag)
fid=fopen(path,'w'); fprintf(fid,'row_type,name,comsol_id,adjacent_domains,subdomain,material_id,k_W_mK,role\n'); sd=[1 2 2 1 3 3];
for i=1:6, fprintf(fid,'domain,%s,%d,%d,%d,mat%d,%.15g,material\n',names{i},did(i),did(i),sd(i),i,k(i)); end
role={'left_temperature','right_temperature','top_convection','bottom_heat_flux','front_insulation','back_insulation','interface_1_2','interface_1_3','interface_2_3'};
for i=1:9, for b=bid{i}, adj=mphgetadj(model,'geom1','domain','boundary',b); adjText=char(strjoin(string(adj(:)'),';')); fprintf(fid,'boundary,%s,%d,"%s",,,,%s\n',btag{i},b,adjText,role{i}); end,end
fclose(fid);
end
