% Mesh-B: read-only source MPH -> independent local-refined DDM MPHTXT.
% No mphsave/model.save is called by this workflow.
addpath('D:\AI agent\codex\DDM_Schur\scripts');
addpath('D:\Program Files\COMSOL\COMSOL61\Multiphysics\mli');
mphstart('localhost',20361);
sourceMph='E:\TSV PDN4_all_clear_ddm32.mph';
config='D:\AI agent\codex\DDM_Schur\scripts\tsv_pdn4_ddm32_mesh_b_partition.json';
output='E:\tsv_pdn4_mesh_b_local_refined';
export_comsol_ddm_case(sourceMph,config,output);
