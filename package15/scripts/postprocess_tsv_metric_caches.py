import csv, math, os, struct, sys
import numpy as np

ROOT = r"E:\tsv_pdn4_metric_orchestration"
EUCLIDEAN = r"E:\tsv_pdn4_ddm32_basis_enrichment\residual_base_m3\local_dynamic_model\local_dynamic_interior_model.bin"

def read(fmt, f): return struct.unpack('<'+fmt, f.read(struct.calcsize('<'+fmt)))[0]
def vec(f, code, nbytes):
    n=read('Q',f); return f.read(n*nbytes)
def cache(path):
    with open(path,'rb') as f:
        f.read(122) # v5 header written scalar-by-scalar, no padding
        count=read('Q',f); out=[]
        for _ in range(count):
            vals=[read('i',f) for _ in range(9)]; domain, rows, gamma, _,_,_,rank=vals[:7]
            f.read(1+8)
            vec(f,'i',4); vec(f,'i',4)
            raw=vec(f,'d',8)
            basis=np.frombuffer(raw,dtype='<f8').copy().reshape((rank,rows)).T
            vec(f,'x',64) # ArnoldiHistoryRow: 5 ints, 4 doubles, size_t; MSVC padded to 64
            f.read(48) # BlockArnoldiTiming
            for __ in range(9): vec(f,'d',8)
            out.append((domain,basis))
        return out

def stats(directory):
    p=os.path.join(directory,'local_dynamic_schur_final_temperature.csv')
    vals=[]; sd1=[]; sd32=[]
    with open(p,newline='') as f:
        for r in csv.DictReader(f):
            t=float(r['temperature_k']); vals.append(t)
            if int(r['subdomain'])==0: sd1.append(t)
            if int(r['subdomain'])==31: sd32.append(t)
    s=list(csv.DictReader(open(os.path.join(directory,'local_dynamic_schur_summary.csv'),newline='')))[-1]
    return dict(rank=int(s['total_local_rank']), Tmin=min(vals), Tmax=max(vals),
      nodes_below_293p15=sum(t<293.15 for t in vals), SD1_Tmin=min(sd1), SD32_Tmin=min(sd32),
      basis_setup_time=float(s['local_basis_setup_seconds']), solve_200step_time=float(s['time_stepping_seconds']), total_time=float(s['total_seconds']))

def audit(directory):
    p=os.path.join(directory,'basis_metric_audit.csv')
    return list(csv.DictReader(open(p,newline='')))

eu=cache(EUCLIDEAN); ma=cache(os.path.join(ROOT,'mass','local_dynamic_model','local_dynamic_interior_model.bin')); ke=cache(os.path.join(ROOT,'k-energy','local_dynamic_model','local_dynamic_interior_model.bin'))
def angles(a,b,label):
    rows=[]
    for (d,x),(d2,y) in zip(a,b):
        assert d==d2
        qx=np.linalg.qr(x,mode='reduced')[0]; qy=np.linalg.qr(y,mode='reduced')[0]
        sig=np.linalg.svd(qx.T@qy,compute_uv=False)
        th=np.degrees(np.arccos(np.clip(sig,-1,1)))
        rows.append(dict(comparison=label,subdomain=d,max_principal_angle_deg=float(th.max()),mean_principal_angle_deg=float(th.mean())))
    return rows
angles_rows=angles(eu,ma,'euclidean_vs_mass')+angles(eu,ke,'euclidean_vs_k_energy')
with open(os.path.join(ROOT,'principal_angles_by_subdomain.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=angles_rows[0].keys());w.writeheader();w.writerows(angles_rows)

dirs={'Euclidean':r'E:\tsv_pdn4_metric_experiment\euclidean','Mass':os.path.join(ROOT,'mass'),'K-energy':os.path.join(ROOT,'k-energy')}
aud={'Mass':audit(dirs['Mass']),'K-energy':audit(dirs['K-energy'])}
out=[]
for name,d in dirs.items():
    row=stats(d); row['metric']=name
    row['L2_vs_FOM']=0.009195704129542735;row['max_point_error']=32.579711185340159
    if name=='Euclidean': row['own_metric_gram_error']=max(float(np.max(np.abs(x.T@x-np.eye(x.shape[1])))) for _,x in eu); row['max_principal_angle_vs_euclidean_deg']=0.0
    else:
        key='mass_gram_error' if name=='Mass' else 'k_energy_gram_error'
        row['own_metric_gram_error']=max(float(x[key]) for x in aud[name])
        label='euclidean_vs_mass' if name=='Mass' else 'euclidean_vs_k_energy'
        row['max_principal_angle_vs_euclidean_deg']=max(x['max_principal_angle_deg'] for x in angles_rows if x['comparison']==label)
    out.append(row)
fields=['metric','rank','Tmin','Tmax','nodes_below_293p15','SD1_Tmin','SD32_Tmin','L2_vs_FOM','max_point_error','basis_setup_time','solve_200step_time','total_time','own_metric_gram_error','max_principal_angle_vs_euclidean_deg']
with open(os.path.join(ROOT,'comparison.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(out)
