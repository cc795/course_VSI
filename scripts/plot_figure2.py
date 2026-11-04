#!/usr/bin/env python3
"""Plot actual Athena++ snapshots; no synthetic wave patterns."""
import argparse, glob, os, re, sys, json
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/vsi-mpl')
sys.path.insert(0, str(Path(os.environ.get('ATHENA_DIR','/workspace/athena'))/'vis/python'))
import athena_read
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from plot_style import configure_times
configure_times()

def load(directory, number):
    files=sorted(glob.glob(str(Path(directory)/f'vsi.block*.out1.{number:05d}.vtk')))
    if not files: raise ValueError('No matching snapshot files')
    blocks=[]; times=[]
    for f in files:
        with open(f,'rb') as stream:
            header=stream.read(180).decode('ascii',errors='replace')
        times.append(float(re.search(r'time=([\deE+.\-]+)',header)[1]))
        rf,tf,_,d=athena_read.vtk(f)
        blocks.append((rf,tf,d))
    if len(set(times))!=1: raise ValueError('Mixed output times')
    blocks.sort(key=lambda b:b[0][0])
    for a,b in zip(blocks,blocks[1:]):
        if not np.isclose(a[0][-1],b[0][0]): raise ValueError('Missing/overlapping radial blocks')
    rf=np.concatenate([b[0][:-1] for b in blocks]+[blocks[-1][0][-1:]])
    tf=blocks[0][1]
    if any(not np.array_equal(tf,b[1]) for b in blocks): raise ValueError('Unexpected theta decomposition')
    d={key:np.concatenate([b[2][key][0] for b in blocks],axis=1) for key in blocks[0][2]}
    return rf,tf,d,times[0]

def main():
    p=argparse.ArgumentParser();p.add_argument('directory');p.add_argument('--number',type=int,default=None)
    p.add_argument('--output',default='figure2.png');p.add_argument('--require-target',action='store_true');a=p.parse_args()
    if a.number is None:
        files=list(Path(a.directory).glob('vsi.block0.out1.*.vtk'))
        if not files: raise ValueError('No snapshots found')
        a.number=max(int(f.name.split('.')[-2]) for f in files)
    rf,tf,d,time=load(a.directory,a.number)
    if a.require_target and not np.isclose(time,500*2*np.pi,rtol=2e-6):
        raise ValueError(f'Not the Figure 2 target time: {time/(2*np.pi)} orbits')
    r=.75*np.diff(rf**4)/np.diff(rf**3)
    t=np.diff(np.sin(tf)-tf*np.cos(tf))/(-np.diff(np.cos(tf)))
    rr,tt=np.meshgrid(r,t); R=rr*np.sin(tt); z=rr*np.cos(tt)
    # Athena VTK velocity components are spherical orthonormal components.
    v=d['vel']; vr=v[...,0]*np.sin(tt)+v[...,1]*np.cos(tt)
    vz=v[...,0]*np.cos(tt)-v[...,1]*np.sin(tt)
    cs=np.sqrt(.0025/R)
    for key,value in d.items():
        if not np.isfinite(value).all(): raise ValueError(f'Nonfinite {key}')
    if (d['rho']<=0).any() or (d['press']<=0).any(): raise ValueError('Nonpositive rho/P')
    error=np.max(np.abs(d['press']/(d['rho']*cs**2)-1))
    print(json.dumps({'time_P0':time/(2*np.pi),'shape':list(R.shape),'max_temperature_relative_error':float(error),
        'vz_cs_rms_R1_2':float(np.sqrt(np.mean((vz/cs)[(R>=1)&(R<=2)]**2)))},indent=2))
    Rf=rf[None,:]*np.sin(tf[:,None]);zf=rf[None,:]*np.cos(tf[:,None])
    fig=plt.figure(figsize=(9,6))
    grid=fig.add_gridspec(2,2,width_ratios=[3,1],height_ratios=[.055,1],hspace=.15,wspace=.08)
    ax=fig.add_subplot(grid[1,0]);bx=fig.add_subplot(grid[1,1],sharey=ax)
    cax=fig.add_subplot(grid[0,0]);bx.tick_params(labelleft=False)
    im=ax.pcolormesh(Rf,zf,vz/cs,cmap='RdBu_r',norm=SymLogNorm(linthresh=1e-4,vmin=-.3,vmax=.3),rasterized=True)
    ax.set(xlim=(1,2),ylim=(-.6,.6),xlabel=r'$R/R_0$',ylabel=r'$z/R_0$')
    mask=(R>=1.7)&(R<=2)&(np.indices(R.shape)[0]%9==0)&(np.indices(R.shape)[1]%2==0)
    bx.quiver(R[mask],z[mask],vr[mask]/cs[mask],vz[mask]/cs[mask],angles='xy',scale_units='xy',scale=3,width=.008)
    bx.set(xlim=(1.7,2),xlabel=r'$R/R_0$')
    fig.colorbar(im,cax=cax,orientation='horizontal',label=r'$v_z/c_s$')
    cax.xaxis.set_ticks_position('top');cax.xaxis.set_label_position('top')
    fig.suptitle(f'Athena++ simulation, t = {time/(2*np.pi):.3f} P0',y=1.05)
    fig.savefig(a.output,dpi=220,bbox_inches='tight');plt.close(fig)
if __name__=='__main__': main()
