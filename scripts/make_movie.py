#!/usr/bin/env python3
"""Render all contiguous integer-orbit snapshots and concatenate cached H.264 clips."""
from pathlib import Path
import argparse, json, re, subprocess, os, fcntl
import numpy as np
from plot_figure2 import load, plt, SymLogNorm
root=Path(__file__).resolve().parents[1]

def discover():
    selected={}
    for directory in [root/'results/movie_early',root/'results/figure2']:
        for file in sorted(directory.glob('vsi.block0.out1.*.vtk')):
            number=int(file.name.split('.')[-2])
            files=list(directory.glob(f'vsi.block*.out1.{number:05d}.vtk'))
            if len(files)!=4: continue
            try:
                with file.open('rb') as f: header=f.read(180).decode(errors='replace')
                time=float(re.search(r'time=([\deE+.\-]+)',header)[1])
                orbit=round(time/(2*np.pi))
                if abs(time/(2*np.pi)-orbit)>.002: continue
                # Ignore files while Athena is still writing them.
                if any(f.stat().st_size<1500000 for f in files): continue
                selected[orbit]=(directory,number,time)
            except (ValueError,OSError): continue
    contiguous=[]
    while len(contiguous) in selected: contiguous.append(selected[len(contiguous)])
    return contiguous

def render(directory,number,time,n,png):
    rf,tf,d,tactual=load(directory,number)
    assert np.isclose(time,tactual)
    assert len(rf)==705 and len(tf)==449
    assert all(np.isfinite(a).all() for a in d.values())
    assert (d['rho']>0).all() and (d['press']>0).all()
    r=.75*np.diff(rf**4)/np.diff(rf**3)
    t=np.diff(np.sin(tf)-tf*np.cos(tf))/(-np.diff(np.cos(tf)))
    rr,tt=np.meshgrid(r,t);R=rr*np.sin(tt);z=rr*np.cos(tt);cs=np.sqrt(.0025/R)
    v=d['vel'];u=(v[...,0]*np.sin(tt)+v[...,1]*np.cos(tt))/cs
    w=(v[...,0]*np.cos(tt)-v[...,1]*np.sin(tt))/cs
    mask=(R>=1.7)&(R<=2)&(np.indices(R.shape)[0]%10==0)&(np.indices(R.shape)[1]%2==0)
    fig=plt.figure(figsize=(12,8),dpi=120)
    gs=fig.add_gridspec(2,2,width_ratios=[3,1],height_ratios=[.055,1],left=.08,right=.96,bottom=.10,top=.82,hspace=.18,wspace=.08)
    ax=fig.add_subplot(gs[1,0]);bx=fig.add_subplot(gs[1,1],sharey=ax);cax=fig.add_subplot(gs[0,0])
    im=ax.pcolormesh(rf[None,:]*np.sin(tf[:,None]),rf[None,:]*np.cos(tf[:,None]),w,
       cmap='RdBu_r',norm=SymLogNorm(linthresh=1e-4,vmin=-.3,vmax=.3),rasterized=True)
    ax.set(xlim=(1,2),ylim=(-.6,.6),xlabel=r'$R/R_0$',ylabel=r'$z/R_0$')
    bx.quiver(R[mask],z[mask],u[mask],w[mask],angles='xy',scale_units='xy',scale=3,width=.006,minlength=.4)
    bx.set(xlim=(1.7,2),xlabel=r'$R/R_0$');bx.tick_params(labelleft=False)
    cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-.1,-.001,0,.001,.1]);cb.set_label(r'$v_z/c_s$')
    cax.xaxis.set_ticks_position('top');cax.xaxis.set_label_position('top')
    fig.suptitle(rf'Athena++ VSI simulation: $t\simeq{n}P_0$',y=.965,fontsize=20)
    fig.text(.5,.90,'Svanberg et al. (2022), Section 3 | 704 x 448',ha='center',fontsize=11)
    fig.text(.5,.025,f'Actual t = {time/(2*np.pi):.6f} P0 | Fixed colour/arrow scales | No temporal interpolation',ha='center',fontsize=10)
    fig.savefig(png);plt.close(fig)
    return {'orbit':n,'time_P0':time/(2*np.pi),'directory':str(directory),'snapshot_number':number,
            'vz_cs_rms_R1_2':float(np.sqrt(np.mean(w[(R>=1)&(R<=2)]**2)))}

def main():
    p=argparse.ArgumentParser();p.add_argument('--min-last',type=int,default=0);a=p.parse_args()
    lock_dir=root/'results/movie_figure2';lock_dir.mkdir(exist_ok=True)
    lock=(lock_dir/'.movie.lock').open('w')
    fcntl.flock(lock,fcntl.LOCK_EX)
    data=discover()
    if not data or len(data)-1<a.min_last: raise ValueError(f'Only {len(data)} contiguous snapshots available')
    out=root/'results/movie_figure2';out.mkdir(exist_ok=True)
    frames=out/'frames';frames.mkdir(exist_ok=True)
    clips=out/'clips';clips.mkdir(exist_ok=True)
    metadata=[]
    for n,(directory,number,time) in enumerate(data):
        png=frames/f'frame_{n:04d}.png';clip=clips/f'clip_{n:04d}.mp4';meta=frames/f'frame_{n:04d}.json'
        if not clip.exists():
            info=render(directory,number,time,n,png)
            subprocess.run(['ffmpeg','-v','error','-y','-loop','1','-framerate','30','-i',str(png),'-t','1','-c:v','libx264','-threads','1','-crf','18','-preset','fast','-pix_fmt','yuv420p',str(clip)],check=True)
            meta.write_text(json.dumps(info,indent=2)+'\n')
        metadata.append(json.loads(meta.read_text()))
    concat=out/'clips.txt';concat.write_text(''.join(f"file '{clips/f'clip_{n:04d}.mp4'}'\n" for n in range(len(data))))
    tmp=out/'figure2_evolution.tmp.mp4';target=out/'figure2_evolution.mp4'
    subprocess.run(['ffmpeg','-v','error','-y','-f','concat','-safe','0','-i',str(concat),'-c','copy','-movflags','+faststart',str(tmp)],check=True)
    os.replace(tmp,target)
    (out/'movie_metadata.json').write_text(json.dumps({'seconds_per_snapshot':1,'frames':metadata,'last_orbit':len(data)-1},indent=2)+'\n')
    print(f'{target}: orbits 0..{len(data)-1}, duration {len(data)} seconds',flush=True)
if __name__=='__main__': main()
