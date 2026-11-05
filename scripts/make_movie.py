#!/usr/bin/env python3
"""Render all contiguous integer-orbit snapshots and concatenate cached H.264 clips."""
from pathlib import Path
import argparse, json, re, subprocess, os, fcntl
import numpy as np
from plot_figure2 import load, plt, SymLogNorm
from matplotlib.colors import LinearSegmentedColormap
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
    v=d['vel']
    w=(v[...,0]*np.cos(tt)-v[...,1]*np.sin(tt))/cs
    background='#080d1c'
    cmap=LinearSegmentedColormap.from_list('neon_diverging',
       ['#00f0ff','#007acc','#193b74',background,'#57151c','#c91d25','#ff0000'])
    fig=plt.figure(figsize=(12,8),dpi=120,facecolor=background)
    gs=fig.add_gridspec(2,1,height_ratios=[.045,1],left=.10,right=.95,bottom=.11,top=.82,hspace=.22)
    ax=fig.add_subplot(gs[1,0],facecolor=background)
    cax=fig.add_subplot(gs[0,0],facecolor=background)
    im=ax.pcolormesh(rf[None,:]*np.sin(tf[:,None]),rf[None,:]*np.cos(tf[:,None]),w,
       cmap=cmap,norm=SymLogNorm(linthresh=1e-4,vmin=-.3,vmax=.3),rasterized=True)
    ax.set(xlim=(1,2),ylim=(-.6,.6),xlabel=r'$R/R_0$',ylabel=r'$z/R_0$')
    cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-.1,-.001,0,.001,.1])
    cb.set_label(r'$v_z/c_s$',color='#e9f2ff',fontsize=15)
    cax.xaxis.set_ticks_position('top');cax.xaxis.set_label_position('top')
    for panel in (ax,cax):
        panel.tick_params(colors='#c6d7ed',labelsize=12)
        panel.xaxis.label.set_color('#e9f2ff');panel.yaxis.label.set_color('#e9f2ff')
        panel.xaxis.label.set_size(16);panel.yaxis.label.set_size(16)
        for spine in panel.spines.values(): spine.set_color('#42516b')
    cb.outline.set_edgecolor('#42516b')
    fig.suptitle(rf'$t={n}P_0$',y=.955,fontsize=28,color='#e9f2ff')
    fig.savefig(png,facecolor=background);plt.close(fig)
    return {'orbit':n,'time_P0':time/(2*np.pi),'directory':str(directory),'snapshot_number':number,
            'vz_cs_rms_R1_2':float(np.sqrt(np.mean(w[(R>=1)&(R<=2)]**2)))}

def main():
    p=argparse.ArgumentParser();p.add_argument('--min-last',type=int,default=0)
    p.add_argument('--duration',type=float,default=40.0,help='Total movie duration in seconds')
    a=p.parse_args()
    if not np.isfinite(a.duration) or a.duration<=0: p.error('--duration must be positive and finite')
    lock_dir=root/'results/movie_figure2';lock_dir.mkdir(exist_ok=True)
    lock=(lock_dir/'.movie.lock').open('w')
    fcntl.flock(lock,fcntl.LOCK_EX)
    data=discover()
    if not data or len(data)-1<a.min_last: raise ValueError(f'Only {len(data)} contiguous snapshots available')
    out=root/'results/movie_figure2';out.mkdir(exist_ok=True)
    frames=out/'frames_neon_times_red_v4';frames.mkdir(exist_ok=True)
    clips=out/'clips_neon_times_red_v4_double_speed';clips.mkdir(exist_ok=True)
    metadata=[]
    for n,(directory,number,time) in enumerate(data):
        png=frames/f'frame_{n:04d}.png';clip=clips/f'clip_{n:04d}.mp4';meta=frames/f'frame_{n:04d}.json'
        if not clip.exists():
            if png.exists() and meta.exists():
                info=json.loads(meta.read_text())
            else:
                info=render(directory,number,time,n,png)
            tmp_clip=clip.with_suffix('.tmp.mp4')
            subprocess.run(['ffmpeg','-v','error','-y','-loop','1','-framerate','24','-i',str(png),'-frames:v','3','-c:v','libx264','-threads','1','-crf','18','-preset','fast','-pix_fmt','yuv420p',str(tmp_clip)],check=True)
            os.replace(tmp_clip,clip)
            meta.write_text(json.dumps(info,indent=2)+'\n')
        metadata.append(json.loads(meta.read_text()))
    concat=out/'clips.txt';concat.write_text(''.join(f"file '{clips/f'clip_{n:04d}.mp4'}'\n" for n in range(len(data))))
    tmp=out/'figure2_evolution.tmp.mp4';target=out/'figure2_evolution.mp4'
    subprocess.run(['ffmpeg','-v','error','-xerror','-y','-f','concat','-safe','0','-i',str(concat),'-c','copy','-movflags','+faststart',str(tmp)],check=True)
    retimed=out/'figure2_evolution.retimed.tmp.mp4'
    speed=a.duration/(len(data)*.125)
    filters=f'setpts={speed:.17g}*(PTS-STARTPTS),tpad=stop_mode=clone:stop_duration=1,fps=30'
    subprocess.run(['ffmpeg','-v','error','-xerror','-y','-i',str(tmp),'-vf',filters,
                    '-t',str(a.duration),'-an','-c:v','libx264','-threads','2',
                    '-crf','18','-preset','fast','-pix_fmt','yuv420p',
                    '-movflags','+faststart',str(retimed)],check=True)
    os.replace(retimed,target)
    tmp.unlink()
    (out/'movie_metadata.json').write_text(json.dumps({'seconds_per_snapshot':a.duration/len(data),'duration_seconds':a.duration,'video_fps':30,'timing_note':'seconds_per_snapshot is the mean display duration; actual holds are quantized to video frames.','style':'neon_single_panel_times_red_v4','frames':metadata,'last_orbit':len(data)-1},indent=2)+'\n')
    print(f'{target}: orbits 0..{len(data)-1}, duration {a.duration:g} seconds',flush=True)
if __name__=='__main__': main()
