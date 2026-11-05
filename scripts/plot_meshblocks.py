#!/usr/bin/env python3
"""Plot the actual 2D spherical-polar MeshBlock decomposition from the input."""
from pathlib import Path
import os,re,argparse,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/vsi-mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,Rectangle,Patch
from plot_style import configure_times
configure_times()
root=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser();p.add_argument('--input',type=Path,default=root/'inputs/athinput.figure2');a=p.parse_args()
    values={};section=None
    for line in a.input.read_text().splitlines():
        line=line.split('#',1)[0].strip()
        if re.fullmatch(r'<[^>]+>',line):section=line[1:-1];values.setdefault(section,{})
        elif '=' in line:
            k,v=line.split('=',1);values[section][k.strip()]=v.strip()
    m=values['mesh'];b=values['meshblock']
    nr,nt=int(m['nx1']),int(m['nx2']);br,bt=int(b['nx1']),int(b['nx2'])
    assert nr%br==0 and nt==bt and int(m['nx3'])==1
    rlo,rhi=float(m['x1min']),float(m['x1max']);q=float(m.get('x1rat',1))
    r=np.linspace(rlo,rhi,nr+1) if q==1 else rlo+(rhi-rlo)*np.expm1(np.arange(nr+1)*np.log(q))/np.expm1(nr*np.log(q))
    t=np.linspace(float(m['x2min']),float(m['x2max']),nt+1)
    edges=np.arange(0,nr+1,br);bounds=r[edges]
    colors=['#26C6CB','#478DDD','#A07CE0','#ED9B65'];bg='#0B1220';ink='#E8F0FC';muted='#ADBED5'
    fig=plt.figure(figsize=(13,9),dpi=180,facecolor=bg)
    gs=fig.add_gridspec(2,1,height_ratios=[3.1,1],left=.085,right=.97,bottom=.09,top=.84,hspace=.5)
    ax=fig.add_subplot(gs[0],facecolor=bg);bx=fig.add_subplot(gs[1],facecolor=bg)
    for n,(lo,hi) in enumerate(zip(bounds[:-1],bounds[1:])):
        poly=np.column_stack([np.r_[lo*np.sin(t),hi*np.sin(t[::-1])],np.r_[lo*np.cos(t),hi*np.cos(t[::-1])]])
        ax.add_patch(Polygon(poly,facecolor=colors[n%4],edgecolor=colors[n%4],alpha=.25,lw=1.2))
        mid=.5*(lo+hi)
        if n<2:
            ax.annotate(f'MeshBlock {n}',xy=(mid,mid*.255),xytext=([2.7,7.5][n],[3.7,4.9][n]),
                        color=ink,ha='center',fontsize=12,fontweight='bold',
                        arrowprops={'arrowstyle':'-','color':ink,'lw':.9})
        else:
            ax.text(mid,0,f'MeshBlock {n}',color=ink,ha='center',va='center',fontsize=13,fontweight='bold')
        bx.add_patch(Rectangle((n*br,0),br,nt,facecolor=colors[n%4],edgecolor=ink,lw=1.3,alpha=.65))
        bx.text((n+.5)*br,nt/2,f'MeshBlock {n}\n{br} × {bt} cells',color=ink,ha='center',va='center',fontsize=13)
    for rr in r[::8]:ax.plot(rr*np.sin(t),rr*np.cos(t),color=muted,lw=.3,alpha=.28)
    for th in t[::16]:ax.plot(r*np.sin(th),r*np.cos(th),color=muted,lw=.3,alpha=.28)
    for rr in bounds:ax.plot(rr*np.sin(t),rr*np.cos(t),color=ink,lw=1.2,alpha=.85)
    for th in [t[0],t[-1]]:ax.plot(r*np.sin(th),r*np.cos(th),color=ink,lw=1.2,alpha=.85)
    ax.set(xlim=(0,31),ylim=(-8.5,8.5),xlabel=r'$R/R_0$',ylabel=r'$z/R_0$')
    ax.set_title('Physical mesh: spherical-polar cells in the R–z plane',color=ink,fontsize=15,pad=13)
    bx.set(xlim=(0,nr),ylim=(nt,0),xticks=edges,yticks=[0,nt//2,nt],xlabel='Radial cell index i',ylabel='Polar cell index j')
    bx.set_title('Logical mesh: all four blocks have the same number of cells',color=ink,fontsize=14,pad=10)
    for panel in [ax,bx]:
        panel.tick_params(colors=muted,labelsize=11)
        for spine in panel.spines.values():spine.set_color('#43536D')
        panel.xaxis.label.set_color(ink);panel.yaxis.label.set_color(ink)
        panel.xaxis.label.set_size(13);panel.yaxis.label.set_size(13)
    handles=[Patch(facecolor=c,label=f'MeshBlock {n}: r = {bounds[n]:.3f}–{bounds[n+1]:.3f}') for n,c in enumerate(colors)]
    legend=fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.53,.91),ncol=4,frameon=False,fontsize=10)
    for text in legend.get_texts():text.set_color(ink)
    fig.suptitle('MeshBlock layout',color=ink,fontsize=24,y=.985)
    fig.text(.53,.93,f'{nr} × {nt} active cells   |   {len(bounds)-1} MeshBlocks   |   adjacent radial cell-width ratio {q}',ha='center',color=muted,fontsize=12)
    fig.text(.53,.025,'Cell lines in the physical view are sampled every 8 radial and 16 polar cells for readability.',ha='center',color=muted,fontsize=10)
    for ext in ['png','svg','pdf']:
        path=root/'media'/f'meshblocks.{ext}'
        fig.savefig(path,facecolor=bg)
        if ext=='svg':path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)
    (root/'media/meshblocks.json').write_text(json.dumps({'mesh':[nr,nt,1],'meshblock':[br,bt,1],'radial_block_boundaries':bounds.tolist(),'source_input':'inputs/athinput.figure2'},indent=2)+'\n')
    print('Radial block boundaries:',bounds.tolist())
if __name__=='__main__':main()
