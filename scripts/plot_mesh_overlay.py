#!/usr/bin/env python3
"""Overlay sampled cell edges and actual MeshBlock boundaries on a real snapshot."""
from pathlib import Path
import numpy as np
from plot_figure2 import load,plt,SymLogNorm
root=Path(__file__).resolve().parents[1]

def main():
    directory=root/'results/figure2'
    numbers=[int(p.name.split('.')[-2]) for p in directory.glob('vsi.block0.out1.*.vtk')]
    rf,tf,d,time=load(directory,max(numbers))
    r=.75*np.diff(rf**4)/np.diff(rf**3);t=np.diff(np.sin(tf)-tf*np.cos(tf))/(-np.diff(np.cos(tf)))
    rr,tt=np.meshgrid(r,t);R=rr*np.sin(tt);cs=np.sqrt(.0025/R)
    vz=d['vel'][...,0]*np.cos(tt)-d['vel'][...,1]*np.sin(tt)
    fig,ax=plt.subplots(figsize=(13,7),dpi=180)
    fig.subplots_adjust(left=.08,right=.96,bottom=.15,top=.78)
    im=ax.pcolormesh(rf[None,:]*np.sin(tf[:,None]),rf[None,:]*np.cos(tf[:,None]),vz/cs,
                    cmap='RdBu_r',norm=SymLogNorm(linthresh=1e-4,vmin=-.3,vmax=.3),rasterized=True)
    for rface in rf[::8]:ax.plot(rface*np.sin(tf),rface*np.cos(tf),color='black',lw=.3,alpha=.65)
    for theta in tf[::16]:ax.plot(rf*np.sin(theta),rf*np.cos(theta),color='black',lw=.3,alpha=.65)
    boundaries=rf[[0,176,352,528,704]]
    for rface in boundaries:ax.plot(rface*np.sin(tf),rface*np.cos(tf),color='black',lw=1.4)
    for theta in [tf[0],tf[-1]]:ax.plot(rf*np.sin(theta),rf*np.cos(theta),color='black',lw=1.4)
    for n,(lo,hi) in enumerate(zip(boundaries[:-1],boundaries[1:])):
        mid=.5*(lo+hi)
        ax.annotate(f'MeshBlock {n}',xy=(mid,mid*.255),
                    xytext=([2.4,6,11,22][n],[2.8,3.8,5.4,8.1][n]),
                    ha='center',fontsize=12,fontweight='bold',
                    arrowprops={'arrowstyle':'-','color':'black','lw':.8})
    ax.set(xlim=(.95,30),ylim=(-8.5,9.4),xlabel=r'$R/R_0$',ylabel=r'$z/R_0$')
    fig.suptitle(rf'MeshBlocks and flow at $t={time/(2*np.pi):.0f}P_0$',fontsize=22,y=.98)
    colorax=fig.add_axes([.17,.855,.70,.025])
    cb=fig.colorbar(im,cax=colorax,orientation='horizontal',ticks=[-.1,-.001,0,.001,.1]);cb.set_label(r'$v_z/c_s$')
    colorax.xaxis.set_ticks_position('top');colorax.xaxis.set_label_position('top')
    fig.text(.52,.055,'Thick lines: 4 MeshBlock boundaries; thin grid lines: every 8 radial and 16 polar cells.',ha='center',fontsize=11)
    fig.text(.52,.02,'Full mesh: 704 × 448 cells. Each block: 176 × 448 cells. No mesh refinement.',ha='center',fontsize=11)
    for ext in ['png','pdf']:fig.savefig(root/'media'/f'mesh_overlay.{ext}')
    plt.close(fig)
if __name__=='__main__':main()
