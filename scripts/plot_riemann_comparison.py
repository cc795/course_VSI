#!/usr/bin/env python3
"""Wave-structure schematic and a computed stationary-contact FV comparison."""
from pathlib import Path
import os,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/vsi-mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Polygon
root=Path(__file__).resolve().parents[1]
font_manager.fontManager.addfont('/usr/share/fonts/opentype/noto/NotoSerifCJK-Regular.ttc')
plt.rcParams.update({'font.family':'Noto Serif CJK JP','mathtext.fontset':'stix','axes.unicode_minus':False})
gamma=1.4

def flux(U):
    rho=U[:,0];v=U[:,1]/rho;p=(gamma-1)*(U[:,2]-.5*rho*v*v)
    assert np.min(rho)>0 and np.min(p)>0
    return np.column_stack((rho*v,rho*v*v+p,(U[:,2]+p)*v)),rho,v,p

def riemann(L,R,kind):
    FL,rl,vl,pl=flux(L);FR,rr,vr,pr=flux(R)
    al=np.sqrt(gamma*pl/rl);ar=np.sqrt(gamma*pr/rr)
    SL=np.minimum(vl-al,vr-ar);SR=np.maximum(vl+al,vr+ar)
    if kind=='HLLE':
        middle=(SR[:,None]*FL-SL[:,None]*FR+(SL*SR)[:,None]*(R-L))/(SR-SL)[:,None]
    else:
        SM=(pr-pl+rl*vl*(SL-vl)-rr*vr*(SR-vr))/(rl*(SL-vl)-rr*(SR-vr))
        ps=pl+rl*(SL-vl)*(SM-vl)
        def star(U,r,v,p,S):
            rs=r*(S-v)/(S-SM)
            es=((S-v)*U[:,2]-p*v+ps*SM)/(S-SM)
            return np.column_stack((rs,rs*SM,es))
        FLs=FL+SL[:,None]*(star(L,rl,vl,pl,SL)-L)
        FRs=FR+SR[:,None]*(star(R,rr,vr,pr,SR)-R)
        middle=np.where((SM>=0)[:,None],FLs,FRs)
    return np.where((SL>=0)[:,None],FL,np.where((SR<=0)[:,None],FR,middle))

def integrate(kind):
    N=200;dx=1/N;x=(np.arange(N)+.5)*dx
    rho=np.where(x<.5,1.,.125);U=np.column_stack((rho,np.zeros(N),np.full(N,1/(gamma-1))))
    time=0.;steps=0
    while time<.05:
        _,r,v,p=flux(U);dt=min(.4*dx/np.max(np.abs(v)+np.sqrt(gamma*p/r)),.05-time)
        ext=np.vstack((U[0],U,U[-1]));F=riemann(ext[:-1],ext[1:],kind)
        U-=dt/dx*(F[1:]-F[:-1]);time+=dt;steps+=1
    _,r,v,p=flux(U)
    mixed=(r>.125+.1*.875)&(r<1-.1*.875)
    return x,r,{'steps':steps,'L1_density_error':float(np.mean(abs(r-rho))),'transition_width_10_90':float(np.sum(mixed)*dx),'max_pressure_error':float(np.max(abs(p-1))),'max_velocity':float(np.max(abs(v)))}

def fan(ax,kind,color):
    left,right=-1.05,1.05;sm=.18
    if kind=='HLLE':
        ax.add_patch(Polygon([[0,0],[left,1],[right,1]],facecolor=color,alpha=.17))
        ax.text(0,.63,r'$U_{\rm HLL}$',ha='center',fontsize=20,color=color)
    else:
        ax.add_patch(Polygon([[0,0],[left,1],[sm,1]],facecolor=color,alpha=.14))
        ax.add_patch(Polygon([[0,0],[sm,1],[right,1]],facecolor=color,alpha=.27))
        ax.plot([0,sm],[0,1],color=color,lw=2.7)
        ax.text(sm+.07,1.03,r'$S_M$',fontsize=14,color=color)
        ax.text(-.38,.63,r'$U_L^*$',ha='center',fontsize=20,color=color)
        ax.text(.60,.63,r'$U_R^*$',ha='center',fontsize=20,color=color)
        ax.annotate('接触波',xy=(sm*.40,.4),xytext=(.60,.22),fontsize=12,color=color,
                    arrowprops={'arrowstyle':'->','color':color})
    for speed,name in [(left,'S_L'),(right,'S_R')]:
        ax.plot([0,speed],[0,1],color='#35445B',lw=2)
        ax.text(speed,1.035,'$'+name+'$',ha='center',fontsize=14,color='#35445B')
    ax.text(-1.15,.55,r'$U_L$',fontsize=17,ha='center');ax.text(1.18,.55,r'$U_R$',fontsize=17,ha='center')
    ax.annotate('',xy=(1.45,0),xytext=(-1.4,0),arrowprops={'arrowstyle':'->','color':'#63738A'})
    ax.annotate('',xy=(0,1.16),xytext=(0,0),arrowprops={'arrowstyle':'->','color':'#63738A'})
    ax.text(1.42,-.1,'x',fontsize=13);ax.text(.04,1.12,'t',fontsize=13)
    ax.set(xlim=(-1.5,1.5),ylim=(-.16,1.2));ax.axis('off')
    ax.set_title(kind+('：两个外侧波速界限' if kind=='HLLE' else '：额外恢复中间接触波'),fontsize=15,fontweight='bold',color=color,pad=13)

def main():
    x,hll,eh=integrate('HLLE');_,hllc,ec=integrate('HLLC')
    assert eh['L1_density_error']>0 and ec['L1_density_error']<1e-12
    assert max(eh['max_pressure_error'],ec['max_pressure_error'])<1e-10
    red='#C94038';teal='#008B9B'
    fig=plt.figure(figsize=(12,9),dpi=180)
    gs=fig.add_gridspec(2,2,height_ratios=[1,1.08],left=.09,right=.96,bottom=.16,top=.86,hspace=.50,wspace=.18)
    fan(fig.add_subplot(gs[0,0]),'HLLE',red);fan(fig.add_subplot(gs[0,1]),'HLLC',teal)
    ax=fig.add_subplot(gs[1,:]);exact=np.where(x<.5,1.,.125)
    ax.step(x,exact,where='mid',color='#172235',lw=4,alpha=.22,label='精确解')
    ax.plot(x,hll,color=red,lw=2.7,label='HLLE')
    ax.step(x,hllc,where='mid',color=teal,lw=2,label='HLLC')
    ax.plot(x,hllc,'o',color=teal,ms=3,mfc='white',markevery=2)
    ax.set(xlim=(.36,.64),ylim=(.05,1.1),xlabel='位置 x',ylabel=r'密度 $\rho$')
    ax.set_title('数值小实验：相同压强和速度，仅密度发生跳变',fontsize=15,pad=13)
    ax.grid(alpha=.15);ax.legend(loc='upper right',frameon=False,ncol=3,fontsize=11)
    ax.annotate('HLLE：密度跳变被展宽',xy=(.48,np.interp(.48,x,hll)),xytext=(.385,.58),color=red,fontsize=12,
                arrowprops={'arrowstyle':'->','color':red})
    ax.annotate('HLLC：保留静止接触间断',xy=(.5,.60),xytext=(.535,.38),color=teal,fontsize=12,
                arrowprops={'arrowstyle':'->','color':teal})
    fig.suptitle('HLLE 与 HLLC：为什么接触间断的清晰程度不同？',fontsize=21,y=.97)
    fig.text(.525,.915,'上：近似黎曼波结构示意。下：实际计算的一维静止接触间断。',ha='center',fontsize=12,color='#4B5B70')
    fig.text(.525,.095,r'初态：$\rho_L=1,\ \rho_R=0.125,\ v_L=v_R=0,\ P_L=P_R=1$；$\gamma=1.4$。',ha='center',fontsize=11)
    fig.text(.525,.055,'两者采用同一网格、CFL 和外侧波速估计；一阶有限体积，200 格，t=0.05。',ha='center',fontsize=10,color='#4B5B70')
    fig.text(.525,.022,'本例的 HLLC 可保持静止接触；移动接触或高阶实际模拟仍可能存在数值扩散。',ha='center',fontsize=10,color='#4B5B70')
    for ext in ['png','pdf']:fig.savefig(root/'media'/f'hlle_hllc.{ext}')
    plt.close(fig)
    report={'example':'stationary_contact_1D_Euler','grid_cells':200,'gamma':gamma,'final_time':.05,'CFL':.4,'reconstruction':'piecewise constant','wave_bounds':'common Davis left/right bounds','HLLE':eh,'HLLC':ec}
    (root/'media/hlle_hllc.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
