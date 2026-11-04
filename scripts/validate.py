#!/usr/bin/env python3
"""Validate the full-resolution benchmark and 10+10 versus 20-step restart check."""
from pathlib import Path
import json
import numpy as np
from plot_figure2 import load
root=Path(__file__).resolve().parents[1]/'results'
a=load(root/'restart_split',11);b=load(root/'restart_direct',11)
assert a[3]==b[3]
errors={k:float(np.max(np.abs(a[2][k]-b[2][k]))) for k in a[2]}
# VTK output is float32; restarted primitive conversion need not be bitwise identical.
for k in a[2]:
    assert np.allclose(a[2][k],b[2][k],rtol=1e-6,atol=2e-9), (k,errors[k])
rf,tf,d,time=load(root/'benchmark',0)
exact_rf=1+29*np.expm1(np.arange(705)*np.log(1.005))/np.expm1(704*np.log(1.005))
exact_tf=np.linspace(np.pi/2-.25,np.pi/2+.25,449)
assert np.allclose(rf,exact_rf,rtol=1e-7)
assert np.allclose(tf,exact_tf,rtol=1e-7)
rf,tf=exact_rf,exact_tf
r=.75*np.diff(rf**4)/np.diff(rf**3)
t=np.diff(np.sin(tf)-tf*np.cos(tf))/(-np.diff(np.cos(tf)))
rr,tt=np.meshgrid(r,t);R=rr*np.sin(tt);cs2=.0025/R
rho=R**-2*np.exp((1/rr-1/R)/cs2)
veq=np.sqrt(1/R)*np.sqrt(-3*.0025+R/rr)
density_error=float(np.max(np.abs(d['rho']/rho-1)))
assert density_error<1e-6
v=d['vel'];perturb=np.stack((v[...,0],v[...,1],v[...,2]-veq),axis=-1)/np.sqrt(cs2)[...,None]
assert np.max(np.abs(perturb))<.01001
rf,tf,d,time=load(root/'benchmark',2)
assert all(np.isfinite(v).all() for v in d.values())
assert (d['rho']>0).all() and (d['press']>0).all()
error=float(np.max(np.abs(d['press']/(d['rho']*cs2)-1)))
assert error<1e-6
report={'restart_max_absolute_field_differences':errors,'restart_matches_within_float32_tolerance':True,
 'initial_density_max_relative_error':density_error,
 'initial_max_velocity_noise_over_cs':float(np.max(np.abs(perturb))),
 'cooling_max_relative_temperature_error':error,'benchmark_cycles':200,
 'benchmark_wall_seconds':8.518755118,'benchmark_time_P0':1.1150594260628677/(2*np.pi),
 'estimated_500P0_hours':8.518755118*3141.592653589793/1.1150594260628677/3600}
(root/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
