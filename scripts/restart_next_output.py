#!/usr/bin/env python3
"""Schedule anchor for a pinned, double-precision Athena++ restart."""
import ctypes, math, struct, sys
class RegionSize(ctypes.Structure):
    _fields_=[('bounds_and_ratios',ctypes.c_double*9),('cells',ctypes.c_int*3)]
with open(sys.argv[1],'rb') as f:
    header=f.read(40960)
    marker=b'<par_end>\n'
    end=header.index(marker)+len(marker)
    # restart.cpp: int nbtotal, int root_level, RegionSize, double time.
    f.seek(end+2*ctypes.sizeof(ctypes.c_int)+ctypes.sizeof(RegionSize))
    time=struct.unpack('@d',f.read(8))[0]
if not math.isfinite(time) or not 0<=time<=500*2*math.pi+1e-6:
    raise ValueError('Unexpected restart time/format')
# Athena++ ForwardNextTime advances this anchor to the next integer orbit.
print(format(math.floor(time/(2*math.pi))*2*math.pi,'.17g'))
