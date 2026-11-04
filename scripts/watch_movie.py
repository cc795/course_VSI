#!/usr/bin/env python3
"""Update the movie when the next complete orbit snapshot becomes available."""
from pathlib import Path
import subprocess,time,json
root=Path(__file__).resolve().parents[1]
from make_movie import discover
previous=-1
while True:
    last=len(discover())-1
    if last>=11 and last!=previous:
        result=subprocess.run(['python3',str(root/'scripts/make_movie.py'),'--min-last','11'])
        if result.returncode==0: previous=last
    if previous>=500: break
    time.sleep(30)
