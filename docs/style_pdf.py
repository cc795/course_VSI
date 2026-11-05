#!/usr/bin/env python3
"""Restore the user's blue code selections after PDF reflow."""
from pathlib import Path
from collections import defaultdict
import json,re
import fitz
root=Path(__file__).resolve().parents[1]
source=(root/'src/svanberg.cpp').read_text().splitlines()
settings=json.loads((root/'docs/pdf_highlights.json').read_text())
pdf=fitz.open(root/'readme.pdf');by_line=defaultdict(list);gutters={}
current=None;cpp=True
for pn,page in enumerate(pdf):
    events=[]
    for block in page.get_text('rawdict')['blocks']:
        for line in block.get('lines',[]):
            for span in line['spans']:
                if span['color'] not in (0x20344B,0x6A8297):continue
                text=''.join(c['c'] for c in span['chars'])
                events.append((round(span['origin'][1],1),span['origin'][0],span,text))
    for _,_,span,text in sorted(events,key=lambda e:(e[0],e[1])):
        if span['color']==0x6A8297 and text.strip().isdigit():
            number=int(text.strip())
            if current==528 and number==1:cpp=False
            if not cpp:continue
            current=number;gutters[number]=(pn,fitz.Rect(span['bbox']))
        elif span['color']==0x20344B and cpp and current is not None:
            for c in span['chars']:
                if not c['c'].isspace():
                    by_line[current].append((pn,round(c['origin'][1],1),c['c'],fitz.Rect(c['bbox'])))

def add(chars):
    rows=defaultdict(list)
    for pn,y,_,r in chars:rows[(pn,y)].append(r)
    for (pn,_),rects in rows.items():
        rect=fitz.Rect(rects[0])
        for r in rects[1:]:rect|=r
        page=pdf[pn]
        annot=page.add_highlight_annot(fitz.Quad(rect.tl,rect.tr,rect.bl,rect.br))
        annot.set_colors(stroke=settings['color']);annot.set_opacity(settings['opacity']);annot.update()

for number in settings['code_lines']:
    found=''.join(c[2] for c in by_line[number])
    expected=re.sub(r'\s+','',source[number-1])
    assert found==expected,(number,found,expected)
    if found:add(by_line[number])
    else:
        pn,r=gutters[number]
        add([(pn,0,'',r)])
for selection in settings['spans']:
    number=selection['line'];chars=by_line[number];text=''.join(c[2] for c in chars)
    expected=re.sub(r'\s+','',source[number-1]);assert text==expected,(number,text,expected)
    needle=re.sub(r'\s+','',selection['text']);start=text.index(needle)
    add(chars[start:start+len(needle)])
assert not any(a.type[1]=='StrikeOut' for page in pdf for a in page.annots() or [])
path=root/'readme.styled.tmp.pdf';pdf.save(path,garbage=4,deflate=True);pdf.close();path.replace(root/'readme.pdf')
print('Restored blue selections on',len(settings['code_lines']),'code lines and',len(settings['spans']),'partial spans.')
