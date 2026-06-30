#!/usr/bin/env python3
"""Generate compact dependency-light SVG figures for the June 2026 ECEE results."""
from __future__ import annotations
from pathlib import Path
from xml.sax.saxutils import escape
import math
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / "figures"
OUT.mkdir(exist_ok=True)

def header(w,h,title):
    return [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
            '<rect width="100%" height="100%" fill="white"/>',
            f'<text x="{w/2:.1f}" y="24" text-anchor="middle" font-family="sans-serif" font-size="16">{escape(title)}</text>']

def footer(): return ['</svg>']

def axes(parts,w,h,xlabel,ylabel):
    l,r,t,b=75,w-25,45,h-65
    parts += [f'<line x1="{l}" y1="{b}" x2="{r}" y2="{b}" stroke="black"/>',
              f'<line x1="{l}" y1="{b}" x2="{l}" y2="{t}" stroke="black"/>',
              f'<text x="{(l+r)/2}" y="{h-18}" text-anchor="middle" font-family="sans-serif" font-size="12">{escape(xlabel)}</text>',
              f'<text x="18" y="{(t+b)/2}" text-anchor="middle" font-family="sans-serif" font-size="12" transform="rotate(-90 18 {(t+b)/2})">{escape(ylabel)}</text>']
    return l,r,t,b

def scale(v,a,b,u,z):
    return u+(v-a)*(z-u)/max(b-a,1e-12)

def score_color(s):
    s=max(0,min(1,float(s))); g=int(235-150*s); b=int(245-35*s)
    return f'rgb({int(55+40*(1-s))},{g},{b})'

df=pd.read_csv(HERE/'nasa_ecee_usable_battery_candidates_compact.csv')
top=pd.read_csv(HERE/'nasa_ecee_usable_battery_top50_compact.csv')

# mass-radius
w,h=900,650; p=header(w,h,'Usable rocky-battery search space: mass-radius')
l,r,t,b=axes(p,w,h,'Planet radius / Earth','Planet mass / Earth')
xmin,xmax=df.planet_radius_earth.min()*.95,df.planet_radius_earth.max()*1.05
ymin,ymax=0,df.planet_mass_earth.max()*1.08
for val in range(0,int(ymax)+1,5):
    y=scale(val,ymin,ymax,b,t); p += [f'<line x1="{l-4}" y1="{y:.1f}" x2="{r}" y2="{y:.1f}" stroke="#ddd"/>',f'<text x="{l-8}" y="{y+4:.1f}" text-anchor="end" font-family="sans-serif" font-size="10">{val}</text>']
for _,q in df.iterrows():
    x=scale(q.planet_radius_earth,xmin,xmax,l,r); y=scale(q.planet_mass_earth,ymin,ymax,b,t)
    p.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{score_color(q.battery_coherence_usable_score)}" stroke="#333" stroke-width=".4"/>')
for _,q in df.nlargest(10,'battery_coherence_usable_score').iterrows():
    x=scale(q.planet_radius_earth,xmin,xmax,l,r); y=scale(q.planet_mass_earth,ymin,ymax,b,t)
    p.append(f'<text x="{x+5:.1f}" y="{y-5:.1f}" font-family="sans-serif" font-size="9">{escape(q["name"])}</text>')
p += footer(); (OUT/'nasa_ecee_mass_radius_candidates.svg').write_text('\n'.join(p))

# core-to-stellar
w,h=900,650; p=header(w,h,'Outer/low-flux usable-battery candidates')
l,r,t,b=axes(p,w,h,'Insolation / Earth (log)','Core-to-stellar ratio / Earth (log)')
core_proxy=(df.planet_mass_earth.pow(1.2)/df.planet_radius_earth.pow(2))/df.insolation_earth.clip(lower=1e-4)
xvals=df.insolation_earth.clip(lower=1e-4).map(math.log10); yvals=core_proxy.clip(lower=1e-4).map(math.log10)
xmin,xmax=xvals.min()-.1,xvals.max()+.1; ymin,ymax=yvals.min()-.1,yvals.max()+.1
for _,q in df.iterrows():
    x=scale(math.log10(max(q.insolation_earth,1e-4)),xmin,xmax,l,r); y=scale(math.log10(max((q.planet_mass_earth**1.2/q.planet_radius_earth**2)/max(q.insolation_earth,1e-4),1e-4)),ymin,ymax,b,t)
    p.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="4" fill="{score_color(q.usable_battery_score)}" stroke="#333" stroke-width=".4"/>')
for _,q in df.nlargest(10,'usable_battery_score').iterrows():
    x=scale(math.log10(max(q.insolation_earth,1e-4)),xmin,xmax,l,r); y=scale(math.log10(max((q.planet_mass_earth**1.2/q.planet_radius_earth**2)/max(q.insolation_earth,1e-4),1e-4)),ymin,ymax,b,t)
    p.append(f'<text x="{x+5:.1f}" y="{y-5:.1f}" font-family="sans-serif" font-size="9">{escape(q["name"])}</text>')
p += footer(); (OUT/'nasa_ecee_core_to_stellar_candidates.svg').write_text('\n'.join(p))

# top25 bars
plot=top.nlargest(25,'battery_coherence_usable_score').sort_values('battery_coherence_usable_score')
w,h=900,720; p=header(w,h,'NASA Archive ECEM/ECEE usable-battery candidates - top 25')
l,r,t,b=105,w-35,45,h-55
rowh=(b-t)/len(plot)
for i,(_,q) in enumerate(plot.iterrows()):
    y=b-(i+1)*rowh+2; bw=(r-l)*float(q.battery_coherence_usable_score)
    p += [f'<text x="{l-8}" y="{y+rowh*.65:.1f}" text-anchor="end" font-family="sans-serif" font-size="9">{escape(q["name"])}</text>',
          f'<rect x="{l}" y="{y:.1f}" width="{bw:.1f}" height="{max(2,rowh-4):.1f}" fill="{score_color(q.battery_coherence_usable_score)}"/>',
          f'<text x="{l+bw+4:.1f}" y="{y+rowh*.65:.1f}" font-family="sans-serif" font-size="9">{q.battery_coherence_usable_score:.3f}</text>']
p += [f'<line x1="{l}" y1="{b}" x2="{r}" y2="{b}" stroke="black"/>',f'<text x="{(l+r)/2}" y="{h-17}" text-anchor="middle" font-family="sans-serif" font-size="12">Battery coherence usable score</text>']+footer()
(OUT/'nasa_ecee_usable_battery_top25.svg').write_text('\n'.join(p))
