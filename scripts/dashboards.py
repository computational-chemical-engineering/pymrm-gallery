#!/usr/bin/env python3
"""Generate the two maintainer dashboards from the case queue.

    python scripts/dashboards.py

Writes docs/dashboards/papers-needed.html and docs/dashboards/needs-input.html.
Both are self-contained (images inlined) so they can be published as artifacts.
"""
from __future__ import annotations

import base64
import html
import json
from collections import Counter
from datetime import date
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
QUEUE = ROOT / "queue_cases"
OUT = ROOT / "docs" / "dashboards"

SECTION_NAME = {
    "A": "Foundations — transport closures",
    "B": "Particles & pellets",
    "C": "Kinetics",
    "D": "Fixed-bed reactors",
    "E": "Fluidised beds",
    "F": "Gas–liquid reactors",
    "G": "Trickle beds & multiphase",
    "H": "Membranes",
    "I": "Environmental & automotive",
    "J": "Adjacent fields",
}

CSS = """
:root{
  --bg:#F5F6F8; --panel:#FFF; --ink:#15181D; --muted:#5A6472; --line:#E0E4EA;
  --accent:#9A5B18; --accent-soft:#FBF0E3; --ok:#2F7D5A; --stop:#98362F;
  --shadow:0 1px 2px rgba(16,20,28,.05),0 8px 22px rgba(16,20,28,.05);
  --on-accent:#FFF;
}
@media (prefers-color-scheme:dark){
 :root{--bg:#0F1216;--panel:#161A20;--ink:#E6E9EE;--muted:#95A0AE;--line:#252B34;
   --accent:#E0A264;--accent-soft:#2A2015;--ok:#68C39A;--stop:#E08078;
   --shadow:0 1px 2px rgba(0,0,0,.4),0 8px 22px rgba(0,0,0,.3);--on-accent:#1A130A;}
}
:root[data-theme="dark"]{--bg:#0F1216;--panel:#161A20;--ink:#E6E9EE;--muted:#95A0AE;
  --line:#252B34;--accent:#E0A264;--accent-soft:#2A2015;--ok:#68C39A;--stop:#E08078;
  --shadow:0 1px 2px rgba(0,0,0,.4),0 8px 22px rgba(0,0,0,.3);--on-accent:#1A130A;}
:root[data-theme="light"]{--bg:#F5F6F8;--panel:#FFF;--ink:#15181D;--muted:#5A6472;
  --line:#E0E4EA;--accent:#9A5B18;--accent-soft:#FBF0E3;--ok:#2F7D5A;--stop:#98362F;
  --shadow:0 1px 2px rgba(16,20,28,.05),0 8px 22px rgba(16,20,28,.05);--on-accent:#FFF;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
 font:16px/1.6 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
 -webkit-font-smoothing:antialiased}
.wrap{max-width:1100px;margin:0 auto;padding:0 20px 90px}
code,.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
header.top{padding:52px 0 22px;border-bottom:1px solid var(--line)}
.eyebrow{font-size:12px;letter-spacing:.15em;text-transform:uppercase;color:var(--accent);font-weight:650}
h1{font-size:clamp(27px,4.2vw,40px);line-height:1.12;letter-spacing:-.02em;margin:11px 0 12px;text-wrap:balance}
.lede{font-size:17px;color:var(--muted);max-width:66ch;margin:0}
.lede strong{color:var(--ink);font-weight:620}
.stats{display:flex;flex-wrap:wrap;gap:10px;margin:24px 0 0;padding:0;list-style:none}
.stat{background:var(--panel);border:1px solid var(--line);border-radius:11px;
 padding:11px 15px;min-width:104px;box-shadow:var(--shadow)}
.stat b{display:block;font-size:23px;line-height:1.15;font-variant-numeric:tabular-nums;letter-spacing:-.01em}
.stat span{font-size:12px;color:var(--muted);letter-spacing:.02em}
.stat.hot b{color:var(--accent)}
h2{font-size:14px;letter-spacing:.09em;text-transform:uppercase;color:var(--muted);
 margin:38px 0 12px;font-weight:700}
.card{background:var(--panel);border:1px solid var(--line);border-radius:13px;
 box-shadow:var(--shadow);margin:0 0 12px;overflow:hidden}
.card-h{display:flex;flex-wrap:wrap;align-items:baseline;gap:9px;padding:15px 18px 0}
.cid{font-family:ui-monospace,Menlo,Consolas,monospace;font-weight:650;font-size:12.5px;
 background:var(--accent-soft);color:var(--accent);border-radius:6px;padding:3px 8px}
.ctitle{font-size:16.5px;font-weight:640;letter-spacing:-.01em}
.pill{font-size:11px;font-weight:620;border:1px solid var(--line);border-radius:999px;
 padding:2px 9px;color:var(--muted)}
.pill.p1{border-color:var(--accent);color:var(--accent)}
.card-b{padding:9px 18px 17px;font-size:14.6px;color:var(--muted)}
.card-b .ref{color:var(--ink)}
.q{color:var(--ink);font-size:15.4px;font-weight:600;margin:10px 0 8px;max-width:76ch}
.detail{white-space:pre-wrap;font-size:14.3px;margin:0 0 4px;max-width:82ch}
a{color:var(--accent)}
.shots{display:grid;gap:10px;margin:13px 0 2px}
.shots img{width:100%;height:auto;border:1px solid var(--line);border-radius:9px;background:#fff}
.empty{background:var(--panel);border:1px dashed var(--line);border-radius:13px;
 padding:34px 22px;text-align:center;color:var(--muted)}
.empty b{display:block;color:var(--ink);font-size:18px;margin-bottom:6px}
table{width:100%;border-collapse:collapse;font-size:14.4px}
th{text-align:left;font-size:11.5px;letter-spacing:.07em;text-transform:uppercase;
 color:var(--muted);padding:9px 12px;border-bottom:1px solid var(--line);font-weight:700}
td{padding:9px 12px;border-bottom:1px solid var(--line);vertical-align:top}
tr:last-child td{border-bottom:0}
.scroll{overflow-x:auto}
.tools{display:flex;gap:9px;flex-wrap:wrap;margin:16px 0 0}
.btn{background:var(--accent);color:var(--on-accent);border:0;border-radius:9px;
 padding:10px 16px;font:inherit;font-weight:650;font-size:14px;cursor:pointer}
.btn.ghost{background:transparent;color:var(--accent);border:1px solid var(--line)}
pre#dump{white-space:pre-wrap;background:var(--bg);border:1px solid var(--line);
 border-radius:9px;padding:13px;font-family:ui-monospace,Menlo,Consolas,monospace;
 font-size:12.4px;max-height:300px;overflow:auto;margin:14px 0 0}
textarea{width:100%;min-height:66px;resize:vertical;background:var(--bg);color:var(--ink);
 border:1px solid var(--line);border-radius:9px;padding:10px 11px;font:inherit;font-size:14.3px}
footer{color:var(--muted);font-size:13px;margin-top:36px;border-top:1px solid var(--line);padding-top:16px}
@media (prefers-reduced-motion:reduce){*{transition:none!important}}
@media (max-width:640px){.wrap{padding:0 14px 70px}.card-h,.card-b{padding-left:14px;padding-right:14px}}
"""


def esc(s):
    return html.escape(str(s or ""))


def entries():
    out = []
    for p in sorted(QUEUE.glob("*.yaml")):
        try:
            out.append(yaml.safe_load(p.read_text(encoding="utf-8")))
        except Exception:
            pass
    return [e for e in out if e]


def img_uri(path: Path):
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode()


def ref_line(e):
    r = e.get("reference") or {}
    bits = []
    if r.get("authors"):
        a = r["authors"]
        bits.append(a[0] + (" et al." if len(a) > 2 else (" & " + a[1] if len(a) == 2 else "")))
    if r.get("year"):
        bits.append(f"({r['year']})")
    if r.get("container"):
        bits.append(r["container"])
    return ", ".join(bits) if bits else (e.get("catalog_reference") or "")


def doi_link(e):
    d = (e.get("reference") or {}).get("doi")
    if not d:
        return ""
    d = str(d).replace("https://doi.org/", "")
    return f' &middot; <a href="https://doi.org/{esc(d)}" target="_blank" rel="noopener">doi:{esc(d)}</a>'


# ---------------------------------------------------------------- papers page
def papers_page(es):
    need = [e for e in es if e["status"] == "needs-paper"]
    unclaimed = [e for e in es if e["status"] == "unclaimed"]
    done = [e for e in es if e["status"] == "published"]
    p1 = [e for e in need if e.get("priority") == "P1"]

    P = ['<div class="wrap">', """<header class="top">
 <div class="eyebrow">pymrm-gallery &middot; literature</div>
 <h1>Papers I need from you</h1>
 <p class="lede">Each of these is a case an agent reached and could not proceed on, because the
 source is neither on disk nor reachable as open access. <strong>Everything else is being worked
 automatically.</strong> Drop PDFs into <code>~/papers/pymrm-gallery/</code> — filename does not
 matter, they get matched by DOI and metadata.</p>
 <ul class="stats">"""]
    for lab, n, hot in (("blocked on a PDF", len(need), True), ("of those, P1", len(p1), True),
                        ("not yet reached", len(unclaimed), False), ("published", len(done), False),
                        ("catalogued", len(es), False)):
        P.append(f'<li class="stat{" hot" if hot and n else ""}"><b>{n}</b><span>{lab}</span></li>')
    P.append("</ul>")
    if need:
        P.append('<div class="tools"><button class="btn" id="copy">Copy the list</button>'
                 '<button class="btn ghost" id="dois">Copy DOIs only</button></div>')
    P.append("</header>")

    if not need:
        P.append('<div class="empty" style="margin-top:30px"><b>Nothing is waiting on you.</b>'
                 'No case has hit a paper wall yet. This page fills in as agents work through the '
                 'catalogue — check back, or watch the other dashboard for decisions.</div>')
    else:
        for sec in sorted({e["section"] for e in need}):
            grp = sorted([e for e in need if e["section"] == sec],
                         key=lambda e: (e.get("priority", "P3"), e["id"]))
            P.append(f'<h2>{esc(sec)} &mdash; {esc(SECTION_NAME.get(sec, ""))} '
                     f'<span style="color:var(--muted);font-weight:600">({len(grp)})</span></h2>')
            for e in grp:
                pr = e.get("priority", "")
                P.append('<div class="card"><div class="card-h">'
                         f'<span class="cid">{esc(e["id"])}</span>'
                         f'<span class="ctitle">{esc(e["title"])}</span>'
                         + (f'<span class="pill {"p1" if pr=="P1" else ""}">{esc(pr)}</span>' if pr else "")
                         + '</div><div class="card-b">'
                         f'<span class="ref">{esc(ref_line(e))}</span>{doi_link(e)}')
                b = e.get("blocker") or {}
                if b.get("detail"):
                    P.append(f'<p class="detail">{esc(b["detail"])}</p>')
                P.append("</div></div>")
    P.append('<pre id="dump" hidden></pre>')
    P.append(f'<footer>Generated {date.today()} from <code>queue_cases/</code>. '
             f'{len(es)} catalogued cases.</footer></div>')

    payload = [dict(id=e["id"], title=e["title"], ref=ref_line(e),
                    doi=(e.get("reference") or {}).get("doi", ""),
                    priority=e.get("priority", "")) for e in need]
    js = """
const NEED=%s;
function show(t){const d=document.getElementById('dump');d.hidden=false;d.textContent=t;
  navigator.clipboard&&navigator.clipboard.writeText(t).catch(()=>{});}
const cp=document.getElementById('copy'), dz=document.getElementById('dois');
if(cp)cp.addEventListener('click',()=>show(NEED.map(e=>
  `${e.id}  ${e.title}\\n    ${e.ref}${e.doi?'\\n    doi:'+e.doi:''}`).join('\\n\\n')||'(none)'));
if(dz)dz.addEventListener('click',()=>show(NEED.filter(e=>e.doi).map(e=>e.doi).join('\\n')||'(no DOIs known)'));
""" % json.dumps(payload)
    return ("<title>Papers needed — pymrm-gallery</title>\n<style>" + CSS + "</style>\n"
            + "\n".join(P) + "\n<script>" + js + "</script>\n")


# ----------------------------------------------------------------- input page
def input_page(es):
    blocked = [e for e in es if e["status"] == "needs-input"]
    ready = [e for e in es if e["status"] == "ready"]
    prog = [e for e in es if e["status"] == "in-progress"]
    done = [e for e in es if e["status"] == "published"]

    P = ['<div class="wrap">', """<header class="top">
 <div class="eyebrow">pymrm-gallery &middot; decisions</div>
 <h1>Waiting on your judgement</h1>
 <p class="lede">Cases an agent built up to the point where it needed something only you can give —
 a visual check on a figure extraction, or a call about scope or a suspect constant.
 <strong>Each one is halted, not abandoned</strong>: answer any subset and the rest keep moving.
 Use the button at the end to copy your answers back.</p>
 <ul class="stats">"""]
    for lab, n, hot in (("need you", len(blocked), True), ("built, awaiting merge", len(ready), False),
                        ("in progress", len(prog), False), ("published", len(done), False)):
        P.append(f'<li class="stat{" hot" if hot and n else ""}"><b>{n}</b><span>{lab}</span></li>')
    P.append("</ul></header>")

    if not blocked:
        P.append('<div class="empty" style="margin-top:30px"><b>Nothing needs you right now.</b>'
                 'Agents halt here when they hit a judgement call — a figure extraction to eyeball, '
                 'a constant that looks wrong, a scope decision. Nothing is stuck silently.</div>')
    else:
        for e in sorted(blocked, key=lambda e: (e.get("priority", "P3"), e["id"])):
            b = e.get("blocker") or {}
            kind = b.get("kind", "decision")
            P.append('<div class="card"><div class="card-h">'
                     f'<span class="cid">{esc(e["id"])}</span>'
                     f'<span class="ctitle">{esc(e["title"])}</span>'
                     f'<span class="pill">{esc(kind)}</span></div><div class="card-b">'
                     f'<span class="ref">{esc(ref_line(e))}</span>{doi_link(e)}')
            if b.get("question"):
                P.append(f'<p class="q">{esc(b["question"])}</p>')
            if b.get("detail"):
                P.append(f'<p class="detail">{esc(b["detail"])}</p>')
            shots = []
            for a in (b.get("artifacts") or []):
                ap = (ROOT / a) if not str(a).startswith("/") else Path(a)
                if ap.is_file() and ap.suffix.lower() == ".png":
                    shots.append(ap)
            if shots:
                P.append('<div class="shots">')
                for sp in shots:
                    P.append(f'<img src="{img_uri(sp)}" alt="{esc(e["id"])} review">')
                P.append("</div>")
            P.append(f'<textarea id="a_{esc(e["id"]).replace(".","_")}" '
                     f'placeholder="Your answer for {esc(e["id"])} — free text is fine, and '
                     f'anything I did not think to ask is the most useful thing here."></textarea>')
            P.append("</div></div>")
        P.append('<div class="card"><div class="card-b">'
                 '<div class="tools"><button class="btn" id="copy">Copy answers</button></div>'
                 '<pre id="dump">Nothing yet.</pre></div></div>')
    P.append(f'<footer>Generated {date.today()} from <code>queue_cases/</code>. '
             f'Figures are shown only to check a data extraction; no page image is committed to the '
             f'public repository.</footer></div>')

    ids = [e["id"] for e in blocked]
    js = """
const IDS=%s;
const b=document.getElementById('copy');
if(b)b.addEventListener('click',()=>{
  const L=['CASE DECISIONS — pymrm-gallery',''];
  IDS.forEach(id=>{const t=document.getElementById('a_'+id.replace('.','_'));
    if(t&&t.value.trim()){L.push('## '+id);t.value.trim().split('\\n').forEach(x=>L.push('  '+x));L.push('');}});
  const out=L.length>2?L.join('\\n'):'(nothing filled in yet)';
  document.getElementById('dump').textContent=out;
  navigator.clipboard&&navigator.clipboard.writeText(out).catch(()=>{});
});
""" % json.dumps(ids)
    return ("<title>Needs your input — pymrm-gallery</title>\n<style>" + CSS + "</style>\n"
            + "\n".join(P) + "\n<script>" + js + "</script>\n")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    es = entries()
    (OUT / "papers-needed.html").write_text(papers_page(es), encoding="utf-8")
    (OUT / "needs-input.html").write_text(input_page(es), encoding="utf-8")
    c = Counter(e["status"] for e in es)
    print(f"{len(es)} cases: " + ", ".join(f"{k}={v}" for k, v in sorted(c.items())))
    for f in ("papers-needed.html", "needs-input.html"):
        print(f"   docs/dashboards/{f}  {(OUT/f).stat().st_size//1024} kB")
