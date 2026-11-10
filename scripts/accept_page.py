#!/usr/bin/env python3
"""Run the refresh acceptance checks on one page and print a short summary.

    python scripts/accept_page.py pages/B2.1-voorhies-coking-law
    python scripts/accept_page.py pages/B2.1-voorhies-coking-law --quick
    python scripts/accept_page.py pages/B2.1-voorhies-coking-law --only lint
    python scripts/accept_page.py pages/A3.4-wakao-funazkri --make-ref <scratch>/a34_old
    python scripts/accept_page.py pages/A3.4-wakao-funazkri --quick \
        --only sync,kernels,port,stable --ref <scratch>/a34_old

One line per check (PASS, FAIL, SKIP or LOOK), at most a few detail lines under a
failure, about 30 lines in all. Everything the checks print goes to the log file
(``--log``, default ``accept_<ID>.log`` in the system temp directory), so an agent
reads the summary and opens the log only at a failure, by line range.

Checks, in order (``docs/playbooks`` and the gallery-refresh skill call them the
acceptance steps):

sync      ``jupytext --sync index.py``
kernels   ``scripts/run_kernels.py`` (three kernels); with ``--quick`` one pinned run
          instead, for iterating. The final acceptance is never ``--quick``.
          Prints the page's claim count and its break-table line, and checks the
          break-row count: FAIL when a row the page declares
          (``claims.define_rows``) recorded nothing, LOOK when it declares none.
agreement ``scripts/check_agreement.py --base <base>``: this page's lines that are
          not OK or NEW
claims    ``scripts/check_claims.py``, ``pytest tests/`` and the page's
          ``test_model.py`` offline
plugin    ``check_notebook_math.py`` on index.ipynb and ``check_model.py`` on
          model.py (agreement.json is restored afterwards)
metadata  ``scripts/check_metadata.py``
stable    a second ``jupytext --sync`` leaves index.py unchanged
dashes    em dashes in the page files, its models.yaml entry and the learnings rows
          naming it, and pymrm version names in the page code and meta.yaml (FAIL);
          history words there (LOOK: "the previous step" is fine, "the previous
          version of this page" is not)
lint      the refresh checklist's "Before you report" scans: Correction quotes
          against the old page (FAIL on a miss) and untagged README numbers (LOOK:
          each hit must be tagged, quoted or an identifier)
outputs   only with ``--only``: the claim count, break-table line and break-row
          count from the executed index.ipynb and agreement.json as they stand,
          without running anything
port      only with ``--only`` and ``--ref``: a port onto the new layout must compute
          exactly what the old page computed. FAIL unless every metric of the
          page's agreement.json equals the reference run's value exactly (the same
          pinned kernel and venv); LOOK for printed lines that differ from the
          reference run's (cell numbers ignored). check_agreement is no substitute:
          a format-1 baseline allows 5 % drift and comes from an older pymrm.

``--make-ref DIR`` runs the OLD page (its build_page.py) once on the pinned Haswell
kernel and keeps the notebook and agreement.json in DIR as ``Haswell.ipynb`` and
``Haswell.json``; the page directory must be clean in git and is restored after.
Run it before the port starts.

Exit status 1 when any check fails.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLUGIN = Path.home() / "Code/pymrm_suite/pymrm/plugins/pymrm/skills/conventions/scripts"
STEPS = ("sync", "kernels", "agreement", "claims", "plugin", "metadata", "stable",
         "dashes", "lint")
EXTRA = ("outputs", "port")
HISTORY = re.compile(r"\bearlier\b|\bprevious|\bused to\b|before the refresh")
VERSION = re.compile(r"pymrm [0-9]")


class Runner:
    def __init__(self, log: Path):
        self.log = log.open("w")
        self.failed = False
        self.lines = 0

    def run(self, cmd: list[str], cwd: Path = ROOT, env: dict | None = None
            ) -> tuple[int, str]:
        """Run a command, log its output, return (status, output)."""
        import os
        full_env = dict(os.environ, **(env or {}))
        self.log.write(f"\n$ {' '.join(map(str, cmd))}  (in {cwd})\n")
        self.log.flush()
        proc = subprocess.run(cmd, cwd=cwd, env=full_env, capture_output=True,
                              text=True)
        out = proc.stdout + proc.stderr
        self.log.write(out)
        self.log.write(f"[exit {proc.returncode}]\n")
        self.log.flush()
        return proc.returncode, out

    def report(self, step: str, status: str, text: str = "",
               details: list[str] | None = None) -> None:
        if status == "FAIL":
            self.failed = True
        print(f"{status:<4} {step:<9} {text}".rstrip())
        for d in (details or [])[:4]:
            print(f"       {d[:160]}")
        if details and len(details) > 4:
            print(f"       ... {len(details) - 4} more in the log")
        sys.stdout.flush()


def tail(out: str, n: int = 3) -> list[str]:
    """The last lines that look like an error, else the last lines."""
    lines = [ln for ln in out.strip().splitlines() if ln.strip()]
    errs = [ln for ln in lines if re.search(r"Error|FAIL|error:|Traceback|assert", ln)]
    return (errs or lines)[-n:]


def models_entry(cat_id: str) -> str:
    text = (ROOT / "models.yaml").read_text()
    m = re.search(rf"^  - id: {re.escape(cat_id)}\n.*?(?=^  - id: |\Z)", text,
                  re.M | re.S)
    return m.group(0) if m else ""


def learnings_rows(cat_id: str) -> list[str]:
    path = ROOT / "docs/refresh-learnings.md"
    if not path.is_file():
        return []
    return [ln for ln in path.read_text().splitlines()
            if ln.startswith("|") and f"`{cat_id}`" in ln]


def notebook_text(nb_path: Path) -> str:
    try:
        nb = json.loads(nb_path.read_text())
    except (OSError, ValueError):
        return ""
    parts = []
    for cell in nb.get("cells", []):
        for o in cell.get("outputs", []):
            if "text" in o:
                parts.append("".join(o["text"]))
            elif "data" in o and "text/plain" in o["data"]:
                parts.append("".join(o["data"]["text/plain"]))
    return "\n".join(parts)


def summary_of_outputs(page: Path, nb_path: Path | None = None) -> list[str]:
    text = notebook_text(nb_path or page / "index.ipynb")
    out = []
    m = re.search(r"Claims for [^\n]*\n((?:  \S[^\n]*\n)+)", text)
    if m:
        out.append(f"{len(m.group(1).splitlines())} claims printed by claims.report()")
    m = re.search(r"[^\n]*claims are failed by at least one row[^\n]*", text)
    if m:
        out.append(m.group(0).strip())
    return out


def check_rows(r: Runner, agreement: Path) -> None:
    """Report the break-row count a run recorded in ``agreement.json``.

    ``Claims.report`` raises for a declared row that recorded no claim, so a
    mismatch here means the file is left by a run that failed; a page that declares no rows cannot
    show a row that silently dropped out of its table.
    """
    try:
        data = json.loads(agreement.read_text())
    except (OSError, ValueError):
        r.report("rows", "FAIL", f"no readable {agreement.name}")
        return
    if not isinstance(data, dict):
        r.report("rows", "FAIL", f"{agreement.name} is not a JSON object")
        return
    if data.get("format") != 2:
        r.report("rows", "SKIP", "agreement.json is format 1 (report_agreement)")
        return
    rows = data.get("rows")
    if not rows:
        r.report("rows", "LOOK", "break rows not declared: call claims.define_rows(labels) "
                 "before running them, so a row that records nothing fails the report")
        return
    ok = rows["defined"] == rows["recorded"]
    r.report("rows", "PASS" if ok else "FAIL",
             f"break rows: {rows['defined']} defined, {rows['recorded']} recorded")


def git_show(path: str, base: str) -> str:
    proc = subprocess.run(["git", "show", f"{base}:{path}"], cwd=ROOT,
                          capture_output=True, text=True)
    return proc.stdout if proc.returncode == 0 else ""


def lint(page: Path, base: str) -> tuple[list[str], list[str]]:
    """Correction quotes missing from the old page; untagged README numbers."""
    readme = (page / "README.md").read_text() if (page / "README.md").is_file() else ""
    rel = page.relative_to(ROOT).as_posix()
    old = ""
    nb = git_show(f"{rel}/index.ipynb", base)
    if nb:
        try:
            old = " ".join("".join(c["source"]) for c in json.loads(nb)["cells"])
        except ValueError:
            pass
    old += " ".join(git_show(f"{rel}/{f}", base)
                    for f in ("README.md", "build_page.py", "meta.yaml"))
    old += git_show("models.yaml", base)

    def norm(s: str) -> str:
        return re.sub(r"\s+", " ", s.replace("*", ""))

    old = norm(old)
    misses = []
    if "## Corrections" in readme:
        corr = readme.split("## Corrections", 1)[1].split("\n## ")[0]
        for q in re.findall(r'"([^"]+)"|“([^”]+)”', corr):
            q = norm("".join(q))
            if q not in old:
                misses.append(f"quote not in the old page: {q[:110]}")

    text = re.sub(r'"[^"]*"|“[^”]*”|`[^`]*`|\]\([^)]*\)|doi\S*',
                  lambda m: "\n" * m.group().count("\n"), readme)
    num = (r'(?<![\w.\-])\d+(?:\.\d+)?(?:e-?\d+)?(?:\s*%|x)?(?![\w.]*\d)'
           r'(?!\s*%?\s*x?\s*<!--\s*claim)')
    untagged = []
    for i, line in enumerate(text.splitlines(), 1):
        hits = [h for h in re.findall(num, line) if re.search(r"[.e%x]", h.rstrip("."))]
        if hits:
            untagged.append(f"README.md:{i} {hits}")
    return misses, untagged


def make_ref(page: Path, out: Path) -> int:
    """Run the old page once on the pinned kernel; keep the result in ``out``."""
    rel = page.relative_to(ROOT).as_posix()
    if not (page / "build_page.py").is_file():
        print(f"FAIL {rel}: no build_page.py, so no old page to run")
        return 1
    dirty = subprocess.run(["git", "status", "--porcelain", "--", rel], cwd=ROOT,
                           capture_output=True, text=True).stdout.strip()
    if dirty:
        print(f"FAIL {rel} has uncommitted changes; --make-ref restores the page "
              "directory afterwards, so it must start clean")
        return 1
    sys.path.insert(0, str(ROOT / "scripts"))
    import run_kernels                                           # noqa: E402
    out.mkdir(parents=True, exist_ok=True)
    try:
        run_kernels.build(page)
        ok, secs, msg = run_kernels.execute(page, "Haswell", out)
    finally:
        subprocess.run(["git", "checkout", "--", rel], cwd=ROOT, check=True)
        subprocess.run(["git", "clean", "-fdq", "--", rel], cwd=ROOT, check=True)
    print(f"{'PASS' if ok else 'FAIL'} reference run of the old page, {secs:.0f} s, "
          f"in {out} {msg}".rstrip())
    return 0 if ok else 1


def port_check(r: Runner, page: Path, ref: Path | None, executed: Path) -> None:
    """Exact comparison of a ported page with the reference run of the old page."""
    if not ref or not (ref / "Haswell.json").is_file():
        r.report("port", "FAIL", "needs --ref DIR made by --make-ref before the port")
        return
    old = json.loads((ref / "Haswell.json").read_text()).get("metrics", {})
    try:
        new = json.loads((page / "agreement.json").read_text()).get("metrics", {})
    except (OSError, ValueError):
        new = {}
    diffs = []
    for k, v in old.items():
        w = new.get(k)
        w = w.get("value") if isinstance(w, dict) else w
        same = (w == v) or (isinstance(v, float) and isinstance(w, float)
                            and v != v and w != w)
        if not same:
            diffs.append(f"{k}: old {v!r}, port {w!r}")
    diffs += [f"{k}: not in the old page" for k in new if k not in old]
    r.report("port", "FAIL" if diffs else "PASS",
             f"{len(old)} metrics compared exactly with the old page's run", diffs)
    sys.path.insert(0, str(ROOT / "scripts"))
    import run_kernels                                           # noqa: E402

    def lines(nb: Path) -> list[str]:
        return [re.sub(r"^cell \d+: ", "", ln) for ln in run_kernels.printed(nb)]

    if executed.is_file():
        a, b = lines(ref / "Haswell.ipynb"), lines(executed)
        gone = [f"old only: {ln}" for ln in a if ln not in b]
        new_ = [f"port only: {ln}" for ln in b if ln not in a]
        if gone or new_:
            r.log.write("\nport, printed lines:\n" + "\n".join(gone + new_) + "\n")
            r.report("port", "LOOK", f"{len(gone)} printed line(s) of the old page "
                     f"missing, {len(new_)} new: each must be a moved cell or the "
                     "bootstrap, not a changed result", gone + new_)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("page", type=Path)
    ap.add_argument("--quick", action="store_true",
                    help="one pinned run instead of three kernels (iteration only)")
    ap.add_argument("--base", default="origin/main")
    ap.add_argument("--only", help="comma-separated subset of: "
                    + ",".join(STEPS + EXTRA))
    ap.add_argument("--log", type=Path)
    ap.add_argument("--ref", type=Path, help="reference run directory for the port step")
    ap.add_argument("--make-ref", type=Path, metavar="DIR",
                    help="run the old page once and keep the result in DIR")
    args = ap.parse_args()

    page = (ROOT / args.page).resolve() if not args.page.is_absolute() else args.page
    if args.make_ref:
        return make_ref(page, args.make_ref)
    if not (page / "index.py").is_file():
        print(f"FAIL {page}: no index.py (not a refreshed page)")
        return 1
    cat_id = page.name.split("-", 1)[0]
    rel = page.relative_to(ROOT).as_posix()
    log = args.log or Path(tempfile.gettempdir()) / f"accept_{cat_id}.log"
    steps = args.only.split(",") if args.only else list(STEPS)
    unknown = set(steps) - set(STEPS + EXTRA)
    if unknown:
        ap.error(f"unknown step(s): {', '.join(sorted(unknown))}")
    r = Runner(log)
    py = sys.executable
    t0 = time.time()
    print(f"accept {rel} ({'quick' if args.quick else 'full'}), log {log}")

    executed = page / "index.ipynb"     # the executed notebook the checks read
    if "sync" in steps:
        code, out = r.run([py, "-m", "jupytext", "--sync", "index.py"], cwd=page)
        r.report("sync", "PASS" if code == 0 else "FAIL", "", None if code == 0 else tail(out))

    if "kernels" in steps:
        if args.quick:
            sys.path.insert(0, str(ROOT / "scripts"))
            import run_kernels                                   # noqa: E402
            run_kernels.build(page)
            # the kernel's own stderr (a transport warning) goes to the log
            import os
            sys.stderr.flush()
            saved_err = os.dup(2)
            os.dup2(r.log.fileno(), 2)
            out_dir = Path(tempfile.mkdtemp())
            try:
                ok, secs, msg = run_kernels.execute(page, "Haswell", out_dir)
            finally:
                os.dup2(saved_err, 2)
                os.close(saved_err)
            r.log.write(f"\nquick run: ok={ok} {secs:.0f} s {msg}\n")
            r.report("kernels", "PASS" if ok else "FAIL",
                     f"one pinned run, {secs:.0f} s (quick: not the acceptance)",
                     (summary_of_outputs(page, out_dir / "Haswell.ipynb") if ok
                      else [msg]))
            executed = out_dir / "Haswell.ipynb"
            if ok:
                check_rows(r, out_dir / "Haswell.json")
        else:
            code, out = r.run([py, "scripts/run_kernels.py", rel])
            executed = page / "index.ipynb"
            details = summary_of_outputs(page) if code == 0 else tail(out, 4)
            r.report("kernels", "PASS" if code == 0 else "FAIL",
                     "three kernels, claims and printed text agree" if code == 0 else "",
                     details)
            if code == 0:
                check_rows(r, page / "agreement.json")

    if "port" in steps:
        port_check(r, page, args.ref, executed)

    if "agreement" in steps:
        code, out = r.run([py, "scripts/check_agreement.py", "--base", args.base])
        mine = [ln for ln in out.splitlines() if page.name in ln]
        bad = [ln for ln in mine if not re.match(r"\s*(OK|NEW)\b", ln)]
        others = [ln for ln in out.splitlines()
                  if ln.startswith("FAIL ") and page.name not in ln]
        total = re.search(r"(\d+) regression\(s\)", out)
        n_reg = int(total.group(1)) if total else -1
        try:
            v1 = json.loads(git_show(f"{rel}/agreement.json", "HEAD") or "{}"
                            ).get("format") != 2
        except ValueError:
            v1 = False
        # Before the commit, a format-1 baseline shows drift and "disappeared" lines
        # for every moved or replaced metric (commit-and-merge.md): each is
        # explained, not hidden, so this page's failures are LOOK there.
        status = "PASS" if code == 0 else ("LOOK" if v1 and not others else "FAIL")
        r.report("agreement", status,
                 f"{n_reg} regression(s) in all pages; {len(bad)} line(s) for this page "
                 "to explain" + (" (format-1 baseline)" if v1 else ""), others + bad)

    if "claims" in steps:
        fails = []
        code, out = r.run([py, "scripts/check_claims.py"])
        if code:
            fails.append("check_claims: " + " / ".join(tail(out, 2)))
        code, out = r.run([py, "-m", "pytest", "-q", "tests/"])
        if code:
            fails.append("pytest tests/: " + " / ".join(tail(out, 2)))
        if (page / "test_model.py").is_file():
            code, out = r.run([py, "-m", "pytest", "-q", f"{rel}/test_model.py"],
                              env={"GALLERY_OFFLINE": "1"})
            if code:
                fails.append("test_model.py: " + " / ".join(tail(out, 2)))
        r.report("claims", "FAIL" if fails else "PASS",
                 "check_claims, tests/, test_model.py offline", fails)

    if "plugin" in steps:
        fails = []
        if not PLUGIN.is_dir():
            # a missing check is a failure in the acceptance, a skip in a subset
            r.report("plugin", "SKIP" if args.only else "FAIL",
                     f"no plugin scripts at {PLUGIN}")
        else:
            code, out = r.run([py, str(PLUGIN / "check_notebook_math.py"),
                               str(page / "index.ipynb")])
            if code:
                fails.append("check_notebook_math: " + " / ".join(tail(out, 2)))
            model = page / "model.py"
            if model.is_file():
                agr = page / "agreement.json"
                saved = agr.read_bytes() if agr.is_file() else None
                try:
                    code, out = r.run([py, str(PLUGIN / "check_model.py"), str(model)],
                                      cwd=page)
                finally:
                    if saved is not None:
                        agr.write_bytes(saved)
                if code:
                    fails.append("check_model: " + " / ".join(tail(out, 2)))
            r.report("plugin", "FAIL" if fails else "PASS",
                     "check_notebook_math" + (", check_model" if model.is_file() else ""),
                     fails)

    if "metadata" in steps:
        code, out = r.run([py, "scripts/check_metadata.py"])
        r.report("metadata", "PASS" if code == 0 else "FAIL", "", None if code == 0 else tail(out))

    if "stable" in steps:
        before = (page / "index.py").read_bytes()
        code, out = r.run([py, "-m", "jupytext", "--sync", "index.py"], cwd=page)
        same = (page / "index.py").read_bytes() == before
        r.report("stable", "PASS" if code == 0 and same else "FAIL",
                 "" if same else "jupytext --sync changed index.py: edit index.py, "
                 "never the notebook alone")

    if "dashes" in steps:
        hits, words = [], []
        for f in sorted(page.glob("*")):
            if f.suffix in {".py", ".md", ".yaml"}:
                for i, ln in enumerate(f.read_text().splitlines(), 1):
                    if "\u2014" in ln:
                        hits.append(f"{f.name}:{i} em dash")
                    if f.suffix == ".py" or f.name == "meta.yaml":
                        if VERSION.search(ln):
                            hits.append(f"{f.name}:{i} pymrm version name")
                        m = HISTORY.search(ln)
                        if m:
                            words.append(f"{f.name}:{i} {m.group(0)!r}: "
                                         f"{ln.strip()[:90]}")
        if "\u2014" in models_entry(cat_id):
            hits.append(f"models.yaml {cat_id} entry: em dash")
        for ln in learnings_rows(cat_id):
            if "\u2014" in ln:
                hits.append(f"learnings row {ln.split('|')[1].strip()}: em dash")
        r.log.write("\ndashes and history words:\n" + "\n".join(hits + words) + "\n")
        r.report("dashes", "FAIL" if hits else "PASS",
                 "em dashes, pymrm version names", hits)
        if words:
            # "the previous step" is fine, "the previous version of this page" is not
            r.report("dashes", "LOOK", f"{len(words)} history word(s): each must be "
                     "about the computation, never the page's or pymrm's history", words)

    if "lint" in steps:
        misses, untagged = lint(page, args.base)
        r.log.write("\nlint:\n" + "\n".join(misses + untagged) + "\n")
        r.report("lint", "FAIL" if misses else "PASS",
                 f"Correction quotes: {len(misses)} missing", misses)
        if untagged:
            r.report("lint", "LOOK", f"{len(untagged)} README line(s) with untagged "
                     "numbers: each must be tagged, quoted or an identifier", untagged)

    if "outputs" in steps:
        found = summary_of_outputs(page)
        r.report("outputs", "PASS" if len(found) == 2 else "LOOK",
                 "from the executed index.ipynb" if found else
                 "no claims.report() or break-table line in the notebook outputs", found)
        check_rows(r, page / "agreement.json")

    print(f"{'FAIL' if r.failed else 'PASS'} overall in {time.time() - t0:.0f} s; "
          f"full output in {log}")
    return 1 if r.failed else 0


if __name__ == "__main__":
    sys.exit(main())
