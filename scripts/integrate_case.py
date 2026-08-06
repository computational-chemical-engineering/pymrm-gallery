"""Integrate a verified staged case into the gallery, end to end.

Usage: python scripts/integrate_case.py <CASE_ID> <BEFORE_ID> [--dry-run]

Does, in order (each step verified before the next):
  1. sanity: queue_cases/<ID>/page exists with meta.yaml + models_entry.yaml;
     a review/verification.md is present (warn loudly if not — the standing
     rule is a verifier on every ready page);
  2. splice models_entry into models.yaml before <BEFORE_ID> via
     splice_entry.splice (which handles the three entry shapes and flips
     ready/in-progress to published) — after checking the case does not
     already exist there (upgrade-in-place entries must be handled by hand);
  3. flip page/meta.yaml status to published; flip queue_cases/<ID>.yaml to
     published with page: pages/<ID>-<slug>/ (slug read from models.yaml);
  4. copy queue_cases/<ID>/page -> pages/<ID>-<slug>/;
  5. gates: check_metadata.py, run_pages.py --changed, check_agreement.py.

Git is deliberately NOT touched: the coordinator pulls, stages named paths,
checks the staged diff for deletions, commits and pushes. See
docs/playbooks/coordinator.md.
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from splice_entry import splice  # noqa: E402


def die(msg):
    sys.exit(f"integrate_case: {msg}")


def run_gate(args, timeout=900):
    print(f"\n$ {' '.join(args)}")
    r = subprocess.run(args, cwd=ROOT, capture_output=True, text=True,
                       timeout=timeout)
    tail = "\n".join((r.stdout + r.stderr).strip().splitlines()[-6:])
    print(tail)
    return r.returncode, tail


def main():
    if len(sys.argv) < 3:
        die(__doc__)
    case, before = sys.argv[1], sys.argv[2]
    dry = "--dry-run" in sys.argv

    qdir = ROOT / "queue_cases" / case
    page = qdir / "page"
    meta_p = page / "meta.yaml"
    entry_p = qdir / "models_entry.yaml"
    for p in (page, meta_p, entry_p):
        if not p.exists():
            die(f"missing {p.relative_to(ROOT)}")
    ver = qdir / "review" / "verification.md"
    if not ver.exists():
        print(f"WARNING: no {ver.relative_to(ROOT)} — the standing rule is a "
              "verifier on every ready page. Proceeding only if you know why.")

    models = (ROOT / "models.yaml").read_text()
    existing = len(re.findall(rf"^  - id: {re.escape(case)}$", models, re.M))
    if existing:
        die(f"models.yaml already has {existing} entry for {case} — "
            "upgrade-in-place needed; do it by hand (see C2.10 incident).")

    if dry:
        print(f"[dry-run] would splice {case} before {before} and integrate")
        return

    splice(case, before)

    d = yaml.safe_load((ROOT / "models.yaml").read_text())
    m = [e for e in d["models"] if e["id"] == case]
    if len(m) != 1 or m[0].get("status") != "published":
        die(f"post-splice check failed: {len(m)} entries, "
            f"status {m[0].get('status') if m else None}")
    slug = m[0]["slug"]
    dest_rel = f"pages/{case}-{slug}/"
    if m[0].get("page") not in (dest_rel, dest_rel.rstrip("/")):
        die(f"models.yaml page '{m[0].get('page')}' != '{dest_rel}' — fix the "
            "entry's page: field first (directory name must equal the slug).")

    t = meta_p.read_text()
    t2 = re.sub(r"^status: (ready|draft|planned|in-progress).*$",
                "status: published", t, count=1, flags=re.M)
    meta_p.write_text(t2)
    if "status: published" not in t2.splitlines()[0:20].__str__() and \
       "status: published" not in t2:
        die("could not flip page/meta.yaml status")

    cy = ROOT / "queue_cases" / f"{case}.yaml"
    t = cy.read_text()
    if "\nstatus: ready\n" in t:
        t = t.replace("\nstatus: ready\n",
                      f"\nstatus: published\npage: {dest_rel}\n", 1)
        cy.write_text(t)
    elif "status: published" in t:
        print("case yaml already published")
    else:
        die(f"{cy.name}: no 'status: ready' line to flip — inspect by hand")

    dest = ROOT / dest_rel.rstrip("/")
    if dest.exists():
        die(f"{dest_rel} already exists — refusing to overwrite")
    shutil.copytree(page, dest)
    print(f"copied page -> {dest_rel}")

    py = sys.executable
    rc, _ = run_gate([py, "scripts/check_metadata.py"])
    if rc:
        die("check_metadata failed — fix before proceeding")
    rc, _ = run_gate([py, "scripts/run_pages.py", "--changed"])
    if rc:
        die("run_pages --changed failed")
    rc, out = run_gate([py, "scripts/check_agreement.py"], timeout=300)
    if rc or " 0 regression" not in out:
        print("check_agreement did not report 0 regressions — review above; "
              "an intended correction is legitimate, silence is not.")

    print(f"\n{case} integrated as {dest_rel}. Now: git pull --rebase, stage "
          "named paths (models.yaml, the case yaml, the page dir), check the "
          "staged diff for deletions, commit, push.")


if __name__ == "__main__":
    main()
