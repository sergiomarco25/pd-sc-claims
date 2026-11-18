"""Start a new claim-store project for another field with one command (claims/COMMUNITY_PLAN.md step 2.6).

Usage: python3 -I scripts/63_new_project.py DIR --name NAME --title TITLE [--description TEXT] [--package PKG]
                                            [--env-prefix PREFIX] [--repository URL] [--author NAME] [--git] [--force]

Writes DIR from template/ (placeholders {{name}}, {{title}}, ... filled in), with:
  config/      project.toml (paths, pipeline, release rules, facts about the field), schema.toml (claim format),
               contrasts.toml, harmonise.toml (data schema, cell classes, markers), server_instructions.md
  claims/      an example claim file, empty findings / reference knowledge / questions, SCHEMA.md, RATING.md (the
               engine's rating rules, copied), and the claim schemas built from the new schema.toml
  inventory/   papers.tsv and datasets.tsv with this project's columns and an example row
  src/         a copy of the engine (src/claimstore, with ENGINE.txt naming the source commit) and a package for the
               field (--package, default the name) with an example dataset loader
  tests/       the new project's checks (config, contrasts, schema, claim files), run by its CI workflow
  LICENSE, LICENSE-DATA.md, CITATION.cff, README.md
  .claude/skills/claim-replication/   the method step by step as a Claude Code skill (copied)
then runs the new project's checks once. With --git it also makes the first commit. Standard library only (the
checks need pyyaml, jsonschema and pandas; they are skipped with a note when missing).
"""
import argparse
import re
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "template"
COPIED = ["LICENSE", "LICENSE-DATA.md", "claims/RATING.md"]
SKILL = ".claude/skills/claim-replication"  # copied whole, like the engine


def git_commit():
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"], capture_output=True, text=True,
                              check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def fill(text, values):
    out = re.sub(r"\{\{(\w+)\}\}", lambda m: values.get(m.group(1), m.group(0)), text)
    left = re.findall(r"\{\{\w+\}\}", out)
    if left:
        raise ValueError(f"unfilled placeholders: {sorted(set(left))}")
    return out


def generate(out, values, force=False):
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        if not force:
            raise SystemExit(f"{out} exists and is not empty (use --force to replace it)")
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)
    for src in sorted(p for p in TEMPLATE.rglob("*") if p.is_file() and "__pycache__" not in p.parts):
        rel = fill(src.relative_to(TEMPLATE).as_posix(), values)
        dst = out / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(fill(src.read_text(), values))
    for rel in COPIED:
        (out / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / rel, out / rel)
    shutil.copytree(ROOT / SKILL, out / SKILL, ignore=shutil.ignore_patterns("__pycache__"))
    eng = out / "src" / "claimstore"
    eng.mkdir(parents=True, exist_ok=True)
    for f in sorted(p for ext in ("*.py", "*.R", "*.html") for p in (ROOT / "src" / "claimstore").glob(ext)):
        shutil.copy2(f, eng / f.name)
    (eng / "ENGINE.txt").write_text(f"claimstore engine copied from {values['source']} at commit "
                                    f"{values['engine_commit']} on {values['date']}.\nUpdate it from there; "
                                    "do not edit it in this project.\n")
    (out / "claims" / "schema").mkdir(parents=True, exist_ok=True)
    (out / "claims" / "schema" / "build_schemas.py").write_text(
        '"""Write claims/schema/*.json from config/schema.toml (python3 -I claims/schema/build_schemas.py)."""\n'
        "import sys\nfrom pathlib import Path\n\nHERE = Path(__file__).resolve().parent\n"
        'sys.path.insert(0, str(HERE.parents[1] / "src"))\n'
        "from claimstore.schema import *  # noqa: E402,F401,F403\nfrom claimstore.schema import write  # noqa: E402\n\n"
        'if __name__ == "__main__":\n    for p in write(HERE):\n        print(p)\n')
    # the claim schemas of the new project, built by its own engine copy from its own config
    subprocess.run([sys.executable, "-I", str(out / "claims" / "schema" / "build_schemas.py")], check=True,
                   capture_output=True, text=True)
    return out


def run_checks(out):
    missing = [m for m in ("pytest", "yaml", "jsonschema", "pandas") if subprocess.run(
        [sys.executable, "-c", f"import {m}"], capture_output=True).returncode]
    if missing:
        print(f"checks skipped (missing: {', '.join(missing)}); run `pytest tests` in {out}")
        return None
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "tests"], cwd=out,
                       capture_output=True, text=True)
    print(r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr[-500:])
    return r.returncode == 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("dir")
    ap.add_argument("--name", required=True, help="short id, e.g. ad-sc-claims (lower case, digits, '-')")
    ap.add_argument("--title", required=True, help="e.g. 'AD single-cell claim store'")
    ap.add_argument("--description", default="")
    ap.add_argument("--package", help="Python package for the field (default: from the name)")
    ap.add_argument("--env-prefix", help="environment variable prefix (default: from the name, upper case)")
    ap.add_argument("--repository", default="", help="public repository URL, if any")
    ap.add_argument("--author", default="")
    ap.add_argument("--git", action="store_true", help="git init and commit the new project")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    if not re.fullmatch(r"[a-z][a-z0-9-]*", a.name):
        sys.exit("--name: lower case letters, digits and '-', starting with a letter")
    package = a.package or a.name.replace("-", "_")
    if not re.fullmatch(r"[a-z_][a-z0-9_]*", package) or package == "claimstore":
        sys.exit("--package: a Python identifier other than claimstore")
    values = {"name": a.name, "title": a.title.replace('"', "'"),
              "description": (a.description or a.title).replace('"', "'"), "package": package,
              "env_prefix": a.env_prefix or re.sub(r"\W", "_", a.name).upper(), "repository": a.repository,
              "author": a.author or "", "date": str(date.today()), "engine_commit": git_commit(),
              "source": "the claim-store reference project"}
    out = generate(a.dir, values, a.force)
    print(f"{out}: new project {a.name!r} (package {package}, env prefix {values['env_prefix']})")
    ok = run_checks(out)
    if a.git:
        subprocess.run(["git", "init", "-q", "-b", "main"], cwd=out, check=True)
        subprocess.run(["git", "add", "-A"], cwd=out, check=True)
        subprocess.run(["git", "commit", "-q", "-m", f"New claim-store project from the template ({values['engine_commit']})"],
                       cwd=out, check=True)
        print("committed")
    if ok is False:
        sys.exit(1)


if __name__ == "__main__":
    main()
