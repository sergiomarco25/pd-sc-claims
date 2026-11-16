"""Question loop, analysis side: run every question with status `specified` and enter its result as a finding.

Usage: python -I scripts/55_answer_questions.py [question_id ...]
For each question in claims/questions.yaml whose status is `specified` (or the ids given):
  1. run its spec with scripts/54_question_analysis.py -> results/questions/<qid>.tsv
  2. add_finding() with the question's finding draft and an evidence reference to that table
  3. set the question's status to `tested` and record the finding id
Then rebuild and rate: 51b -> 52 -> 51b -> 51c. The atlas's question log is updated by Claude (ArtifactData), see
.claude/skills/atlas-questions/SKILL.md.
"""
import importlib.util
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.store import QUESTIONS, add_finding, _dump  # noqa: E402

spec_mod = importlib.util.spec_from_file_location("qa", ROOT / "scripts" / "54_question_analysis.py")
qa = importlib.util.module_from_spec(spec_mod)
spec_mod.loader.exec_module(qa)


def main():
    want = set(sys.argv[1:])
    qs = yaml.safe_load(QUESTIONS.read_text()) or []
    for q in qs:
        if (want and q["id"] not in want) or (not want and q.get("status") != "specified"):
            continue
        spec, draft = q["proposed"]["spec"], dict(q["proposed"]["finding"])
        res = qa.run(spec)
        qa.OUT.mkdir(parents=True, exist_ok=True)
        out = qa.OUT / f"{q['id']}.tsv"
        res.to_csv(out, sep="\t", index=False)
        if res.empty:
            q["status"], q["answer"] = "no_data", "No dataset had enough donors, genes or groups for this analysis."
            print(q["id"], "no data")
            continue
        corr = "correlate" in spec
        evidence = {"table": f"questions/{q['id']}.tsv", "scale": "fisher_z" if corr else "smd",
                    "columns": {"n": "n"} if corr else {}}
        draft.setdefault("themes", "")
        fid = draft.pop("id")
        add_finding(id=fid, evidence=evidence,
                    analysis={"script": "scripts/54_question_analysis.py", "results": f"results/questions/{q['id']}.tsv",
                              "question": q["id"]},
                    notes=f"From question {q['id']}: {q['text']}", **draft)
        q["status"] = "tested"
        q["finding_ids"] = sorted(set(q.get("finding_ids", []) + [fid]))
        print(q["id"], "->", fid, f"({res.dataset.nunique()} datasets, {len(res)} units)")
    _dump(QUESTIONS, qs)


if __name__ == "__main__":
    main()
