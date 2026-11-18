"""claims/RATING.md states the rules the engine's rating (src/claimstore/rating.py) applies: same version, same thresholds (step 2.4).
Standard library only."""
import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _script_constants():
    tree = ast.parse((ROOT / "src" / "claimstore" / "rating.py").read_text())
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "RATING_VERSION":
            out["version"] = ast.literal_eval(node.value)
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "THRESH":
            out["thresh"] = {k.arg: ast.literal_eval(k.value) if not isinstance(k.value, ast.BinOp)
                             else eval(compile(ast.Expression(k.value), "", "eval")) for k in node.value.keywords}
    return out


def test_rating_doc_matches_script():
    doc = (ROOT / "claims" / "RATING.md").read_text()
    c = _script_constants()
    assert re.search(rf"^# .*version {re.escape(c['version'])}\s*$", doc, re.M), "RATING.md title version != script"
    table = doc.split("## 6. Thresholds", 1)[1]
    for name, value in c["thresh"].items():
        row = next((ln for ln in table.splitlines() if re.search(rf"\b{name}\b", ln)), None)
        assert row, f"threshold {name} missing from RATING.md section 6"
        shown = "2/3" if abs(value - 2 / 3) < 1e-9 else f"{value:g}"
        assert shown in row, f"{name}: script {value}, RATING.md row: {row}"


if __name__ == "__main__":
    test_rating_doc_matches_script()
    print("PASS")
