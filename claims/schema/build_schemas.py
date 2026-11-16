"""Write the machine-checkable claim schemas (COMMUNITY_PLAN step 2.3) from one field list.

Outputs JSON Schema (draft-07) files next to this script:
  claim.schema.json     one published claim in claims/<paper_id>.yaml (claims/SCHEMA.md describes each field)
  finding.schema.json   one finding of this study in claims/findings.yaml (same fields plus origin and evidence)

Validated in CI by tests/test_claim_schema.py, which also checks that the vocabularies match the ones
scripts/05_merge_claims.py enforces. Change a field here, run this script, commit both files.
Usage: python3 -I claims/schema/build_schemas.py      (standard library only)
"""
import json
import re
import sys
from pathlib import Path

SCHEMA_VERSION = "1.0.0"
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from pdrep import config as cfg  # noqa: E402

BASE = cfg.PUBLIC_REPOSITORY.rstrip("/") + "/blob/main/claims/schema/"  # the schemas' public home

VOCAB = {
    "centrality": ["headline", "main", "secondary"],
    "modality": ["snRNA", "scRNA", "multiome-RNA", "snATAC", "multiome-ATAC", "spatial-Visium", "spatial-Xenium",
                 "spatial-GeoMx", "spatial-other", "imaging"],
    "region": ["SN", "VTA", "midbrain", "putamen", "caudate", "striatum", "PFC", "DLPFC", "cingulate", "temporal", "MTG",
               "motor", "visual", "parietal", "amygdala", "hippocampus", "cerebellum", "DMNX", "GPi", "multiple"],
    "cell_type": ["DA neurons", "Neurons (other)", "Excitatory neurons", "Inhibitory neurons", "MSNs", "Astrocytes",
                  "Microglia", "Oligodendrocytes", "OPCs", "Endothelial", "Pericytes", "Fibroblasts", "T cells",
                  "Immune (other)", "All"],
    "claim_type": ["composition", "de_gene", "de_pathway", "state_existence", "state_enrichment", "genetic_enrichment",
                   "communication", "pathology_association", "trajectory", "chromatin", "spatial_colocalisation",
                   "clinical_association", "other"],
    "direction": ["up", "down", "enriched", "depleted", "exists", "associated_pos", "associated_neg", "NA"],
    "testable_with": ["rna", "atac", "spatial", "metadata", "not_testable"],
    "test_template": ["composition_donor", "de_pseudobulk", "geneset_shift", "state_signature", "genetic_enrichment",
                      "communication_liana", "pathology_correlation", "trajectory", "chromatin_da",
                      "spatial_neighbourhood", "not_testable"],
}
LIST_FIELDS = {"modality", "region", "cell_type"}  # one or more values separated by ";" or ","
REQUIRED = ["id", "paper", "statement", "location", "centrality", "modality", "region", "cell_type", "claim_type",
            "direction", "comparison", "unit_of_analysis", "testable_with", "test_template"]
TEXT = {
    "id": "<paper_id>-C<nn>, numbered in order of appearance",
    "paper": "paper id (inventory/papers.tsv); the file name without .yaml",
    "statement": "the conclusion in one sentence, close to the authors' wording",
    "quote": "short verbatim quote supporting it (at most 40 words)",
    "location": "figure / panel / section",
    "centrality": "headline (abstract or title) / main (Results heading or main figure) / secondary",
    "modality": "assay(s) the claim rests on",
    "region": "brain region(s) the claim is about",
    "cell_type": "broad cell class(es)",
    "cell_state": "named subpopulation or signature, or empty",
    "claim_type": "kind of claim; decides the test",
    "direction": "direction asserted",
    "comparison": "groups compared, e.g. 'PD vs control', 'Braak 5-6 vs 3-4', 'within PD'",
    "genes": "genes needed to test the claim, as reported",
    "gene_source": "where the full gene list is, or 'text'",
    "reported_stat": "effect size / statistic / p or FDR as reported, and n",
    "unit_of_analysis": "what the authors' test treated as independent (cells, donors, samples, ROIs)",
    "method": "test used by the authors",
    "testable_with": "data type needed to re-test",
    "test_template": "re-test template (src/pdrep/tests/engine.py)",
    "requires": "metadata or features needed beyond disease labels",
    "notes": "caveats, ambiguities",
}


def list_pattern(values):
    alt = "|".join(re.escape(v) for v in values)
    return rf"^({alt})(\s*[;,]\s*({alt}))*$"


def field(name):
    p = {"description": TEXT[name]}
    if name == "genes":
        p.update(type="array", items={"type": "string", "minLength": 1})
    elif name in LIST_FIELDS:
        p.update(type="string", pattern=list_pattern(VOCAB[name]))
    elif name in VOCAB:
        p.update(type="string", enum=VOCAB[name])
    else:
        p.update(type="string")
    if name in REQUIRED and p.get("type") == "string" and "enum" not in p and "pattern" not in p:
        p["minLength"] = 1
    return p


def claim_schema():
    props = {k: field(k) for k in TEXT}
    props["id"]["pattern"] = r"^[a-z][a-z0-9]*-C[0-9]{2}$"
    return {
        "$schema": "http://json-schema.org/draft-07/schema#",
        "$id": BASE + "claim.schema.json",
        "title": "Published claim",
        "description": f"One conclusion of a published study, as extracted for re-testing (schema {SCHEMA_VERSION}). "
                       "Field meanings: claims/SCHEMA.md.",
        "version": SCHEMA_VERSION,
        "type": "object",
        "required": REQUIRED,
        "properties": props,
        "additionalProperties": False,
    }


def finding_schema():
    s = claim_schema()
    props = s["properties"]
    props["id"] = {"description": "F-<short-name>", "type": "string", "pattern": r"^F-[A-Za-z0-9][A-Za-z0-9-]*$"}
    props["paper"] = {"description": "always this-study", "type": "string", "enum": ["this-study"]}
    props["test_template"] = {"description": "'finding': estimates come from the finding's own analysis (evidence), "
                                            "not from a re-test template",
                              "type": "string", "enum": VOCAB["test_template"] + ["finding"]}
    props.update({
        "origin": {"description": "record origin", "type": "string", "enum": ["finding"]},
        "themes": {"description": "biological themes, ';'-separated", "type": "string"},
        "contrast": {"description": "contrast code (src/pdrep/contrasts.py)", "type": "string"},
        "null_claim": {"description": "true if the finding asserts the absence of an effect", "type": "boolean"},
        "evidence": {"description": "where the per-dataset estimates are (pdrep.store.resolve_evidence)",
                     "type": "object", "required": ["table"],
                     "properties": {"table": {"type": "string", "minLength": 1}}},
        "analysis": {"description": "how the finding was produced (script, spec, date)", "type": "object"},
    })
    s.update({"$id": BASE + "finding.schema.json", "title": "Finding of this study",
              "description": f"A result of this study's own analyses, in the claim format (schema {SCHEMA_VERSION}).",
              "required": REQUIRED + ["origin", "evidence"]})
    return s


if __name__ == "__main__":
    for name, s in (("claim", claim_schema()), ("finding", finding_schema())):
        out = HERE / f"{name}.schema.json"
        out.write_text(json.dumps(s, indent=2, ensure_ascii=False) + "\n")
        print(out)
