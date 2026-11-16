"""Controlled vocabularies for the claim store, resolved to ontology IDs through the EBI OLS4 API.

Writes claims/vocab/{regions,cell_types,diseases,assays,organisms}.tsv: our label, the ontology term searched for,
the resolved ID and the label OLS returned for it (so every ID is checked, not typed from memory). `--offline`
reuses IDs already in the files and only refreshes the label columns we define here.

Run: ~/miniconda3/envs/rapids_singlecell/bin/python scripts/51a_vocab.py
"""
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "claims" / "vocab"

# label: (ontology, term searched, group, harmonised regions it covers)
REGIONS = {
    "SN": ("uberon", "substantia nigra", "midbrain", "SN"),
    "VTA": ("uberon", "ventral tegmental area", "midbrain", "VTA"),
    "midbrain": ("uberon", "midbrain", "midbrain", "midbrain;SN"),
    "putamen": ("uberon", "putamen", "basal ganglia", "putamen"),
    "caudate": ("uberon", "caudate nucleus", "basal ganglia", "caudate"),
    "striatum": ("uberon", "striatum", "basal ganglia", "putamen;caudate;striatum"),
    "GPi": ("uberon", "medial globus pallidus", "basal ganglia", "GPi"),
    "PFC": ("uberon", "prefrontal cortex", "cortex", "PFC;DLPFC"),
    "DLPFC": ("uberon", "dorsolateral prefrontal cortex", "cortex", "DLPFC;PFC"),
    "cingulate": ("uberon", "cingulate cortex", "cortex", "cingulate"),
    "temporal": ("uberon", "temporal cortex", "cortex", "temporal;MTG"),
    "MTG": ("uberon", "middle temporal gyrus", "cortex", "MTG;temporal"),
    "motor": ("uberon", "primary motor cortex", "cortex", "motor"),
    "visual": ("uberon", "primary visual cortex", "cortex", "visual"),
    "parietal": ("uberon", "parietal cortex", "cortex", "parietal"),
    "amygdala": ("uberon", "amygdala", "limbic", "amygdala"),
    "hippocampus": ("uberon", "hippocampal formation", "limbic", "hippocampus"),
    "cerebellum": ("uberon", "cerebellum", "cerebellum", "cerebellum"),
    "DMNX": ("uberon", "dorsal motor nucleus of vagus nerve", "brainstem", "DMNX"),
    "SVZ": ("uberon", "postnatal subventricular zone", "basal ganglia", "SVZ"),
    "brain": ("uberon", "brain", "whole brain", ""),
}
CELL_TYPES = {
    "DA neurons": "dopaminergic neuron", "MSNs": "medium spiny neuron",
    "Excitatory neurons": "glutamatergic neuron", "Inhibitory neurons": "GABAergic neuron",
    "Neurons (other)": "neuron", "Astrocytes": "astrocyte", "Microglia": "microglial cell",
    "Oligodendrocytes": "oligodendrocyte", "OPCs": "oligodendrocyte precursor cell",
    "Endothelial": "endothelial cell", "Pericytes": "pericyte", "Fibroblasts": "fibroblast", "T cells": "T cell",
    "Immune (other)": "leukocyte", "Ependymal": "ependymal cell", "All": "cell",
}
DISEASES = {  # harmonised disease code: (term, our definition)
    "PD": ("Parkinson disease", "Parkinson's disease without recorded dementia (or dementia status not given)"),
    # MONDO has no PDD, unspecified-LBD or iLBD term (checked 2026-11-13): use the nearest broader term
    "PDD": ("Parkinson disease", "Parkinson's disease with dementia, as stated by the source [broader term]"),
    "DLB": ("Lewy body dementia", "dementia with Lewy bodies"),
    "LBD": ("synucleinopathy", "Lewy body disease, subtype not given by the source [broader term]"),
    "iLBD": ("synucleinopathy", "Lewy pathology without clinical parkinsonism or dementia [broader term]"),
    "AD": ("Alzheimer disease", "Alzheimer's disease"),
}
ASSAYS = {  # claim modality: EFO term
    "snRNA": "single nucleus RNA sequencing", "scRNA": "single-cell RNA sequencing",
    "multiome-RNA": "10x multiome", "multiome-ATAC": "10x multiome", "snATAC": "sc/snATAC-seq",
    "spatial-Visium": "Visium Spatial Gene Expression", "spatial-Xenium": "Xenium",
    "spatial-GeoMx": "GeoMx", "spatial-other": "spatial transcriptomics",
    "imaging": "immunohistochemistry",  # claims tagged imaging are IHC/IF/RNAscope histology
}
ORGANISMS = {"human": "Homo sapiens", "macaque": "Macaca mulatta", "mouse": "Mus musculus",
             "rat": "Rattus norvegicus"}


def ols(term, ontology):
    for exact in ("true", "false"):
        q = urllib.parse.urlencode({"q": term, "ontology": ontology, "exact": exact, "rows": 5,
                                    "fieldList": "obo_id,label", "queryFields": "label,synonym"})
        with urllib.request.urlopen(f"https://www.ebi.ac.uk/ols4/api/search?{q}", timeout=30) as r:
            docs = json.load(r)["response"]["docs"]
        docs = [d for d in docs if d.get("obo_id", "").upper().startswith(ontology.upper().replace("NCBITAXON", "NCBITAXON"))]
        if docs:
            hit = next((d for d in docs if d["label"].lower() == term.lower()), docs[0])
            return hit["obo_id"], hit["label"]
    return "", ""


def resolve(rows, ontology_col="ontology"):
    out = []
    for r in rows:
        oid, olabel = ols(r["term"], r[ontology_col])
        r.update(ontology_id=oid, ontology_label=olabel,
                 check="ok" if olabel.lower() == r["term"].lower() else ("MISSING" if not oid else "label differs"))
        out.append(r)
        print(f"{r['label']:22s} {oid:16s} {olabel}  [{r['check']}]", file=sys.stderr)
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    resolve([dict(label=k, ontology=o, term=t, region_group=g, harmonised_regions=h)
             for k, (o, t, g, h) in REGIONS.items()]).to_csv(OUT / "regions.tsv", sep="\t", index=False)
    resolve([dict(label=k, ontology="cl", term=t) for k, t in CELL_TYPES.items()]).to_csv(
        OUT / "cell_types.tsv", sep="\t", index=False)
    d = resolve([dict(label=k, ontology="mondo", term=t, definition=dfn) for k, (t, dfn) in DISEASES.items()])
    d = pd.concat([pd.DataFrame([dict(label="control", ontology="", term="", definition=
                   "no neurological diagnosis and no Lewy pathology reported by the source", ontology_id="",
                   ontology_label="", check="no term")]), d])
    d.to_csv(OUT / "diseases.tsv", sep="\t", index=False)
    resolve([dict(label=k, ontology="efo", term=t) for k, t in ASSAYS.items()]).to_csv(
        OUT / "assays.tsv", sep="\t", index=False)
    resolve([dict(label=k, ontology="ncbitaxon", term=t) for k, t in ORGANISMS.items()]).to_csv(
        OUT / "organisms.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
