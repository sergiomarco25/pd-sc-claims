"""Nalls 2019 sumstats (chr, bp, alleles) -> MAGMA p-value file with 1000G EUR rsIDs (SNP P N).

Usage: python -I scripts/16b_magma_pfile.py <GCST009325.tsv> <g1000_eur.bim> <out.pval>
Match on chr + position + unordered allele pair; N = N_cases + N_controls.
"""
import sys
import pandas as pd

ss, bim, out = sys.argv[1:4]
b = pd.read_csv(bim, sep="\t", header=None, names=["chr", "snp", "cm", "bp", "a1", "a2"], dtype={"chr": str})
b["key"] = b.chr + ":" + b.bp.astype(str) + ":" + b[["a1", "a2"]].apply(lambda r: "".join(sorted(r)), axis=1)
b = b.drop_duplicates("key").set_index("key").snp
s = pd.read_csv(ss, sep="\t", dtype={"chromosome": str})
s["key"] = (s.chromosome + ":" + s.base_pair_location.astype(str) + ":" +
            s[["effect_allele", "other_allele"]].apply(lambda r: "".join(sorted(r)), axis=1))
s["SNP"] = s.key.map(b)
m = s.dropna(subset=["SNP"]).drop_duplicates("SNP")
m.assign(P=m.p_value, N=m.N_cases + m.N_controls)[["SNP", "P", "N"]].to_csv(out, sep="\t", index=False)
print(f"{len(s)} variants, {len(m)} matched to reference ({len(m) / len(s):.1%})")
