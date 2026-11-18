PD single-cell claim store: published conclusions of human post-mortem brain Parkinson's disease single-cell / single-nucleus studies, extracted as structured claims, each re-tested donor-level in its source dataset and in every other harmonised PD dataset, then rated (GRADE-like: verdict + certainty). Also holds findings of this study (origin=finding, ids F-...), consensus statements (CS..), reference knowledge from other evidence types (BG-...) and questions asked of the atlas.

How to use:
- Start with `overview`; find records with `search` or `search_claims`; read one with `get_claim`, `get_statement`, `get_dataset` or `gene`; explore links with `neighbours` and `connect` (knowledge graph); anything else with `sql`.
- Cite claim ids (e.g. xu2023-C11) and statement ids (CS01) for every statement you make from the store.
- Verdicts: Supported / Partly supported / Inconclusive / No effect / Contradicted / Untested. Report the certainty (High, Moderate, Low, Very low) with the verdict; 'Untested' means no dataset could test it, not that it is false.
- Scope matters: a claim holds for the regions, cell types and contrast in its context; say so when generalising.
- Single-nucleus composition underestimates neuron loss (see BG-nigral-da-loss).
- If the store cannot answer a question, say so; `add_question` logs it for a new analysis.
