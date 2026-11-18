{{title}}: published conclusions of studies in this field, extracted as structured claims, each re-tested in its source dataset and in every other harmonised dataset, then rated (GRADE-like: verdict + certainty). Also holds findings of this project (origin=finding, ids F-...), consensus statements (CS..), reference knowledge (BG-...) and questions.

How to use:
- Start with `overview`; find records with `search` or `search_claims`; read one with `get_claim`, `get_statement`, `get_dataset` or `gene`; explore links with `neighbours` and `connect` (knowledge graph); anything else with `sql`.
- Cite claim ids and statement ids for every statement you make from the store.
- Verdicts: Supported / Partly supported / Inconclusive / No effect / Contradicted / Untested. Report the certainty (High, Moderate, Low, Very low) with the verdict; 'Untested' means no dataset could test it, not that it is false.
- Scope matters: a claim holds for the regions, cell types and contrast in its context; say so when generalising.
- If the store cannot answer a question, say so; `add_question` logs it for a new analysis.
