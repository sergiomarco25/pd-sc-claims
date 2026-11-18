# Claim format

One YAML file per paper, `claims/<paper_id>.yaml`, holding a list of claims. Each claim is one conclusion the authors
draw from their data, written so that it can be re-tested in other datasets. Fields, their meaning and their allowed
values are in `config/schema.toml`; the generated JSON Schema (`claims/schema/claim.schema.json`) checks every file in
CI. Ids are `<paper_id>-C<nn>`, numbered in order of appearance; findings of this project go to `claims/findings.yaml`.
