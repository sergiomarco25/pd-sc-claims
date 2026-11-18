"""claimstore: the field-independent engine behind a re-tested claim store (claims/COMMUNITY_PLAN.md step 2.1).

Modules: config (project file), meta (pooling), contrasts (case and reference groups from config/contrasts.toml), store
(findings and their evidence), rating (verdict and certainty), kg (knowledge graph), release (public release and its
redaction), public_findings, release_diff, linked_data, audit, mcp_server and mcp_http (serving the store).
Everything about one field (its datasets, papers, cell types, regions, contrasts, wording) comes from the project's
config files and claim files; the engine names none of it (tests/test_engine_boundary.py). A project adds its own
package for what is specific to it: dataset loaders, harmonisation, plotting and its analysis scripts.
"""
