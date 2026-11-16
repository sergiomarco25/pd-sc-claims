# ggseg brain atlas outlines (vendored)

SVG path outlines used by `scripts/51c_store_browser.py` to draw the atlas brain map.

- Source: python-ggseg 0.1 (https://github.com/ggseg/python-ggseg, PyPI `ggseg`), licence MIT (as declared in its setup.py).
- `dk/`: Desikan-Killiany cortical parcellation, left hemisphere, lateral and medial views (`lateral_left`, `medial_left` are the hemisphere silhouettes).
- `aseg/`: FreeSurfer subcortical segmentation in one coronal and one midsagittal slice (`Coronal`, `Sagittal` are the slice outlines).
- Files copied unchanged on 2026-11-13. The outlines derive from the FreeSurfer fsaverage parcellations (Desikan et al. 2006; Fischl et al. 2002).
