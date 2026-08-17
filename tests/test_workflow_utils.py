import sys
from pathlib import Path

import anndata as ad
import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts"))

from workflow_utils import get_counts, record_provenance


def test_get_counts_accepts_sparse_integer_counts():
    adata = ad.AnnData(sparse.csr_matrix([[0, 1], [2, 0]]))
    adata.layers["counts"] = adata.X.copy()
    assert get_counts(adata, "counts").shape == (2, 2)


def test_get_counts_rejects_normalized_values():
    adata = ad.AnnData(np.array([[0.0, 0.25], [1.5, 0.0]]))
    try:
        get_counts(adata, "counts")
    except ValueError as exc:
        assert "not integer-like" in str(exc)
    else:
        raise AssertionError("Normalized values were accepted as raw counts")


def test_provenance_survives_h5ad_roundtrip(tmp_path):
    adata = ad.AnnData(sparse.csr_matrix([[1, 0], [0, 1]]))
    record_provenance(adata, "unit_test", {"seed": 7}, ["anndata", "missing-package"])
    output = tmp_path / "roundtrip.h5ad"
    adata.write_h5ad(output)
    restored = ad.read_h5ad(output)
    assert len(restored.uns["scrna_workflow"]["runs"]) == 1
    assert "unit_test" in restored.uns["scrna_workflow"]["runs"][0]

