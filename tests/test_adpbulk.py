"""
Testing Suite for adpbulk
"""

import numpy as np
import pandas as pd
import anndata as ad
import pytest
from adpbulk import ADPBulk

SIZE_N = 100
SIZE_M = 100

np.random.seed(42)

def build_adat() -> ad.AnnData:
    """
    creates an anndata for testing
    """
    mat = np.random.random((SIZE_N, SIZE_M))
    raw = np.random.randint(0, 1000, (SIZE_N, SIZE_M))
    obs = pd.DataFrame({
        "cell": [f"b{idx}" for idx in np.arange(SIZE_N)],
        "cA": np.random.choice(np.random.choice(5)+1, SIZE_N),
        "cB": np.random.choice(np.random.choice(5)+1, SIZE_N),
        "cC": np.random.choice(np.random.choice(5)+1, SIZE_N),
        "cD": np.random.choice(np.random.choice(5)+1, SIZE_N),
        }).set_index("cell")
    var = pd.DataFrame({
        "symbol": [f"g{idx}" for idx in np.arange(SIZE_M)],
        "cA": np.random.choice(np.random.choice(5)+1, SIZE_M),
        "cB": np.random.choice(np.random.choice(5)+1, SIZE_M),
        "cC": np.random.choice(np.random.choice(5)+1, SIZE_M),
        "cD": np.random.choice(np.random.choice(5)+1, SIZE_M),
        }).set_index("symbol")
    adat = ad.AnnData(
            X=mat,
            obs=obs,
            var=var)
    adat_raw = ad.AnnData(
            X=raw,
            obs=obs,
            var=var)
    adat.raw = adat_raw
    return adat


def test_init():
    """
    tests whether the ADPBulk object can be instantiated correctly
    """
    adat = build_adat()

    # tests singular group conditions
    for group in adat.obs.columns:
        _ = ADPBulk(adat, groupby=group)

    # tests multiple group conditions
    _ = ADPBulk(adat, groupby=["cA", "cD"])

    assert True


def test_init_raw():
    """
    tests whether the ADPBulk object can be instantiated correctly
    """
    adat = build_adat()

    # tests singular group conditions
    for group in adat.obs.columns:
        _ = ADPBulk(adat, groupby=group, use_raw=True)

    # tests multiple group conditions
    _ = ADPBulk(adat, groupby=["cA", "cD"], use_raw=True)

    assert True


def test_init_missing_group():
    """
    tests whether the ADPBulk object be init incorrectly
    with missing groups
    """
    adat = build_adat()
    with pytest.raises(ValueError):
        _ = ADPBulk(adat, groupby="foobar")
    with pytest.raises(ValueError):
        _ = ADPBulk(adat, groupby=["foobar", "barbaz"])
    with pytest.raises(ValueError):
        _ = ADPBulk(adat, groupby=["foobar", "cA"])
    assert True

def test_init_missing_method():
    """
    tests whether the ADPBulk object be init incorrectly
    with missing method
    """
    adat = build_adat()
    with pytest.raises(ValueError):
        _ = ADPBulk(adat, groupby="cA", method="not_a_method")
    assert True

def test_fit_indices():
    """
    tests whether the indices are properly fit
    """
    adat = build_adat()
    adpb = ADPBulk(adat, "cA")
    adpb.fit()

    assert "cA" in adpb.group_idx.keys()
    for g in adat.obs["cA"].unique():
        assert g in adpb.group_idx["cA"].keys()
        assert len(adpb.group_idx["cA"][g]) ==\
                np.sum(adat.obs["cA"] == g)

def test_transform_singular():
    """
    tests whether the transformation can be performed
    """
    adat = build_adat()
    for condition in adat.obs.columns.values:
        for use_raw in [True, False]:
            for method in ["mean", "median", "sum"]:
                adpb = ADPBulk(adat, condition, method=method, use_raw=use_raw)
                matrix = adpb.fit_transform()
                assert isinstance(matrix, pd.DataFrame)

def test_transform_combination():
    """
    tests whether the transformation can be performed
    """
    adat = build_adat()
    for use_raw in [True, False]:
        for method in ["mean", "median", "sum"]:
            adpb = ADPBulk(adat, ["cA", "cB"], method=method, use_raw=use_raw)
            matrix = adpb.fit_transform()
            assert isinstance(matrix, pd.DataFrame)

def test_order_of_operations():
    """
    test whether the order of observations are properly organized
    """
    adat = build_adat()
    adpb = ADPBulk(adat, "cA")
    with pytest.raises(AttributeError):
        adpb.transform()
    assert True
