# adpbulk

# Summary
Performs pseudobulking of an AnnData object based on columns available in the `.obs` dataframe. This was originally intended to be used to pseudo-bulk single-cell RNA-seq data to higher order combinations of the data as to use existing RNA-seq differential expression tools such as edgeR and DESeq2. An example usage of this would be pseudobulking cells based on their cluster, sample of origin, or CRISPRi guide identity. This is intended to work on both individual categories (i.e. one of the examples) or combinations of categories (two of the three, etc.)

# Installation
```bash
git clone https://github.com/noamteyssier/adpbulk
cd adpbulk
pip install .
pytest -v 
```

# Usage
This package is intended to be used as a python module. 

## Single Category Pseudo-Bulk
```python3
from adpbulk import ADPBulk

# initialize the object
adpb = ADPBulk(adat, "category")

# perform the pseudobulking
pseudobulk_matrix = adpb.fit_transform()

# retrieve the sample meta data (useful for easy incorporation with edgeR)
sample_meta = adpb.get_meta()
```

## Multiple Category Pseudo-Bulk
```python3
from adpbulk import ADPBulk

# initialize the object
adpb = ADPBulk(adat, ["category_a", "category_b"])

# perform the pseudobulking
pseudobulk_matrix = adpb.fit_transform()

# retrieve the sample meta data (useful for easy incorporation with edgeR)
sample_meta = adpb.get_meta()
```

## Pseudo-Bulk using raw counts
Some differential expression software expects the counts to be untransformed counts. SCANPY uses the `.raw` attribute in its AnnData objects to store the initial AnnData object before transformation. If you'd like to perform the pseudo-bulk aggregation using these raw counts you can provide the `use_raw=True` flag. 
```python3
from adpbulk import ADPBulk

# initialize the object w. aggregation on the `.raw` attribute
adpb = ADPBulk(adat, ["category_a", "category_b"], use_raw=True)

# perform the pseudobulking
pseudobulk_matrix = adpb.fit_transform()

# retrieve the sample meta data (useful for easy incorporation with edgeR)
sample_meta = adpb.get_meta()
```

## Alternative Aggregation Options
```python3
from adpbulk import ADPBulk

# initialize the object w. an alternative aggregation option
# aggregation options are: sum, mean, and median
# default aggregation is sum
adpb = ADPBulk(adat, "category", method="mean")

# perform the pseudobulking
pseudobulk_matrix = adpb.fit_transform()

# retrieve the sample meta data (useful for easy incorporation with edgeR)
sample_meta = adpb.get_meta()
```

## Alternative Formatting Options
```python3
from adpbulk import ADPBulk

# initialize the object w. alternative name formatting options
adpb = ADPBulk(adat, ["category_a", "category_b"], name_delim=".", group_delim="::")

# perform the pseudobulking
pseudobulk_matrix = adpb.fit_transform()

# retrieve the sample meta data (useful for easy incorporation with edgeR)
sample_meta = adpb.get_meta()
```


## Example AnnData Function
Here is a function to generate
```python3
import numpy as np
import pandas as pd
import anndata as ad

def build_adat():
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

adat = build_adat()
```
