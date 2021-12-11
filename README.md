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
The simplest use case is to aggregate on a single category. This will aggregate all the observations belonging to the same class within the category and return a pseudo-bulked matrix with dimensions equal to the number of values within the category. 
```python3
from adpbulk import ADPBulk

# initialize the object
adpb = ADPBulk(adat, "category_name")

# perform the pseudobulking
pseudobulk_matrix = adpb.fit_transform()

# retrieve the sample meta data (useful for easy incorporation with edgeR)
sample_meta = adpb.get_meta()
```

## Multiple Category Pseudo-Bulk
A common use case is to aggregate on multiple categories. This will aggregate all observations beloging to the combination of classes within two categories and return a pseudo-bulked matrix with dimensions equal to the number of values of nonzero intersections between categories. 
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
It may also be useful to aggregate using an alternative function besides the sum - this option will allow you to choose between sum, mean, and median as an aggregation function.
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
Here is a function to generate an AnnData object to test the module or to play with the object if unfamiliar.
```python3
import numpy as np
import pandas as pd
import anndata as ad

def build_adat():
    """
    creates an anndata for testing
    """
    # generates random values (mock transformed data)
	mat = np.random.random((SIZE_N, SIZE_M))

	# generates random values (mock raw count data)
    raw = np.random.randint(0, 1000, (SIZE_N, SIZE_M))

	# creates the observations and categories
    obs = pd.DataFrame({
        "cell": [f"b{idx}" for idx in np.arange(SIZE_N)],
        "cA": np.random.choice(np.random.choice(5)+1, SIZE_N),
        "cB": np.random.choice(np.random.choice(5)+1, SIZE_N),
        "cC": np.random.choice(np.random.choice(5)+1, SIZE_N),
        "cD": np.random.choice(np.random.choice(5)+1, SIZE_N),
        }).set_index("cell")

	# creates the variables (genes) and categories
    var = pd.DataFrame({
        "symbol": [f"g{idx}" for idx in np.arange(SIZE_M)],
        "cA": np.random.choice(np.random.choice(5)+1, SIZE_M),
        "cB": np.random.choice(np.random.choice(5)+1, SIZE_M),
        "cC": np.random.choice(np.random.choice(5)+1, SIZE_M),
        "cD": np.random.choice(np.random.choice(5)+1, SIZE_M),
        }).set_index("symbol")
    
	# Creates the `AnnData` object
	adat = ad.AnnData(
            X=mat,
            obs=obs,
            var=var)
    
	# Creates an `AnnData` object to simulate the `.raw` attribute
	adat_raw = ad.AnnData(
            X=raw,
            obs=obs,
            var=var)
    
	# Sets the `.raw` attribute
	adat.raw = adat_raw
    
	return adat

adat = build_adat()
```
