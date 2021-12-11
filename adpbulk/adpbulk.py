"""
Class for Pseudobulking
"""
from typing import Union, List
import itertools as it
import numpy as np
import pandas as pd
import anndata as ad
from tqdm import tqdm


class ADPBulk:
    def __init__(
            self,
            adat: ad.AnnData,
            groupby: Union[List[str], str],
            method: str = "sum",
            name_delim: str = "-",
            group_delim: str = ".",
            use_raw: bool = False):
        """
        Class of Pseudo-Bulking `AnnData` objects based on categorical variables
        found in the `.obs` attribute

        inputs:
            adat: anndata.AnnData
                The `AnnData` object to process
            groupby: Union[List[str], str]
                The categories to group by. Can provide as a single value
                or a list of values.
            method: str
                The method to aggregate with (sum[default], mean, median)
            name_delim: str
                The delimiter to use when grouping multiple categories together.
                example: 'cat1{delim}cat2'
            group_delim: str
                The delimiter to use for adding the value to its category.
                example: 'cat{delim}value'
            use_raw: bool
                Whether to use the `.raw` attribute on the `AnnData` object
        """

        self.agg_methods = {
            "sum": np.sum,
            "mean": np.mean,
            "median": np.median}

        self.adat = adat
        self.groupby = groupby
        self.method = method
        self.name_delim = name_delim
        self.group_delim = group_delim
        self.use_raw = use_raw

        self.group_idx = dict()
        self.groupings = list()
        self.grouping_masks = dict()

        self.meta = pd.DataFrame([])
        self.matrix = pd.DataFrame([])
        self.samples = pd.DataFrame([])

        self._isfit = False
        self._istransform = False

        self._validate()

    def _validate(self):
        """
        validates that the input is as expected
        """
        self._validate_anndata()
        self._validate_groups()
        self._validate_method()
        self._validate_raw()

    def _validate_anndata(self):
        """
        validates that the anndata object is as expected
        """
        if self.adat.X is None:
            raise AttributeError("Provided Matrix is None")
        if self.adat.obs is None:
            raise AttributeError("Provided Obs are None")
        if self.adat.var is None:
            raise AttributeError("Provided Var are None")

    def _validate_groups(self):
        """
        validates that the groups are as expected
        """
        # convert groups to list if provided as str
        if isinstance(self.groupby, str):
            self.groupby = [self.groupby]

        if isinstance(self.groupby, list):
            self.groupby = np.unique(self.groupby)
            for group in self.groupby:
                self._validate_group(group)
        else:
            raise TypeError("Groupby is not a list or str")

    def _validate_group(self, group):
        """
        confirms that provided group is a column in the observations
        """
        if group not in self.adat.obs.columns:
            raise ValueError(f"Provided group {group} not in observations")

    def _validate_method(self):
        """
        confirms that the method is known
        """
        if self.method not in self.agg_methods.keys():
            raise ValueError(
                f"Provided method {self.method} not in known methods {''.join(self.agg_methods)}")

    def _validate_raw(self):
        """
        if the `use_raw` flag is provided will confirm that
        the raw field is present
        """
        if self.use_raw and self.adat.raw is None:
            raise AttributeError(
                "use_raw provided, but no raw field is found in AnnData")


    def _fit_indices(self):
        """
        determines the indices for each of the provided groups
        """
        for group in self.groupby:
            unique_values = np.unique(self.adat.obs[group].values)
            self.group_idx[group] = {
                uv: set(np.flatnonzero(self.adat.obs[group].values == uv))
                    for uv in tqdm(unique_values, desc=f"fitting indices: {group}")}

    def _get_mask(self, pairs: tuple) -> np.ndarray:
        """
        retrieve the indices for the provided values from their respective groups
        calculates the global intersection between the sets
        """
        group_indices = []
        for j, key in enumerate(pairs):
            group_indices.append(self.group_idx[self.groupby[j]][key])
        mask = set.intersection(*group_indices)
        return np.array(list(mask))

    def _get_name(self, pairs: tuple) -> str:
        """
        create a name for the provided values based on their respective groups
        """
        name = self.name_delim.join([
                f"{self.groupby[i]}{self.group_delim}{pairs[i]}" for i in range(len(pairs))])
        return name

    def _get_agg(self, mask: np.ndarray) -> np.ndarray:
        """
        runs the aggregation function with the provided sample mask
        """
        if self.use_raw:
            mat = self.adat.raw.X[mask]
        else:
            mat = self.adat.X[mask]
        return self.agg_methods[self.method](mat, axis=0)

    def _prepare_meta(self, pairs: tuple) -> dict:
        """
        defines the meta values for the pairs
        """
        values = {
            self.groupby[idx]: pairs[idx] for idx in np.arange(len(self.groupby))
            }
        values["SampleName"] = self._get_name(pairs)
        return values

    def _build_groupings(self):
        """
        generate the grouping iterable
        """
        if len(self.groupby) > 1:
            group_keys = [self.group_idx[g].keys() for g in self.groupby]
            self.groupings = list(it.product(*group_keys))
        else:
            self.groupings = self.group_idx[self.groupby[0]]

    def _build_masks(self):
        """
        generate the masks for each of the groupings and generates the
        metadata for each of the groupings
        """
        self.grouping_masks = dict()
        self.meta = []
        for pairs in self.groupings:
            if not isinstance(pairs, tuple):
                pairs = tuple([pairs])
            mask = self._get_mask(pairs)
            if mask.size > 0:
                self.grouping_masks[pairs] = mask
                self.meta.append(self._prepare_meta(pairs))

        if len(self.meta) == 0:
            raise ValueError("No combinations of the provided groupings found")
        self.meta = pd.DataFrame(self.meta)

    def fit(self):
        """
        fits the indices for each of the groups
        """
        self._fit_indices()
        self._build_groupings()
        self._build_masks()
        self._isfit = True

    def transform(self) -> pd.DataFrame:
        """
        performs the aggregation based on the fit indices
        """
        if not self._isfit:
            raise AttributeError("Please fit the object first")

        matrix = []
        for pairs in tqdm(self.groupings, desc="Aggregating Samples"):
            if not isinstance(pairs, tuple):
                pairs = tuple([pairs])
            if pairs in self.grouping_masks:
                matrix.append(self._get_agg(self.grouping_masks[pairs]))

        self.matrix = pd.DataFrame(
            matrix,
            index=self.meta.SampleName.values,
            columns=self.adat.var.index.values)

        self._istransform = True
        return self.matrix

    def fit_transform(self):
        """
        firs the indices and performs the aggregation based on those indices
        """
        self.fit()
        return self.transform()

    def get_meta(self) -> pd.DataFrame:
        """
        return the meta dataframe
        """
        if not self._isfit:
            raise AttributeError("Please fit the object first")
        return self.meta
