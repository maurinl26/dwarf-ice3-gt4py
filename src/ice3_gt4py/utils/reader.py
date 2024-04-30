# -*- coding: utf-8 -*-
import numpy as np
import xarray as xr


class NetCDFReader:
    def __init__(self, filename: str):
        self.filename = filename

    def get_field(self, name: str):
        with xr.open_dataset(self.filename) as ds:
            return ds[name].data
