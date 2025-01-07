# -*- coding: utf-8 -*-
from __future__ import annotations

from datetime import datetime
from functools import partial
from typing import TYPE_CHECKING

from ifs_physics_common.framework.grid import I, J, K
from ifs_physics_common.framework.storage import allocate_data_array
import numpy as np

from ice3_gt4py.initialisation.utils import initialize_field
from ice3_gt4py.utils.reader import NetCDFReader
import logging

if TYPE_CHECKING:
    from typing import Literal, Tuple

    from ifs_physics_common.framework.config import GT4PyConfig
    from ifs_physics_common.framework.grid import ComputationalGrid, DimSymbol
    from ifs_physics_common.utils.typingx import (
        DataArray,
        DataArrayDict,
        NDArrayLikeDict,
    )

KEYS = {
    "exn": "PEXNREF",
    "exnref": "PEXNREF",
    "rhodref": "PRHODREF",
    "pabs": "PPABSM",
    "sigs": "PSIGS",
    "cf_mf": "PCF_MF",
    "rc_mf": "PRC_MF",
    "ri_mf": "PRI_MF",
    "th": "ZRS",
    "rv": "ZRS",
    "rc": "ZRS",
    "rr": "ZRS",
    "ri": "ZRS",
    "rs": "ZRS",
    "rg": "ZRS",
    "cldfr": "PCLDFR_OUT",
    "sigqsat": None,
    "ifr": None,
    "hlc_hrc": "PHLC_HRC_OUT",
    "hlc_hcf": "PHLC_HCF_OUT",
    "hli_hri": "PHLI_HRI_OUT",
    "hli_hcf": "PHLI_HCF_OUT",
    "sigrc": None,
    "ths": "PRS",
    "rcs": "PRS",
    "rrs": "PRS",
    "ris": "PRS",
    "rss": "PRS",
    "rvs": "PRS",
    "rgs": "PRS",
    "ths": "PTHS",
}

KRR_MAPPING = {"h": 0, "v": 1, "c": 2, "r": 3, "i": 4, "s": 5, "g": 6}


######################### Ice Adjust ###########################
ice_adjust_fields_keys = [
    "sigqsat",
    "exnref",  # ref exner pression
    "exn",
    "rhodref",
    "pabs",  # absolute pressure at t
    "sigs",  # Sigma_s at time t
    "cf_mf",  # convective mass flux fraction
    "rc_mf",  # convective mass flux liquid mixing ratio
    "ri_mf",
    "th",
    "rv",
    "rc",
    "rr",
    "ri",
    "rs",
    "rg",
    "th_t",
    "sigqsat",
    "cldfr",
    "ifr",
    "hlc_hrc",
    "hlc_hcf",
    "hli_hri",
    "hli_hcf",
    "sigrc",
    "ths",
    "rvs",
    "rcs",
    "rrs",
    "ris",
    "rss",
    "rgs",
]


def allocate_state_ice_adjust(
    computational_grid: ComputationalGrid, gt4py_config: GT4PyConfig
) -> NDArrayLikeDict:
    """Allocate field to state keys following type (float, int, bool) and dimensions (2D, 3D).

    Args:
        computational_grid (ComputationalGrid): grid indexes
        gt4py_config (GT4PyConfig): gt4py configuration

    Returns:
        NDArrayLikeDict: dictionnary of field with associated keys for field name
    """

    def _allocate(
        grid_id: Tuple[DimSymbol, ...],
        units: str,
        dtype: Literal["bool", "float", "int"],
    ) -> DataArray:
        return allocate_data_array(
            computational_grid, grid_id, units, gt4py_config=gt4py_config, dtype=dtype
        )

    allocate_b_ij = partial(_allocate, grid_id=(I, J), units="", dtype="bool")
    allocate_f = partial(_allocate, grid_id=(I, J, K), units="", dtype="float")
    allocate_h = partial(_allocate, grid_id=(I, J, K - 1 / 2), units="", dtype="float")
    allocate_ij = partial(_allocate, grid_id=(I, J), units="", dtype="float")
    allocate_i_ij = partial(_allocate, grid_id=(I, J), units="", dtype="int")

    return {
        "time": datetime(year=2024, month=1, day=1),
        "sigqsat": allocate_f(),
        "exnref": allocate_f(),  # ref exner pression
        "exn": allocate_f(),
        "rhodref": allocate_f(),
        "pabs": allocate_f(),  # absolute pressure at t
        "sigs": allocate_f(),  # Sigma_s at time t
        "cf_mf": allocate_f(),  # convective mass flux fraction
        "rc_mf": allocate_f(),  # convective mass flux liquid mixing ratio
        "ri_mf": allocate_f(),
        "th": allocate_f(),
        "rv": allocate_f(),
        "rc": allocate_f(),
        "rr": allocate_f(),
        "ri": allocate_f(),
        "rs": allocate_f(),
        "rg": allocate_f(),
        "th_t": allocate_f(),
        "sigqsat": allocate_f(),
        "cldfr": allocate_f(),
        "ifr": allocate_f(),
        "hlc_hrc": allocate_f(),
        "hlc_hcf": allocate_f(),
        "hli_hri": allocate_f(),
        "hli_hcf": allocate_f(),
        "sigrc": allocate_f(),
        # tendencies
        "ths": allocate_f(),
        "rvs": allocate_f(),
        "rcs": allocate_f(),
        "rrs": allocate_f(),
        "ris": allocate_f(),
        "rss": allocate_f(),
        "rgs": allocate_f(),
    }


def get_state_ice_adjust(
    computational_grid: ComputationalGrid,
    *,
    gt4py_config: GT4PyConfig,
    netcdf_reader: NetCDFReader,
) -> DataArrayDict:
    """Create a state with reproductibility data set.

    Args:
        computational_grid (ComputationalGrid): grid
        gt4py_config (GT4PyConfig): config for gt4py
        keys (Dict[keys]): field names

    Returns:
        DataArrayDict: dictionnary of data array containing reproductibility data
    """
    state = allocate_state_ice_adjust(computational_grid, gt4py_config=gt4py_config)
    initialize_state(state, netcdf_reader)

    logging.info("Slicing for reproductibility tests with ice_adjust")
    state = slicing(state)

    return state


def slicing(state: DataArrayDict) -> DataArrayDict:
    logging.info("Slicing equivalent to PHYEX")
    new_state = {}
    for key in state.keys():
        new_state[key] = (
            state[key].isel(x=slice(0, 9472)) if key != "time" else state[key]
        )

    return new_state


def initialize_state(
    state: DataArrayDict,
    netcdreader: NetCDFReader,
) -> None:
    """Initialize fields of state dictionnary with a constant field.

    Args:
        state (DataArrayDict): dictionnary of state
        gt4py_config (GT4PyConfig): configuration of gt4py
    """

    for name, FORTRAN_NAME in KEYS.items():
        logging.info(f"name={name}, FORTRAN_NAME={FORTRAN_NAME}")
        if FORTRAN_NAME is not None:
            if FORTRAN_NAME == "ZRS":
                buffer = netcdreader.get_field(FORTRAN_NAME)[
                    :, :, KRR_MAPPING[name[-1]]
                ]

            if FORTRAN_NAME == "PRS":
                buffer = netcdreader.get_field(FORTRAN_NAME)[
                    :, :, KRR_MAPPING[name[-2]]
                ]

            elif FORTRAN_NAME not in ["ZRS", "PRS"]:
                buffer = netcdreader.get_field(FORTRAN_NAME)

        else:
            dims = netcdreader.get_dims()
            n_IJ, n_K = dims["IJ"], dims["K"]
            buffer = np.zeros((n_IJ, n_K))

        logging.info(f"name = {name}, buffer.shape = {buffer.shape}")
        initialize_field(state[name], buffer)
