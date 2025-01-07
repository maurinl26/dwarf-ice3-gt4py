# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
from datetime import timedelta
import sys
from functools import cached_property
from itertools import repeat
from typing import Dict
from ifs_physics_common.framework.components import ImplicitTendencyComponent
from ifs_physics_common.framework.config import GT4PyConfig
from ifs_physics_common.framework.grid import ComputationalGrid, I, J, K
from ifs_physics_common.framework.storage import managed_temporary_storage
from ifs_physics_common.utils.typingx import NDArrayLikeDict, PropertyDict

from ice3_gt4py.phyex_common.phyex import Phyex

logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
logging.getLogger()


class IceAdjustSplit(ImplicitTendencyComponent):
    """Implicit Tendency Component calling
    ice_adjust : saturation adjustment of temperature and mixing ratios

    ice_adjust stencil is ice_adjust.F90 in PHYEX
    """

    def __init__(
        self,
        computational_grid: ComputationalGrid,
        gt4py_config: GT4PyConfig,
        phyex: Phyex,
        *,
        enable_checks: bool = True,
    ) -> None:
        super().__init__(
            computational_grid, enable_checks=enable_checks, gt4py_config=gt4py_config
        )

        self.externals = phyex.to_externals()
        self.thermo = self.compile_stencil("thermodynamic_fields", self.externals)
        self.condensation = self.compile_stencil("condensation", self.externals)
        self.cloud_fraction = self.compile_stencil("cloud_fraction", self.externals)

        logging.info(f"IceAdjustSplit - Keys")
        logging.info(f"LSUBG_COND : {phyex.nebn.LSUBG_COND}")
        logging.info(f"LSIGMAS :  {phyex.nebn.LSIGMAS}")
        logging.info(f"FRAC_ICE_ADJUST : {phyex.nebn.FRAC_ICE_ADJUST}")
        logging.info(f"CONDENS : {phyex.nebn.CONDENS}")
        logging.info(f"LAMBDA3 : {phyex.nebn.LAMBDA3}")
        logging.info(f"OCOMPUTE_SRC absent")
        logging.info(f"LMFCONV : {phyex.LMFCONV}")
        logging.info(f"LOCND2 absent")
        logging.info(f"LHGT_QS : {phyex.nebn.LHGT_QS}")
        logging.info(f"LSTATNW : {phyex.nebn.LSTATNW}")
        logging.info(f"SUBG_MF_PDF : {phyex.param_icen.SUBG_MF_PDF}")

        logging.info(f"Constants for condensation")
        logging.info(f"RD : {phyex.cst.RD}")
        logging.info(f"RV : {phyex.cst.RV}")

        logging.info(f"Constants for thermodynamic fields")
        logging.info(f"CPD : {phyex.cst.CPD}")
        logging.info(f"CPV : {phyex.cst.CPV}")
        logging.info(f"CL : {phyex.cst.CL}")
        logging.info(f"CI : {phyex.cst.CI}")

        logging.info(f"Constants for cloud fraction")
        logging.info(f"CRIAUTC : {phyex.rain_ice_param.CRIAUTC}")
        logging.info(f"CRIAUTI : {phyex.rain_ice_param.CRIAUTI}")
        logging.info(f"ACRIAUTI : {phyex.rain_ice_param.ACRIAUTI}")
        logging.info(f"BCRIAUTI : {phyex.rain_ice_param.BCRIAUTI}")
        logging.info(f"TT : {phyex.cst.TT}")

    @cached_property
    def _input_properties(self) -> PropertyDict:
        return {
            "sigqsat": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": None,
                "dtype": "float",
            },
            "exn": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PEXNREF",
                "dtype": "float",
            },
            "exnref": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PEXNREF",
                "dtype": "float",
            },
            "rhodref": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PRHODREF",
                "dtype": "float",
            },
            "pabs": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PPABSM",
                "dtype": "float",
            },
            "sigs": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PSIGS",
                "dtype": "float",
            },
            "cf_mf": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PCF_MF",
                "dtype": "float",
            },
            "rc_mf": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PRC_MF",
                "dtype": "float",
            },
            "ri_mf": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PRI_MF",
                "dtype": "float",
            },
            "th": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "ZRS",
                "irr": 0,
                "dtype": "float",
            },
            "rv": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "ZRS",
                "irr": 1,
                "dtype": "float",
            },
            "rc": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "ZRS",
                "irr": 2,
                "dtype": "float",
            },
            "rr": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "ZRS",
                "irr": 3,
                "dtype": "float",
            },
            "ri": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "ZRS",
                "irr": 4,
                "dtype": "float",
            },
            "rs": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "ZRS",
                "irr": 5,
                "dtype": "float",
            },
            "rg": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "ZRS",
                "irr": 6,
                "dtype": "float",
            },
            "cldfr": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": None,
                "dtype": "float",
            },
            "ifr": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": None,
                "dtype": "float",
            },
            "hlc_hrc": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": None,
                "dtype": "float",
            },
            "hlc_hcf": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": None,
                "dtype": "float",
            },
            "hli_hri": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": None,
                "dtype": "float",
            },
            "hli_hcf": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": None,
                "dtype": "float",
            },
            "sigrc": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": None,
                "dtype": "float",
            },
            "ths": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PRS",
                "irr": 0,
                "dtype": "float",
            },
            "rcs": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PRS",
                "irr": 1,
                "dtype": "float",
            },
            "rrs": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PRS",
                "irr": 2,
                "dtype": "float",
            },
            "ris": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PRS",
                "irr": 3,
                "dtype": "float",
            },
            "rss": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PRS",
                "irr": 4,
                "dtype": "float",
            },
            "rvs": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PRS",
                "irr": 5,
                "dtype": "float",
            },
            "rgs": {
                "grid": (I, J, K),
                "units": "",
                "fortran_name": "PRS",
                "irr": 6,
                "dtype": "float",
            },
        }

    @cached_property
    def _tendency_properties(self) -> PropertyDict:
        return {}

    @cached_property
    def _diagnostic_properties(self) -> PropertyDict:
        return {}

    @cached_property
    def _temporaries(self) -> PropertyDict:
        return {
            "lv": {"grid": (I, J, K), "units": ""},
            "ls": {"grid": (I, J, K), "units": ""},
            "cph": {"grid": (I, J, K), "units": ""},
        }

    def array_call(
        self,
        state: NDArrayLikeDict,
        timestep: timedelta,
        out_tendencies: NDArrayLikeDict,
        out_diagnostics: NDArrayLikeDict,
        overwrite_tendencies: Dict[str, bool],
    ) -> None:

        logging.info(f"Timestep : {timestep.total_seconds()}")

        with managed_temporary_storage(
            self.computational_grid,
            *repeat(((I, J, K), "float"), 8),
            gt4py_config=self.gt4py_config,
        ) as (lv, ls, cph, t, rc_out, ri_out, rv_out, t_out):

            state_thermo = {
                key: state[key]
                for key in [
                    "th",
                    "exn",
                    "rv",
                    "rc",
                    "rr",
                    "ri",
                    "rs",
                    "rg",
                ]
            }

            temporaries_thermo = {"cph": cph, "lv": lv, "ls": ls, "t": t}

            logging.info("Launching thermo")
            # logging.debug(f"Thermo __externals__ : {
            #     self.externals["NRR"],
            #     self.externals["CPV"],
            #     self.externals["CPD"],
            #     self.externals["CL"],
            #     self.externals["CI"]}")

            self.thermo(
                **state_thermo,
                **temporaries_thermo,
                origin=(0, 0, 0),
                domain=self.computational_grid.grids[I, J, K].shape,
                validate_args=self.gt4py_config.validate_args,
                exec_info=self.gt4py_config.exec_info,
            )

            state_condensation = {
                **{key: state[key] for key in ["sigqsat", "pabs", "cldfr", "sigs"]},
                **{
                    "rv_in": state["rv"],
                    "ri_in": state["ri"],
                    "rc_in": state["rc"],
                },
            }

            temporaries_condensation = {
                "cph": cph,
                "lv": lv,
                "ls": ls,
                "t": t,
                "rv_out": rv_out,
                "ri_out": ri_out,
                "rc_out": rc_out,
            }

            logging.info("Launching condensation")
            self.condensation(
                **state_condensation,
                **temporaries_condensation,
                origin=(0, 0, 0),
                domain=self.computational_grid.grids[I, J, K].shape,
                validate_args=self.gt4py_config.validate_args,
                exec_info=self.gt4py_config.exec_info,
            )

            # TODO : insert diagnostic on sigrc here
            # Translation note : if needed, sigrc_computation could be inserted here

            state_cloud_fraction = {
                key: state[key]
                for key in [
                    # "t",
                    "rhodref",
                    "exnref",
                    "rc",
                    "ri",
                    "ths",
                    "rvs",
                    "rcs",
                    "ris",
                    "rc_mf",
                    "ri_mf",
                    "cf_mf",
                    "hlc_hrc",
                    "hlc_hcf",
                    "hli_hri",
                    "hli_hcf",
                    "cldfr",
                ]
            }

            # TODO: check the scope of t
            temporaries_cloud_fraction = {
                "lv": lv,
                "ls": ls,
                "cph": cph,
                "t": t,
                "rc_tmp": rc_out,
                "ri_tmp": ri_out,
            }

            logging.info("Launching cloud fraction")
            self.cloud_fraction(
                **state_cloud_fraction,
                **temporaries_cloud_fraction,
                dt=timestep.total_seconds(),
                origin=(0, 0, 0),
                domain=self.computational_grid.grids[I, J, K].shape,
                validate_args=self.gt4py_config.validate_args,
                exec_info=self.gt4py_config.exec_info,
            )
