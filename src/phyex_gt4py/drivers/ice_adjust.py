# -*- coding: utf-8 -*-
from __future__ import annotations
from dataclasses import asdict
from datetime import timedelta
from functools import cached_property
from itertools import repeat
from typing import Dict

from ifs_physics_common.framework.grid import ComputationalGrid
from ifs_physics_common.framework.config import GT4PyConfig
from ifs_physics_common.framework.components import ImplicitTendencyComponent
from ifs_physics_common.utils.typingx import PropertyDict, NDArrayLikeDict
from ifs_physics_common.framework.grid import I, J, K
from ifs_physics_common.framework.storage import managed_temporary_storage
from phyex_gt4py.phyex_common.phyex import Phyex


class IceAdjust(ImplicitTendencyComponent):
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

        externals = {}
        externals.update(asdict(phyex.nebn))
        externals.update(asdict(phyex.cst))
        externals.update(asdict(phyex.param_icen))
        externals.update(
            {
                "nrr": 6,
                "criautc": 0,
                "acriauti": 0,
                "bcriauti": 0,
                "criauti": 0,
                "tstep": 1,
            }
        )

        self.ice_adjust = self.compile_stencil("ice_adjust", externals)

    @cached_property
    def _input_properties(self) -> PropertyDict:
        return {
            "sigqsat": {
                "grid": (I, J, K),
                "units": "",
            },  # coeff applied to qsat variance
            "exnref": {"grid": (I, J, K), "units": ""},  # ref exner pression
            "exn": {"grid": (I, J, K), "units": ""},
            "rhodref": {"grid": (I, J, K), "units": ""},  #
            "pabs": {"grid": (I, J, K), "units": ""},  # absolute pressure at t
            "sigs": {"grid": (I, J, K), "units": ""},  # Sigma_s at time t
            "cf_mf": {"grid": (I, J, K), "units": ""},  # convective mass flux fraction
            "rc_mf": {
                "grid": (I, J, K),
                "units": "",
            },  # convective mass flux liquid mixing ratio
            "ri_mf": {"grid": (I, J, K), "units": ""},
        }

    @cached_property
    def _tendency_properties(self) -> PropertyDict:
        return {
            "ths": {"grid": (I, J, K), "units": ""},
            "rvs": {"grid": (I, J, K), "units": ""},  # PRS(1)
            "rcs": {"grid": (I, J, K), "units": ""},  # PRS(2)
            "ris": {"grid": (I, J, K), "units": ""},  # PRS(4)
            "th": {"grid": (I, J, K), "units": ""},  # ZRS(0)
            "rv": {"grid": (I, J, K), "units": ""},  # ZRS(1)
            "rc": {"grid": (I, J, K), "units": ""},  # ZRS(2)
            "rr": {"grid": (I, J, K), "units": ""},  # ZRS(3)
            "ri": {"grid": (I, J, K), "units": ""},  # ZRS(4)
            "rs": {"grid": (I, J, K), "units": ""},  # ZRS(5)
            "rg": {"grid": (I, J, K), "units": ""},  # ZRS(6)
            "cldfr": {"grid": (I, J, K), "units": ""},
            "ice_cld_wgt": {"grid": (I, J, K), "units": ""},
        }

    @cached_property
    def _diagnostic_properties(self) -> PropertyDict:
        return {
            "icldfr": {"grid": (I, J, K), "units": ""},
            "wcldfr": {"grid": (I, J, K), "units": ""},
            "ssio": {"grid": (I, J, K), "units": ""},
            "ssiu": {"grid": (I, J, K), "units": ""},
            "srcs": {"grid": (I, J, K), "units": ""},
            "ifr": {"grid": (I, J, K), "units": ""},
            "hlc_hrc": {"grid": (I, J, K), "units": ""},
            "hlc_hcf": {"grid": (I, J, K), "units": ""},
            "hli_hri": {"grid": (I, J, K), "units": ""},
            "hli_hcf": {"grid": (I, J, K), "units": ""},
        }

    @cached_property
    def _temporaries(self) -> PropertyDict:
        return {
            "rt": {
                "grid": (I, J, K),
                "units": "",
            },  # work array for total water mixing ratio
            "pv": {"grid": (I, J, K), "units": ""},  # thermodynamics
            "piv": {"grid": (I, J, K), "units": ""},  # thermodynamics
            "qsl": {"grid": (I, J, K), "units": ""},  # thermodynamics
            "qsi": {"grid": (I, J, K), "units": ""},
            "frac_tmp": {"grid": (I, J, K), "units": ""},  # ice fraction
            "cond_tmp": {"grid": (I, J, K), "units": ""},  # condensate
            "a": {"grid": (I, J, K), "units": ""},  # related to computation of Sig_s
            "sbar": {"grid": (I, J, K), "units": ""},
            "sigma": {"grid": (I, J, K), "units": ""},
            "q1": {"grid": (I, J, K), "units": ""},
        }

    def array_call(
        self,
        state: NDArrayLikeDict,
        timestep: timedelta,
        out_tendencies: NDArrayLikeDict,
        out_diagnostics: NDArrayLikeDict,
        overwrite_tendencies: Dict[str, bool],
    ) -> None:

        with managed_temporary_storage(
            self.computational_grid,
            *repeat(((I, J, K), "float"), 21),
            gt4py_config=self.gt4py_config,
        ) as (
            rt,
            pv,
            piv,
            qsl,
            qsi,
            frac_tmp,
            cond_tmp,
            a,
            sbar,
            sigma,
            q1,
            icldfr,
            wcldfr,
            ssio,
            ssiu,
            srcs,
            ifr,
            hlc_hrc,
            hlc_hcf,
            hli_hri,
            hli_hcf,
        ):
            inputs = {name: state[name] for name in self.input_properties}
            tendencies = {
                name: out_tendencies[name] for name in self.tendency_properties
            }
            diagnostics = {
                name: out_diagnostics[name] for name in self.diagnostic_properties
            }
            temporaries = {
                "rt": rt,
                "pv": pv,
                "piv": piv,
                "qsl": qsl,
                "qsi": qsi,
                "frac_tmp": frac_tmp,
                "cond_tmp": cond_tmp,
                "a": a,
                "sbar": sbar,
                "sigma": sigma,
                "q1": q1,
            }

            self.ice_adjust(
                **inputs,
                **tendencies,
                **diagnostics,
                **temporaries,
                dt=timestep.total_seconds(),
                origin=(0, 0, 0),
                domain=self.computational_grid.grids[I, J, K].shape,
                validate_args=self.gt4py_config.validate_args,
                exec_info=self.gt4py_config.exec_info,
            )
