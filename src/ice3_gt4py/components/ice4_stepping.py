# -*- coding: utf-8 -*-
from __future__ import annotations

from datetime import timedelta
from functools import cached_property
from itertools import repeat
from typing import Dict

from ifs_physics_common.framework.components import ImplicitTendencyComponent
from ifs_physics_common.framework.config import GT4PyConfig
from ifs_physics_common.framework.grid import ComputationalGrid, I, J, K
from ifs_physics_common.framework.storage import managed_temporary_storage
from ifs_physics_common.utils.typingx import NDArrayLikeDict, PropertyDict
from ifs_physics_common.utils.f2py import ported_method
import numpy as np


from ice3_gt4py.components.ice4_tendencies import Ice4Tendencies
from ice3_gt4py.phyex_common.phyex import Phyex


class Ice4Stepping(ImplicitTendencyComponent):
    """Component for step computation"""

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

        externals = phyex.to_externals()

        # Stencil collections

        self.ice4_stepping_heat = self.compile_stencil("ice4_stepping_heat", externals)
        self.ice4_limiter = self.compile_stencil("limiter", externals)

        self.ice4_state_update = self.compile_stencil("state_update", externals)
        self.external_tendencies_update = self.compile_stencil(
            "external_tendencies_update", externals
        )
        self.tmicro_init = self.compile_stencil("ice4_stepping_tmicro_init", externals)
        self.tsoft_init = self.compile_stencil("ice4_stepping_tsoft_init", externals)
        self.ldcompute_init = self.compile_stencil(
            "ice4_stepping_ldcompute_init", externals
        )

        # Component for tendency update
        self.ice4_tendencies_heavy = Ice4Tendencies(
            computational_grid=self.computational_grid,
            gt4py_config=self.gt4py_config,
            phyex=phyex,
            ldsoft=False,
        )

        self.ice4_tendencies_light = Ice4Tendencies(
            computational_grid=self.computational_grid,
            gt4py_config=self.gt4py_config,
            phyex=phyex,
            ldsoft=True,
        )

    @cached_property
    def _input_properties(self) -> PropertyDict:
        return {
            "ldmicro": {"grid": (I, J, K), "units": ""},
            "rhodref": {"grid": (I, J, K), "units": ""},
            "pabs_t": {"grid": (I, J, K), "units": ""},
            "exn": {"grid": (I, J, K), "units": ""},
            "cf": {"grid": (I, J, K), "units": ""},
            "sigma_rc": {"grid": (I, J, K), "units": ""},
            "ci_t": {"grid": (I, J, K), "units": ""},
            # "ai": {"grid": (I, J, K), "units": ""},
            # "cj": {"grid": (I, J, K), "units": ""},
            # "ssi": {"grid": (I, J, K), "units": ""},      zssi (mode_ice4_stepping.F90)
            "hlc_hcf": {"grid": (I, J, K), "units": ""},
            # "hlc_lcf": {"grid": (I, J, K), "units": ""},  cf = hlc_hcf + hlc_lcf /// hlc_lcf is the proportion of low cloud fraction in the grid
            "hlc_hrc": {"grid": (I, J, K), "units": ""},
            # "hlc_lrc": {"grid": (I, J, K), "units": ""},  rc = hlc_hrc + hlc_lrc /// Low LWC in the grid
            "hli_hcf": {"grid": (I, J, K), "units": ""},
            # "hli_lcf": {"grid": (I, J, K), "units": ""},  for iceclouds
            "hli_hri": {"grid": (I, J, K), "units": ""},
            # "hli_lri": {"grid": (I, J, K), "units": ""},  for ice content
            # "fr": {"grid": (I, J, K), "units": ""},
            "th_t": {"grid": (I, J, K), "units": ""},
            "ls_fact": {"grid": (I, J, K), "units": ""},
            "lv_fact": {"grid": (I, J, K), "units": ""},
            "t": {"grid": (I, J, K), "units": ""},
            "rv_t": {"grid": (I, J, K), "units": ""},
            "rc_t": {"grid": (I, J, K), "units": ""},
            "rr_t": {"grid": (I, J, K), "units": ""},
            "ri_t": {"grid": (I, J, K), "units": ""},
            "rs_t": {"grid": (I, J, K), "units": ""},
            "rg_t": {"grid": (I, J, K), "units": ""},
        }

    @cached_property
    def _tendency_properties(self) -> PropertyDict:
        return {}

    @cached_property
    def _diagnostic_properties(self) -> PropertyDict:
        return {}

    @cached_property
    def _temporaries(self) -> PropertyDict:
        return {}

    @ported_method(
        from_file="PHYEX/src/common/micro/mode_ice4_stepping.F90",
        from_line=214,
        to_line=438,
    )
    def array_call(
        self,
        state: NDArrayLikeDict,
        timestep: timedelta,
        out_tendencies: NDArrayLikeDict,
        out_diagnostics: NDArrayLikeDict,
        overwrite_tendencies: Dict[str, bool],
    ):

        with managed_temporary_storage(
            self.computational_grid,
            *repeat(((I, J, K), "bool"), 1),
            *repeat(((I, J, K), "float"), 28),
            gt4py_config=self.gt4py_config,
        ) as (
            # masks
            ldcompute,
            # intial mixing ratios
            rc_0r_t,
            rr_0r_t,
            ri_0r_t,
            rs_0r_t,
            rg_0r_t,
            # increments
            theta_b,
            rv_b,
            rc_b,
            rr_b,
            ri_b,
            rs_b,
            rg_b,
            # tnd update
            theta_a_tnd,
            rv_a_tnd,
            rc_a_tnd,
            rr_a_tnd,
            ri_a_tnd,
            rs_a_tnd,
            rg_a_tnd,
            # tendances externes
            theta_ext_tnd,
            rc_ext_tnd,
            rr_ext_tnd,
            ri_ext_tnd,
            rs_ext_tnd,
            rg_ext_tnd,
            # timing
            t_micro,
            delta_t_micro,
        ):
            # TODO : add first stencils
            # TODO : l192 -> l206 if LDAUCV_ADJU
            # TODO : l187 -> l191 if LDSIGMA_RC
            # TODO : l180 -> l186 if LDEXT_TEND

            # Translation note : Ice4Stepping is implemented assuming PARAMI%XTSTEP_TS = 0
            #                   l225 to l229 omitted
            #                   l334 to l341 omitted
            #                   l174 to l178 omitted

            ############## t_micro_init ################
            dt = timestep.total_seconds()

            state_tmicro_init = {"ldmicro": state["ldmicro"], "t_micro": t_micro}
            self.tmicro_init(**state_tmicro_init)

            ############# First iteration ##############

            ############################################
            outerloop_counter = 0
            lsoft = False

            # Translation note : l230 to l 237 in Fortran
            self.ldcompute_init(ldcompute, t_micro)

            # Iterations limiter
            while np.any(ldcompute[...]):

                # loop breaker

                ####### ice4_stepping_heat #############
                state_stepping_heat = {
                    key: state[key]
                    for key in [
                        "rv_t",
                        "rc_t",
                        "rr_t",
                        "ri_t",
                        "rs_t",
                        "rg_t",
                        "exn",
                        "th_t",
                        "ls_fact",
                        "lv_fact",
                        "t",
                    ]
                }

                self.ice4_stepping_heat(**state_stepping_heat)

                ####### tendencies #######
                state_ice4_tendencies = {
                    **{
                        key: state[key]
                        for key in [
                            "rhodref",
                            "exn",
                            "ls_fact",
                            "lv_fact",
                            "rv_t",
                            "cf",
                            "sigma_rc",
                            "ci_t",
                            "ai",
                            "cj",
                            "ssi",
                            "t",
                            "th_t",
                            "rv_t",
                            "rc_t",
                            "rr_t",
                            "ri_t",
                            "rs_t",
                            "rg_t",
                            "hlc_hcf",
                            "hlc_lcf",
                            "hlc_hrc",
                            "hlc_lrc",
                            "hli_hcf",
                            "hli_lcf",
                            "hli_hri",
                            "hli_lri",
                            "fr",
                        ]
                    },
                    **{"pres": state["pabs_t"]},
                    **{
                        "ldcompute": ldcompute,
                        "theta_tnd": theta_a_tnd,
                        "rv_tnd": rv_a_tnd,
                        "rc_tnd": rc_a_tnd,
                        "rr_tnd": rr_a_tnd,
                        "ri_tnd": ri_a_tnd,
                        "rs_tnd": rs_a_tnd,
                        "rg_tnd": rg_a_tnd,
                        "theta_increment": theta_b,
                        "rv_increment": rv_b,
                        "rc_increment": rc_b,
                        "rr_increment": rr_b,
                        "ri_increment": ri_b,
                        "rs_increment": rs_b,
                        "rg_increment": rg_b,
                    },
                }

                self.ice4_tendencies_light.array_call(
                    ldsoft=lsoft,
                    state=state_ice4_tendencies,
                    timestep=dt,
                    out_diagnostics={},
                    out_tendencies={},
                    overwrite_tendencies={},
                )

                # Translation note : l277 to l283 omitted, no external tendencies in AROME

                ######### ice4_step_limiter ############################
                state_limiter = {
                    **{
                        key: state[key]
                        for key in [
                            "exn",
                            "theta_t",
                            "rc_t",
                            "rr_t",
                            "ri_t",
                            "rs_t",
                            "rg_t",
                            "delta_t_soft",
                            "t_soft",
                        ]
                    },
                    **{
                        "theta_a_tnd": theta_a_tnd,
                        "theta_b": theta_b,
                        "theta_ext_tnd": theta_ext_tnd,
                        "rc_a_tnd": rc_a_tnd,
                        "rr_a_tnd": rr_a_tnd,
                        "ri_a_tnd": ri_a_tnd,
                        "rs_a_tnd": rs_a_tnd,
                        "rg_a_tnd": rg_a_tnd,
                        "rc_ext_tnd": rc_ext_tnd,
                        "rr_ext_tnd": rr_ext_tnd,
                        "ri_ext_tnd": ri_ext_tnd,
                        "rs_ext_tnd": rs_ext_tnd,
                        "rg_ext_tnd": rg_ext_tnd,
                        "rc_b": rc_b,
                        "rr_b": rr_b,
                        "ri_b": ri_b,
                        "rs_b": rs_b,
                        "rg_b": rg_b,
                        "rc_0r_t": rc_0r_t,
                        "rr_0r_t": rr_0r_t,
                        "ri_0r_t": ri_0r_t,
                        "rs_0r_t": rs_0r_t,
                        "rg_0r_t": rg_0r_t,
                        "delta_t_micro": delta_t_micro,
                        "t_micro": t_micro,
                        "ldcompute": ldcompute,
                    },
                }

                self.ice4_limiter(**state_limiter)

                # l394 to l404
                # 4.7 new values for next iteration
                ############### ice4_state_update ######################
                state_state_update = {
                    **{
                        key: state[key]
                        for key in [
                            "th_t",
                            "rc_t",
                            "rr_t",
                            "ri_t",
                            "rs_t",
                            "rg_t",
                            "ci_t",
                            "ldmicro",
                        ]
                    },
                    **{
                        "theta_a_tnd": theta_a_tnd,
                        "rc_a_tnd": rc_a_tnd,
                        "rr_a_tnd": rr_a_tnd,
                        "ri_a_tnd": ri_a_tnd,
                        "rs_a_tnd": rs_a_tnd,
                        "rg_a_tnd": rg_a_tnd,
                        "theta_b": theta_b,
                        "rc_b": rc_b,
                        "rr_b": rr_b,
                        "ri_b": ri_b,
                        "rs_b": rs_b,
                        "rg_b": rg_b,
                        "delta_t_micro": delta_t_micro,
                        "t_micro": t_micro,
                    },
                }
                self.ice4_state_update(**state_state_update)

                # TODO : next loop
                lsoft = True
                outerloop_counter += 1

            ################ external_tendencies_update ############
            # Translation note : l440 to l452 removed
            # No external tendencies in AROME
