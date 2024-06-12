# -*- coding: utf-8 -*-
from __future__ import annotations

import logging
from dataclasses import asdict
from datetime import timedelta
from functools import cached_property
from itertools import repeat
from typing import Dict
from gt4py.storage import from_array
import numpy as np
from ice3_gt4py.phyex_common.xker_raccs import ker_raccs, ker_raccss, ker_saccrg
from ice3_gt4py.phyex_common.xker_sdryg import ker_sdryg
from ice3_gt4py.phyex_common.xker_rdryg import ker_rdryg

from ifs_physics_common.framework.components import ImplicitTendencyComponent
from ifs_physics_common.framework.config import GT4PyConfig
from ifs_physics_common.framework.grid import ComputationalGrid, I, J, K
from ifs_physics_common.framework.storage import managed_temporary_storage
from ifs_physics_common.utils.typingx import NDArrayLikeDict, PropertyDict

from ice3_gt4py.phyex_common.phyex import Phyex


class Ice4Tendencies(ImplicitTendencyComponent):
    """Implicit Tendency Component calling
    ice_adjust : saturation adjustment of temperature and mixing ratios

    ice_adjust stencil is ice_adjust.F90 in PHYEX
    """

    def __init__(
        self,
        computational_grid: ComputationalGrid,
        gt4py_config: GT4PyConfig,
        phyex: Phyex,
        ldsoft: bool,
        *,
        enable_checks: bool = True,
    ) -> None:
        super().__init__(
            computational_grid, enable_checks=enable_checks, gt4py_config=gt4py_config
        )

        externals = phyex.to_externals()

        self.gaminc_rim1 = phyex.rain_ice_param.GAMINC_RIM1
        self.gaminc_rim2 = phyex.rain_ice_param.GAMINC_RIM2
        self.gaminc_rim4 = phyex.rain_ice_param.GAMINC_RIM4

        self.ker_raccs = phyex.rain_ice_param.ker_raccs
        self.ker_raccss = phyex.rain_ice_param.ker_raccss
        self.ker_saccrg = phyex.rain_ice_param.ker_saccrg
        self.ker_sdryg = phyex.rain_ice_param.ker_sdryg
        self.ker_rdryg = phyex.rain_ice_param.ker_rdryg

        # Tendencies
        self.ice4_nucleation = self.compile_stencil("ice4_nucleation", externals)
        self.ice4_nucleation_post_processing = self.compile_stencil(
            "ice4_nucleation_post_processing", externals
        )

        self.ice4_rrhong = self.compile_stencil("ice4_rrhong", externals)
        self.ice4_rrhong_post_processing = self.compile_stencil(
            "ice4_rrhong_post_processing", externals
        )

        self.ice4_rimltc = self.compile_stencil("ice4_rimltc", externals)
        self.ice4_rimltc_post_processing = self.compile_stencil(
            "ice4_rimltc_post_processing"
        )

        self.ice4_increment_update = self.compile_stencil(
            "ice4_increment_update", externals
        )
        self.ice4_derived_fields = self.compile_stencil(
            "ice4_derived_fields", externals
        )

        # TODO: add ice4_compute_pdf
        self.ice4_slope_parameters = self.compile_stencil(
            "ice4_slope_parameters", externals
        )

        externals.update({"LDSOFT": False})
        self.ice4_slow_heavy = self.compile_stencil("ice4_slow", externals)
        self.ice4_warm_heavy = self.compile_stencil("ice4_warm", externals)
        self.ice4_fast_rs_heavy = self.compile_stencil("ice4_fast_rs", externals)
        self.ice4_fast_rg_heavy = self.compile_stencil("ice4_fast_rg", externals)
        self.ice4_fast_ri_heavy = self.compile_stencil("ice4_fast_ri", externals)

        externals.update({"LDSOFT": True})
        self.ice4_slow_light = self.compile_stencil("ice4_slow", externals)
        self.ice4_warm_light = self.compile_stencil("ice4_warm", externals)
        self.ice4_fast_rs_light = self.compile_stencil("ice4_fast_rs", externals)
        self.ice4_fast_rg_light = self.compile_stencil("ice4_fast_rg", externals)
        self.ice4_fast_ri_light = self.compile_stencil("ice4_fast_ri", externals)

        self.ice4_fast_rg_pre_processing = self.compile_stencil(
            "ice4_fast_rg_pre_processing", externals
        )

        self.ice4_tendencies_update = self.compile_stencil(
            "ice4_tendencies_update", externals
        )

        self.ice4_compute_pdf = self.compile_stencil("ice4_compute_pdf", externals)

    @cached_property
    def _input_properties(self) -> PropertyDict:
        return {
            "ldcompute": {"grid": (I, J, K), "units": ""},
            "pres": {"grid": (I, J, K), "units": ""},
            "rhodref": {"grid": (I, J, K), "units": ""},
            "exn": {"grid": (I, J, K), "units": ""},
            "ls_fact": {"grid": (I, J, K), "units": ""},
            "lv_fact": {"grid": (I, J, K), "units": ""},
            "rv_t": {"grid": (I, J, K), "units": ""},
            "cf": {"grid": (I, J, K), "units": ""},
            "sigma_rc": {"grid": (I, J, K), "units": ""},
            "ci_t": {"grid": (I, J, K), "units": ""},
            "ai": {"grid": (I, J, K), "units": ""},
            "cj": {"grid": (I, J, K), "units": ""},
            "ssi": {"grid": (I, J, K), "units": ""},
            "t": {"grid": (I, J, K), "units": ""},
            "th_t": {"grid": (I, J, K), "units": ""},
            # PVART in f90
            "rv_t": {"grid": (I, J, K), "units": ""},
            "rc_t": {"grid": (I, J, K), "units": ""},
            "rr_t": {"grid": (I, J, K), "units": ""},
            "ri_t": {"grid": (I, J, K), "units": ""},
            "rs_t": {"grid": (I, J, K), "units": ""},
            "rg_t": {"grid": (I, J, K), "units": ""},
            # A in f90
            "theta_tnd": {"grid": (I, J, K), "units": ""},
            "rv_tnd": {"grid": (I, J, K), "units": ""},
            "rc_tnd": {"grid": (I, J, K), "units": ""},
            "rr_tnd": {"grid": (I, J, K), "units": ""},
            "ri_tnd": {"grid": (I, J, K), "units": ""},
            "rs_tnd": {"grid": (I, J, K), "units": ""},
            "rg_tnd": {"grid": (I, J, K), "units": ""},
            # B in f90 TODO : in diagnostics
            "theta_increment": {"grid": (I, J, K), "units": ""},
            "rv_increment": {"grid": (I, J, K), "units": ""},
            "rc_increment": {"grid": (I, J, K), "units": ""},
            "rr_increment": {"grid": (I, J, K), "units": ""},
            "ri_increment": {"grid": (I, J, K), "units": ""},
            "rs_increment": {"grid": (I, J, K), "units": ""},
            "rg_increment": {"grid": (I, J, K), "units": ""},
            # others
            "hlc_hcf": {"grid": (I, J, K), "units": ""},
            "hlc_lcf": {"grid": (I, J, K), "units": ""},
            "hlc_hrc": {"grid": (I, J, K), "units": ""},
            "hlc_lrc": {"grid": (I, J, K), "units": ""},
            "hli_hcf": {"grid": (I, J, K), "units": ""},
            "hli_lcf": {"grid": (I, J, K), "units": ""},
            "hli_hri": {"grid": (I, J, K), "units": ""},
            "hli_lri": {"grid": (I, J, K), "units": ""},
            "fr": {"grid": (I, J, K), "units": ""},
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

    def array_call(
        self,
        ldsoft: bool,
        state: NDArrayLikeDict,
        timestep: timedelta,
        out_tendencies: NDArrayLikeDict,
        out_diagnostics: NDArrayLikeDict,
        # overwrite_tendencies: OptionalDict[str, bool],
    ) -> None:
        with managed_temporary_storage(
            self.computational_grid,
            *repeat(((I, J, K), "float"), 63),
            *repeat(((I, J, K), "int"), 4),
            gt4py_config=self.gt4py_config,
        ) as (
            # mr
            rvheni_mr,
            rrhong_mr,
            rimltc_mr,
            rgsi_mr,
            rsrimcg_mr,
            # slopes
            lbdar,
            lbdar_rf,
            lbdas,
            lbdag,
            # tnd
            rc_honi_tnd,
            rv_deps_tnd,
            ri_aggs_tnd,
            ri_auts_tnd,
            rv_depg_tnd,
            rs_mltg_tnd,
            rc_mltsr_tnd,
            rs_rcrims_tnd,
            rs_rcrimss_tnd,
            rs_rsrimcg_tnd,
            rs_rraccs_tnd,
            rs_rraccss_tnd,
            rs_rsaccrg_tnd,
            rs_freez1_tnd,
            rs_freez2_tnd,
            rg_rcdry_tnd,
            rg_ridry_tnd,
            rg_rsdry_tnd,
            rg_rrdry_tnd,
            rg_riwet_tnd,
            rg_rswet_tnd,
            rg_freez1_tnd,
            rg_freez2_tnd,
            rc_beri_tnd,
            # transfos
            rgsi,
            rchoni,
            rvdeps,
            riaggs,
            riauts,
            rvdepg,
            rcautr,
            rcaccr,
            rrevav,
            rcberi,
            rsmltg,
            rcmltsr,
            rraccss,  # 13
            rraccsg,  # 14
            rsaccrg,  # 15
            rcrimss,  # 16
            rcrimsg,  # 17
            rsrimcg,  # 18
            ricfrrg,  # 19
            rrcfrig,  # 20
            ricfrr,  # 21
            rcwetg,  # 22
            riwetg,  # 23
            rrwetg,  # 24
            rswetg,  # 25
            rcdryg,  # 26
            ridryg,  # 27
            rrdryg,  # 28
            rsdryg,  # 29
            rgmltr,  # 31
            index_floor,
            index_floor_r,
            index_floor_s,
            index_floor_g,
        ):

            logging.info(f"{state['th_t']}, {state['rhodref']}")

            ############## ice4_nucleation ################
            state_nucleation = {
                **{
                    key: state[key]
                    for key in [
                        "th_t",
                        "rhodref",
                        "exn",
                        "ls_fact",
                        "t",
                        "rv_t",
                        "ci_t",
                        "ssi",
                        "ldcompute",
                    ]
                },
                **{"pabs_t": state["pres"]},
            }

            temporaries_nucleation = {
                "rvheni_mr": rvheni_mr,
            }

            # timestep
            self.ice4_nucleation(**state_nucleation, **temporaries_nucleation)

            ############## ice4_nucleation_post_processing ####################

            state_nucleation_pp = {
                **{
                    key: state[key]
                    for key in [
                        "t",
                        "exn",
                        "ls_fact",
                        "lv_fact",
                        "th_t",
                        "rv_t",
                        "ri_t",
                    ]
                },
            }

            tmps_nucleation_pp = {"rvheni_mr": rvheni_mr}

            # Timestep
            self.ice4_nucleation_post_processing(
                **state_nucleation_pp, **tmps_nucleation_pp
            )

            ########################### ice4_rrhong #################################
            state_rrhong = {
                "ldcompute": state["ldcompute"],
                **{
                    key: state[key]
                    for key in ["t", "exn", "lv_fact", "ls_fact", "th_t", "rr_t"]
                },
            }

            tmps_rrhong = {"rrhong_mr": rrhong_mr}

            self.ice4_rrhong(**state_rrhong, **rrhong_mr)

            ########################### ice4_rrhong_post_processing #################
            state_rrhong_pp = {
                **{
                    key: state[key]
                    for key in [
                        "t",
                        "exn",
                        "lv_fact",
                        "ls_fact",
                        "th_t",
                        "rg_t",
                        "rr_t",
                    ]
                },
            }

            self.ice4_rrhong_post_processing(**state_rrhong_pp, **tmps_rrhong)

            ########################## ice4_rimltc ##################################
            state_rimltc = {
                "ldcompute": state["ldcompute"],
                **{
                    key: state[key]
                    for key in [
                        "t",
                        "exn",
                        "lv_fact",
                        "ls_fact",
                        "th_t",
                        "ri_t",
                    ]
                },
            }

            tmps_rimltc = {"rimltc_mr": rimltc_mr}

            self.ice4_rimltc(**state_rimltc, **tmps_rimltc)

            ####################### ice4_rimltc_post_processing #####################

            state_rimltc_pp = {
                **{
                    key: state[key]
                    for key in [
                        "t",
                        "exn",
                        "lv_fact",
                        "ls_fact",
                        "th_t",
                        "rc_t",
                        "ri_t",
                    ]
                },
            }

            self.ice4_rimltc_post_processing(**state_rimltc_pp, **tmps_rimltc)

            ######################## ice4_increment_update ##########################
            state_increment_update = {
                **{key: state[key] for key in ["ls_fact", "lv_fact"]},
                **{
                    key: state[key]
                    for key in [
                        "theta_increment",
                        "rv_increment",
                        "rc_increment",
                        "rr_increment",
                        "ri_increment",
                        "rs_increment",
                        "rg_increment",
                    ]
                },  # PB in F90
            }

            tmps_increment_update = {
                "rvheni_mr": rvheni_mr,
                "rimltc_mr": rimltc_mr,
                "rrhong_mr": rrhong_mr,
                "rsrimcg_mr": rsrimcg_mr,
            }

            self.ice4_increment_update(
                **state_increment_update, **tmps_increment_update
            )

            ######################## ice4_compute_pdf ###############################
            state_compute_pdf = {
                key: state[key]
                for key in [
                    "ldcompute",
                    "rhodref",
                    "rc_t",
                    "ri_t",
                    "cf",
                    "t",
                    "sigma_rc",
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
            }

            self.ice4_compute_pdf(**state_compute_pdf)

            # l263 to l278 omitted because LLRFR is False in AROME

            ######################## ice4_derived_fields ############################
            state_derived_fields = {
                key: state[key]
                for key in [
                    "t",
                    "rhodref",
                    "rv_t",
                    "pres",
                    "ssi",
                    "ka",
                    "dv",
                    "ai",
                    "cj",
                ]
            }

            self.ice4_derived_fields(**state_derived_fields)

            ######################## ice4_slope_parameters ##########################
            state_slope_parameters = {
                **{key: state[key] for key in ["rhodref", "t", "rr_t", "rs_t", "rg_t"]},
            }

            tmps_slopes = {
                "lbdar": lbdar,
                "lbdar_rf": lbdar_rf,
                "lbdas": lbdas,
                "lbdag": lbdag,
            }

            self.ice4_slope_parameters(**state_slope_parameters, **tmps_slopes)

            ######################## ice4_slow ######################################
            state_slow = {
                "ldcompute": state["ldcompute"],
                **{
                    key: state[key]
                    for key in [
                        "rhodref",
                        "t",
                        "ssi",
                        "lv_fact",
                        "ls_fact",
                        "rv_t",
                        "rc_t",
                        "ri_t",
                        "rs_t",
                        "rg_t",
                        "ai",
                        "cj",
                        "hli_hcf",
                        "hli_hri",
                    ]
                },
            }

            tmps_slow = {
                "lbdas": lbdas,
                "lbdag": lbdag,
                "rc_honi_tnd": rc_honi_tnd,
                "rv_deps_tnd": rv_deps_tnd,
                "ri_aggs_tnd": ri_aggs_tnd,
                "ri_auts_tnd": ri_auts_tnd,
                "rv_depg_tnd": rv_depg_tnd,
            }

            # Choice between processes computation (heavy) and light update
            self.ice4_slow_light(
                ldsoft=ldsoft, **state_slow, **tmps_slow
            ) if ldsoft else self.ice4_slow_heavy(
                ldsoft=ldsoft, **state_slow, **tmps_slow
            )

            ######################## ice4_warm ######################################
            state_warm = {
                "ldcompute": state["ldcompute"],
                **{
                    key: state[key]
                    for key in [
                        "rhodref",
                        "lv_fact",
                        "t",  # temperature
                        "th_t",
                        "pres",
                        "ka",  # thermal conductivity of the air
                        "dv",  # diffusivity of water vapour
                        "cj",  # function to compute the ventilation coefficient
                        "hlc_hcf",  # High Cloud Fraction in grid
                        "hlc_lcf",  # Low Cloud Fraction in grid
                        "hlc_hrc",  # LWC that is high in grid
                        "hlc_lrc",  # LWC that is low in grid
                        "rv_t",  # water vapour mixing ratio at t
                        "rc_t",  # cloud water mixing ratio at t
                        "rr_t",  # rain water mixing ratio at t
                        "cf",
                        "rf",
                    ]
                },
            }

            tmps_warm = {
                "lbdar": lbdar,
                "lbdar_rf": lbdar_rf,
                "rcautr": rcautr,
                "rcaccr": rcaccr,
                "rrevav": rrevav,
            }

            self.ice4_warm_light(
                **state_warm, **tmps_warm
            ) if ldsoft else self.ice4_warm_heavy(**state_warm, **tmps_warm)

            ######################## ice4_fast_rs ###################################
            state_fast_rs = {
                "ldcompute": state["ldcompute"],
                **{
                    key: state[key]
                    for key in [
                        "rhodref",
                        "lv_fact",
                        "ls_fact",
                        "pres",  # absolute pressure at t
                        "dv",  # diffusivity of water vapor in the air
                        "ka",  # thermal conductivity of the air
                        "cj",  # function to compute the ventilation coefficient
                        "t",
                        "rv_t",
                        "rc_t",
                        "rr_t",
                        "rs_t",
                    ]
                },
            }

            temporaries_fast_rs = {
                "lbdar": lbdar,
                "lbdar_rf": lbdar_rf,
                "rs_mltg_tnd": rs_mltg_tnd,
                "rc_mltsr_tnd": rc_mltsr_tnd,
                "rs_rcrims_tnd": rs_rcrims_tnd,  # extra dimension 8 in Fortran PRS_TEND
                "rs_rcrimss_tnd": rs_rcrimss_tnd,
                "rs_rsrimcg_tnd": rs_rsrimcg_tnd,
                "rs_rraccs_tnd": rs_rraccs_tnd,
                "rs_rraccss_tnd": rs_rraccss_tnd,
                "rs_rsaccrg_tnd": rs_rsaccrg_tnd,
                "rs_freez1_tnd": rs_freez1_tnd,
                "rs_freez2_tnd": rs_freez2_tnd,
                "riaggs": riaggs,
                "rcrimss": rcrimss,
                "rcrimsg": rcrimsg,
                "rsrimcg": rsrimcg,
                "rraccss": rraccss,
                "rraccsg": rraccsg,
                "rsaccrg": rsaccrg,
                "index_floor": index_floor,
                "index_floor_r": index_floor_r,
                "index_floor_s": index_floor_s,
            }

            gaminc_rim1 = from_array(
                self.gaminc_rim1, backend=self.gt4py_config.backend
            )
            gaminc_rim2 = from_array(
                self.gaminc_rim2, backend=self.gt4py_config.backend
            )
            gaminc_rim4 = from_array(
                self.gaminc_rim4, backend=self.gt4py_config.backend
            )

            ker_raccs = from_array(self.ker_raccs, backend=self.gt4py_config.backend)
            ker_raccss = from_array(
                self.ker_raccss, backend=self.gt4py_config.backendEND
            )
            ker_saccrg = from_array(self.ker_saccrg, backend=self.gt4py_config.backend)

            self.ice4_fast_rs_light(
                gaminc_rim1=gaminc_rim1,
                gaminc_rim2=gaminc_rim2,
                gaminc_rim4=gaminc_rim4,
                ker_raccs=ker_raccs,
                ker_raccss=ker_raccss,
                ker_saccrg=ker_saccrg,
                **state_fast_rs,
                **temporaries_fast_rs,
            ) if ldsoft else self.ice4_fast_rs_heavy(
                gaminc_rim1=gaminc_rim1,
                gaminc_rim2=gaminc_rim2,
                gaminc_rim4=gaminc_rim4,
                ker_raccs=ker_raccs,
                ker_raccss=ker_raccss,
                ker_saccrg=ker_saccrg,
                **state_fast_rs,
                **temporaries_fast_rs,
            )

            ######################## ice4_fast_rg_pre_processing ####################
            state_fast_rg_pp = {
                **{
                    key: state[key]
                    for key in [
                        "rgsi",
                        "rvdepg",
                        "rsmltg",
                        "rraccsg",
                        "rsaccrg",
                        "rcrimsg",
                        "rsrimcg",
                    ]
                },
            }

            tmps_fast_rg_pp = {
                "rgsi_mr": rgsi_mr,
                "rrhong_mr": rrhong_mr,
                "rsrimcg_mr": rsrimcg_mr,
            }

            self.ice4_fast_rg_pre_processing(**state_fast_rg_pp, **tmps_fast_rg_pp)

            ######################## ice4_fast_rg ###################################
            state_fast_rg = {
                "ldcompute": state["ldcompute"],
                **{
                    key: state[key]
                    for key in [
                        "t",
                        "rhodref",
                        "pres",
                        "rv_t",
                        "rr_t",
                        "ri_t",
                        "rg_t",
                        "rc_t",
                        "rs_t",
                        "ci_t",
                        "ka",
                        "dv",
                        "cj",
                    ]
                },
            }

            temporaries_fast_rg = {
                "lbdar": lbdar,
                "lbdas": lbdas,
                "lbdag": lbdag,
                "rg_rcdry_tnd": rg_rcdry_tnd,
                "rg_ridry_tnd": rg_ridry_tnd,
                "rg_rsdry_tnd": rg_rsdry_tnd,
                "rg_rrdry_tnd": rg_rrdry_tnd,
                "rg_riwet_tnd": rg_riwet_tnd,
                "rg_rswet_tnd": rg_rswet_tnd,
                "rg_freez1_tnd": rg_freez1_tnd,
                "rg_freez2_tnd": rg_freez2_tnd,
                "ricfrrg": ricfrrg,
                "rrcfrig": rrcfrig,
                "ricfrr": ricfrr,
                "rgmltr": rgmltr,
                "index_floor_s": index_floor_s,
                "index_floor_g": index_floor_g,
                "index_floor_r": index_floor_r,
            }

            ker_sdryg = from_array(self.ker_sdryg, backend=self.gt4py_config.backend)
            ker_rdryg = from_array(self.ker_rdryg, backend=self.gt4py_config.backend)

            self.ice4_fast_rg_light(
                ker_sdryg=ker_sdryg,
                ker_rdryg=ker_rdryg,
                **state_fast_rg,
                **temporaries_fast_rg,
            ) if ldsoft else self.ice4_fast_rg_heavy(
                ker_sdryg=ker_sdryg,
                ker_rdryg=ker_rdryg,
                **state_fast_rg,
                **temporaries_fast_rg,
            )

            ######################## ice4_fast_ri ###################################
            state_fast_ri = {
                "ldcompute": state["ldcompute"],
                **{
                    key: state[key]
                    for key in [
                        "rhodref",
                        "lv_fact",
                        "ls_fact",
                        "ai",
                        "cj",
                        "ci_t",
                        "ssi",
                        "rc_t",
                        "ri_t",
                    ]
                },
            }

            tmps_fast_ri = {
                "rc_beri_tnd": rc_beri_tnd,
            }

            self.ice4_fast_ri_light(
                **state_fast_ri, **tmps_fast_ri
            ) if ldsoft else self.ice4_fast_ri_heavy(**state_fast_ri, **tmps_fast_ri)

            ######################## ice4_tendencies_update #########################

            state_tendencies_update = {
                **{key: state[key] for key in ["ls_fact", "lv_fact"]},
                **{
                    key: state[key]
                    for key in [
                        "theta_tnd",
                        "rv_tnd",
                        "rc_tnd",
                        "rr_tnd",
                        "ri_tnd",
                        "rs_tnd",
                        "rg_tnd",
                    ]
                },
            }

            tmps_tnd_update = {
                "rvheni_mr": rvheni_mr,
                "rrhong_mr": rrhong_mr,
                "rimltc_mr": rimltc_mr,
                "rsrimcg_mr": rsrimcg_mr,
                "rchoni": rchoni,
                "rvdeps": rvdeps,
                "riaggs": riaggs,
                "riauts": riauts,
                "rvdepg": rvdepg,
                "rcautr": rcautr,
                "rcaccr": rcaccr,
                "rrevav": rrevav,
                "rcberi": rcberi,
                "rsmltg": rsmltg,
                "rcmltsr": rcmltsr,
                "rraccss": rraccss,  # 13
                "rraccsg": rraccsg,  # 14
                "rsaccrg": rsaccrg,  # 15  # Rain accretion onto the aggregates
                "rcrimss": rcrimss,  # 16
                "rcrimsg": rcrimsg,  # 17
                "rsrimcg": rsrimcg,  # 18  # Cloud droplet riming of the aggregates
                "ricfrrg": ricfrrg,  # 19
                "rrcfrig": rrcfrig,  # 20
                "ricfrr": ricfrr,  # 21  # Rain contact freezing
                "rcwetg": rcwetg,  # 22
                "riwetg": riwetg,  # 23
                "rrwetg": rrwetg,  # 24
                "rswetg": rswetg,  # 25  # Graupel wet growth
                "rcdryg": rcdryg,  # 26
                "ridryg": ridryg,  # 27
                "rrdryg": rrdryg,  # 28
                "rsdryg": rsdryg,  # 29  # Graupel dry growth
                "rgmltr": rgmltr,
            }

            self.ice4_tendencies_update(**state_tendencies_update, **tmps_tnd_update)
