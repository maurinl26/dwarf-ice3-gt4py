# -*- coding: utf-8 -*-
from ifs_physics_common.framework.stencil import compile_stencil
import sys
import logging

logging.basicConfig(stream=sys.stdout, level=logging.INFO)

STENCIL_COLLECTIONS = [
    "aro_filter",
    "ice_adjust",
    "ice4_nucleation",
    "nucleation_post_processing",
    "ice4_rimltc",
    "rimltc_post_processing",
    "ice4_increment_update",
    "ice4_derived_fields",
    "ice4_slope_parameters",
    "ice4_slow",
    "ice4_warm",
    "ice4_fast_rs",
    "ice4_fast_rg_pre_processing",
    "ice4_fast_rg",
    "ice4_fast_ri",
    "ice4_tendencies_update",
    "ice4_stepping_tmicro_init",
    "ice4_stepping_tsoft_init",
    "ice4_stepping_heat",
    "step_limiter",
    "mixing_ratio_step_limiter",
    "state_update",
    "external_tendencies_update",
]


def build(phyex, backend, config, stencil_collection):
    """Compile stencils given externals, backends and config

    Args:
        externals (_type_): list of externals
        backend (_type_): backend for gt4py
        config (_type_): gt4py config (including precision)
        stencil_collection (_type_): stencil_collection to compile
    """
    try:
        _ = compile_stencil(stencil_collection, config, externals=phyex.to_externals())
    except:
        logging.info(
            f"Building {stencil_collection}, on {backend} : Compilation Failed"
        )

    else:
        logging.info(
            f"Building {stencil_collection}, on {backend} : Compilation Succeded \n"
        )
