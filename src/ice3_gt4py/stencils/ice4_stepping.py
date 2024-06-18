# -*- coding: utf-8 -*-
from __future__ import annotations

from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import interval, computation, PARALLEL
from ifs_physics_common.framework.stencil import stencil_collection

from ifs_physics_common.utils.f2py import ported_method

from ice3_gt4py.functions.ice_adjust import (
    cph,
    sublimation_latent_heat,
    vaporisation_latent_heat,
)
from ice3_gt4py.functions.temperature import theta2temperature


@ported_method(
    from_file="PHYEX/src/common/micro/mode_ice4_stepping.F90",
    from_line=215,
    to_line=221,
)
@stencil_collection("ice4_stepping_tmicro_init")
def ice4_stepping_tmicro_init(
    t_micro: gtscript.Field["float"], ldmicro: gtscript.Field["bool"]
):
    """Initialise t_soft with value of t_micro after each loop
    on LSOFT condition.

    Args:
        t_micro (gtscript.Field[float]): time for microphsyics loops
        ldmicro (gtscript.Field[bool]): microphsyics activation mask
    """

    from __externals__ import TSTEP

    # 4.4 Temporal loop
    with computation(PARALLEL), interval(...):
        t_micro = 0 if ldmicro else TSTEP


@ported_method(
    from_file="PHYEX/src/common/micro/mode_ice4_stepping.F90",
    from_line=244,
    to_line=254,
)
@stencil_collection("ice4_stepping_heat")
def ice4_stepping_heat(
    rv_t: gtscript.Field["float"],
    rc_t: gtscript.Field["float"],
    rr_t: gtscript.Field["float"],
    ri_t: gtscript.Field["float"],
    rs_t: gtscript.Field["float"],
    rg_t: gtscript.Field["float"],
    exn: gtscript.Field["float"],
    th_t: gtscript.Field["float"],
    ls_fact: gtscript.Field["float"],
    lv_fact: gtscript.Field["float"],
    t: gtscript.Field["float"],
):
    """Compute and convert heat variables before computations

    Args:
        rv_t (gtscript.Field[float]): vapour mixing ratio
        rc_t (gtscript.Field[float]): cloud droplet mixing ratio
        rr_t (gtscript.Field[float]): rain m.r.
        ri_t (gtscript.Field[float]): ice m.r.
        rs_t (gtscript.Field[float]): snow m.r.
        rg_t (gtscript.Field[float]): graupel m.r.
        exn (gtscript.Field[float]): exner pressure
        th_t (gtscript.Field[float]): potential temperature
        ls_fact (gtscript.Field[float]): sublimation latent heat over heat capacity
        lv_fact (gtscript.Field[float]): vapourisation latent heat over heat capacity
        t (gtscript.Field[float]): temperature
    """
    with computation(PARALLEL), interval(...):
        specific_heat = cph(rv_t, rc_t, ri_t, rr_t, rs_t, rg_t)
        t = theta2temperature(th_t, exn)
        ls_fact = sublimation_latent_heat(t) / specific_heat
        lv_fact = vaporisation_latent_heat(t) / specific_heat


@ported_method(
    from_file="PHYEX/src/common/micro/mode_ice4_stepping.F90",
    from_line=230,
    to_line=237,
)
@stencil_collection("ice4_update_ldcompute")
def ice4_update_ldcompute(
    ldcompute: gtscript.Field["bool"], t_micro: gtscript.Field["float"]
):
    """Initialize or update ldcompute mask

    Args:
        ldcompute (gtscript.Field[bool]): mask to compute microphysical species
        t_micro (gtscript.Field[float]): microphysical time-step
    """

    from __externals__ import TSTEP

    with computation(PARALLEL), interval(...):
        ldcompute = t_micro < TSTEP
