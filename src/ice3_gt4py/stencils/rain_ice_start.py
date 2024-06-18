# -*- coding: utf-8 -*-
from __future__ import annotations

from gt4py.cartesian.gtscript import (
    Field,
    computation,
    __INLINED,
    PARALLEL,
    interval,
    __externals__,
)
from ifs_physics_common.framework.stencil import stencil_collection
from ifs_physics_common.utils.f2py import ported_method


@ported_method(
    from_file="PHYEX/src/common/micro/rain_ice.F90", from_line=367, to_line=396
)
@stencil_collection("rain_ice_init")
def rain_ice_init(
    ldmicro: Field["bool"],
    exn: Field["float"],
    ls_fact: Field["float"],
    lv_fact: Field["float"],
    th_t: Field["float"],
    rv_t: Field["float"],
    rc_t: Field["float"],
    rr_t: Field["float"],
    ri_t: Field["float"],
    rs_t: Field["float"],
    rg_t: Field["float"],
):
    """Computes ldmicro mask given minimum mixing ratios per specy.
    Computes conversions from theta, exn to temperature

    Args:
        ldmicro (Field[bool]): mask for microphysical computations
        exn (Field[float]): exner pressure
        ls_fact (Field[float]): sublimation latent heat over capacity
        lv_fact (Field[float]): vaporisation latent heat over capacity
        th_t (Field[float]): potential temperature at t
        rv_t (Field[float]): vapour m.r. at t
        rc_t (Field[float]): cloud droplet m.r. at t
        rr_t (Field[float]): rain m.r. at t
        ri_t (Field[float]): ice m.r. at t
        rs_t (Field[float]): snow m.r.
        rg_t (Field[float]): graupel m.r.
    """

    from __externals__ import (
        C_RTMIN,
        R_RTMIN,
        I_RTMIN,
        S_RTMIN,
        G_RTMIN,
        CPD,
        CPV,
        CI,
        CL,
        TT,
        LSTT,
        LVTT,
    )

    with computation(PARALLEL), interval(...):
        divider = CPD + CPV * rv_t + CL * (rc_t + rr_t) + CI * (ri_t + rs_t + rg_t)
        t = th_t * exn
        ls_fact = (LSTT + (CPV - CI) * (t - TT)) / divider
        lv_fact = (LVTT + (CPV - CL) * (t - TT)) / divider

        ldmicro = (
            rc_t > C_RTMIN
            or rr_t > R_RTMIN
            or ri_t > I_RTMIN
            or rs_t > S_RTMIN
            or rg_t > G_RTMIN
        )


@ported_method(
    from_file="PHYEX/src/common/micro/rain_ice.F90", from_line=424, to_line=444
)
@stencil_collection("initial_values_saving")
def initial_values_saving(
    wr_th: Field["float"],
    wr_v: Field["float"],
    wr_c: Field["float"],
    wr_r: Field["float"],
    wr_i: Field["float"],
    wr_s: Field["float"],
    wr_g: Field["float"],
    th_t: Field["float"],
    rv_t: Field["float"],
    rc_t: Field["float"],
    rr_t: Field["float"],
    ri_t: Field["float"],
    rs_t: Field["float"],
    rg_t: Field["float"],
    evap3d: Field["float"],
    rainfr: Field["float"],
):

    from __externals__ import LWARM

    with computation(PARALLEL), interval(...):
        wr_th = th_t
        wr_v = rv_t
        wr_c = rc_t
        wr_r = rr_t
        wr_i = ri_t
        wr_s = rs_t
        wr_g = rg_t

        # LWARM is True for AROME
        if __INLINED(LWARM):
            evap3d = 0
        rainfr = 0
