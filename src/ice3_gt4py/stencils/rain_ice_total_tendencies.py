# -*- coding: utf-8 -*-
from __future__ import annotations

from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import exp, log, computation, interval, PARALLEL
from ifs_physics_common.framework.stencil import stencil_collection
from ifs_physics_common.utils.f2py import ported_method


@ported_method(
    from_file="PHYEX/src/common/micro/rain_ice.F90", from_line=693, to_line=728
)
@stencil_collection("rain_ice_total_tendencies")
def rain_ice_total_tendencies(
    wr_th: gtscript.Field["float"],
    wr_v: gtscript.Field["float"],
    wr_c: gtscript.Field["float"],
    wr_r: gtscript.Field["float"],
    wr_i: gtscript.Field["float"],
    wr_s: gtscript.Field["float"],
    wr_g: gtscript.Field["float"],
    ls_fact: gtscript.Field["float"],
    lv_fact: gtscript.Field["float"],
    exnref: gtscript.Field["float"],
    ths: gtscript.Field["float"],
    rvs: gtscript.Field["float"],
    rcs: gtscript.Field["float"],
    rrs: gtscript.Field["float"],
    ris: gtscript.Field["float"],
    rss: gtscript.Field["float"],
    rgs: gtscript.Field["float"],
    rvheni: gtscript.Field["float"],
    rv_t: gtscript.Field["float"],
    rc_t: gtscript.Field["float"],
    rr_t: gtscript.Field["float"],
    ri_t: gtscript.Field["float"],
    rs_t: gtscript.Field["float"],
    rg_t: gtscript.Field["float"],
):
    """Update tendencies

    Args:
        wr_th (gtscript.Field[float]): potential temperature initial value
        wr_v (gtscript.Field[float]): vapour initial value
        wr_c (gtscript.Field[float]): cloud droplets initial value
        wr_r (gtscript.Field[float]): rain initial value
        wr_i (gtscript.Field[float]): ice initial value
        wr_s (gtscript.Field[float]): snow initial value
        wr_g (gtscript.Field[float]): graupel initial value
        ls_fact (gtscript.Field[float]): sublimation latent heat over heat capacity
        lv_fact (gtscript.Field[float]): vapourisation latent heat over heat capacity
        exnref (gtscript.Field[float]): reference exner pressure
        ths (gtscript.Field[float]): source (tendency) of potential temperature
        rvs (gtscript.Field[float]): source (tendency) of vapour
        rcs (gtscript.Field[float]): source (tendency) of cloud droplets
        rrs (gtscript.Field[float]): source (tendency) of rain
        ris (gtscript.Field[float]): source (tendency) of ice
        rss (gtscript.Field[float]): source (tendency) of snow
        rgs (gtscript.Field[float]): source (tendency) of graupel
        rvheni (gtscript.Field[float]): _description_
        rv_t (gtscript.Field[float]): vapour m.r. at t
        rc_t (gtscript.Field[float]): droplets m.r. at t
        rr_t (gtscript.Field[float]): rain m.r. at t
        ri_t (gtscript.Field[float]): ice m.r. at t
        rs_t (gtscript.Field[float]): snow m.r. at t
        rg_t (gtscript.Field[float]): graupel m.r. at t
    """

    from __externals__ import INV_TSTEP

    with computation(PARALLEL), interval(...):

        # Translation note ls, lv replaced by ls_fact, lv_fact

        # Hydrometeor tendency
        wr_v = (wr_v - rv_t) * INV_TSTEP
        wr_c = (wr_c - rc_t) * INV_TSTEP
        wr_r = (wr_r - rr_t) * INV_TSTEP
        wr_i = (wr_i - ri_t) * INV_TSTEP
        wr_s = (wr_s - rs_t) * INV_TSTEP
        wr_g = (wr_g - rg_t) * INV_TSTEP

        # Theta tendency
        wr_th = (wr_c + wr_r) * lv_fact + (wr_i + wr_s + wr_g) * ls_fact

        # Tendencies to sources, taking nucleation into account (rv_heni)
        ths += wr_th + rvheni * ls_fact
        rvs += wr_v - rvheni
        rcs += wr_c
        rrs += wr_r
        ris += wr_i + rvheni
        rss += wr_s
        rgs += wr_g
