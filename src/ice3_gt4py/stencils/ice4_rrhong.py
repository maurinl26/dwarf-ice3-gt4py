# -*- coding: utf-8 -*-
from __future__ import annotations

from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import interval, computation, PARALLEL
from ifs_physics_common.framework.stencil import stencil_collection
from ifs_physics_common.utils.f2py import ported_method


@ported_method(from_file="PHYEX/src/common/micro/mode_ice4_rrhong.F90")
@stencil_collection("ice4_rrhong")
def ice4_rrhong(
    ldcompute: gtscript.Field["bool"],
    t: gtscript.Field["float"],
    exn: gtscript.Field["float"],
    lv_fact: gtscript.Field["float"],
    ls_fact: gtscript.Field["float"],
    tht: gtscript.Field["float"],  # theta at time t
    rrhong_mr: gtscript.Field["float"],
    rr_t: gtscript.Field["float"],  # rain water mixing ratio at t
):
    """Compute the spontaneous frezzing source RRHONG

    Args:
        ldcompute (gtscript.Field[bool]): switch to activate microphysical processes on column
        t (gtscript.Field[float]): temperature at t
        exn (gtscript.Field[float]): exner pressure
        lv_fact (gtscript.Field[float]): vaporisation latent heat
        ls_fact (gtscript.Field[float]): sublimation latent heat
        tht (gtscript.Field[float]): potential temperature
        rr_t (gtscript.Field[float]): rain mixing ratio at t
        rrhong_mr (gtscript.Field[float]): mixing ratio for spontaneous freezing source
    """

    from __externals__ import LFEEDBACKT, R_RTMIN, TT

    # 3.3 compute the spontaneous frezzing source: RRHONG
    with computation(PARALLEL), interval(...):
        if t < TT - 35 and rr_t > R_RTMIN and ldcompute:
            rrhong_mr = rr_t

            # limitation for -35 degrees crossing
            if LFEEDBACKT == 1:
                rrhong_mr = min(
                    rrhong_mr, max(0, ((TT - 35) / exn - tht) / (ls_fact - lv_fact))
                )

        else:
            rrhong_mr = 0


@ported_method(
    from_file="PHYEX/src/common/micro/mode_ice4_tendencies.F90",
    from_line=166,
    to_line=171,
)
@stencil_collection("ice4_rrhong_post_processing")
def ice4_rrhong_post_processing(
    t: gtscript.Field["float"],
    exn: gtscript.Field["float"],
    ls_fact: gtscript.Field["float"],
    lv_fact: gtscript.Field["float"],
    tht: gtscript.Field["float"],
    rr_t: gtscript.Field["float"],
    rg_t: gtscript.Field["float"],
    rrhong_mr: gtscript.Field["float"],
):
    """adjust mixing ratio with nucleation increments

    Args:
        t (gtscript.Field[float]): temperature
        exn (gtscript.Field[float]): exner pressure
        ls_fact (gtscript.Field[float]): sublimation latent heat over heat capacity
        tht (gtscript.Field[float]): potential temperature
        rr_t (gtscript.Field[float]): rain m.r.
        rg_t (gtscript.Field[float]): graupel m.r.
        rrhong (gtscript.Field[float]): rain m.r. increment due to homogeneous nucleation
    """

    with computation(PARALLEL), interval(...):
        tht += rrhong_mr * (ls_fact - lv_fact)
        t = tht / exn
        rr_t -= rrhong_mr
        rg_t += rrhong_mr
