# -*- coding: utf-8 -*-
from __future__ import annotations

from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import __INLINED, computation, interval, PARALLEL
from ifs_physics_common.framework.stencil import stencil_collection
from ifs_physics_common.utils.f2py import ported_method


@ported_method(from_file="PHYEX/src/common/micro/mode_ice4_fast_ri.F90")
@stencil_collection("ice4_fast_ri")
def ice4_fast_ri(
    ldcompute: gtscript.Field["bool"],
    rhodref: gtscript.Field["float"],
    lv_fact: gtscript.Field["float"],
    ls_fact: gtscript.Field["float"],
    ai: gtscript.Field["float"],
    cj: gtscript.Field["float"],
    ci_t: gtscript.Field["float"],
    ssi: gtscript.Field["float"],
    rc_t: gtscript.Field["float"],
    ri_t: gtscript.Field["float"],
    rc_beri_tnd: gtscript.Field["float"],
):
    """Computes Bergeron-Findeisen effect RCBERI.

    Evaporation of cloud droplets for deposition over ice-crystals.

    Args:
        lcompute (gtscript.Field[bool]): switch to compute microphysical processes
        lv_fact (gtscript.Field[float]): latent heat of vaporisation
        ls_fact (gtscript.Field[float]): latent heat of sublimation
        ai (gtscript.Field[float]): thermodynamical function
        cj (gtscript.Field[float]): function to compute ventilation factor
        ci_t (gtscript.Field[float]): concentration of ice at t
        ssi (gtscript.Field[float]): supersaturation over ice
        rc_t (gtscript.Field[float]): cloud droplets mixing ratio at t
        ri_t (gtscript.Field[float]): pristine ice mixing ratio at t
        rc_beri_tnd (gtscript.Field[float]): tendency for Bergeron Findeisen effect
    """

    from __externals__ import LDSOFT, C_RTMIN, DI, I_RTMIN, LBEXI, LBI, O0DEPI, O2DEPI

    # 7.2 Bergeron-Findeisen effect: RCBERI
    with computation(PARALLEL), interval(...):
        if __INLINED(not LDSOFT):
            if (
                ssi > 0
                and rc_t > C_RTMIN
                and ri_t > I_RTMIN
                and ci_t > 1e-20
                and ldcompute
            ):
                rc_beri_tnd = min(
                    1e-8, LBI * (rhodref * ri_t / ci_t) ** LBEXI
                )  # lambda_i
                rc_beri_tnd = (
                    (ssi / (rhodref * ai))
                    * ci_t
                    * (O0DEPI / rc_beri_tnd + O2DEPI * cj / rc_beri_tnd ** (DI + 2.0))
                )

            else:
                rc_beri_tnd = 0
