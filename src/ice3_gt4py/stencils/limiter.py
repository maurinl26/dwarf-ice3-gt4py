# -*- coding: utf-8 -*-
from __future__ import annotations

from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import computation, PARALLEL, interval
from ifs_physics_common.framework.stencil import stencil_collection
from ifs_physics_common.utils.f2py import ported_method
from ice3_gt4py.functions.sign import sign


@ported_method(
    from_file="PHYEX/src/common/micro/mode_ice4_stepping.F90",
    from_line=290,
    to_line=388,
)
@stencil_collection("limiter")
def limiter(
    exn: gtscript.Field["float"],
    theta_t: gtscript.Field["float"],
    theta_a_tnd: gtscript.Field["float"],
    theta_b: gtscript.Field["float"],
    theta_ext_tnd: gtscript.Field["float"],
    rc_t: gtscript.Field["float"],
    rr_t: gtscript.Field["float"],
    ri_t: gtscript.Field["float"],
    rs_t: gtscript.Field["float"],
    rg_t: gtscript.Field["float"],
    rc_0r_t: gtscript.Field["float"],
    rr_0r_t: gtscript.Field["float"],
    ri_0r_t: gtscript.Field["float"],
    rs_0r_t: gtscript.Field["float"],
    rg_0r_t: gtscript.Field["float"],
    rc_a_tnd: gtscript.Field["float"],
    rr_a_tnd: gtscript.Field["float"],
    ri_a_tnd: gtscript.Field["float"],
    rs_a_tnd: gtscript.Field["float"],
    rg_a_tnd: gtscript.Field["float"],
    rc_ext_tnd: gtscript.Field["float"],
    rr_ext_tnd: gtscript.Field["float"],
    ri_ext_tnd: gtscript.Field["float"],
    rs_ext_tnd: gtscript.Field["float"],
    rg_ext_tnd: gtscript.Field["float"],
    rc_b: gtscript.Field["float"],
    rr_b: gtscript.Field["float"],
    ri_b: gtscript.Field["float"],
    rs_b: gtscript.Field["float"],
    rg_b: gtscript.Field["float"],
    delta_t_micro: gtscript.Field["float"],
    t_micro: gtscript.Field["float"],
    delta_t_soft: gtscript.Field["float"],
    t_soft: gtscript.Field["float"],
    ldcompute: gtscript.Field["bool"],
):
    from __externals__ import (
        C_RTMIN,
        G_RTMIN,
        I_RTMIN,
        MNH_TINY,
        R_RTMIN,
        S_RTMIN,
        TSTEP,
        TSTEP_TS,
        TT,
        MRSTEP,
    )

    # Adding externals tendencies
    # TODO : clear into another stencil
    with computation(PARALLEL), interval(...):
        theta_a_tnd += theta_ext_tnd
        rc_a_tnd += rc_ext_tnd
        rr_a_tnd += rr_ext_tnd
        ri_a_tnd += ri_ext_tnd
        rs_a_tnd += rs_ext_tnd
        rg_a_tnd += rg_ext_tnd

    # 4.6 Time integration
    with computation(PARALLEL), interval(...):
        delta_t_micro = TSTEP - t_micro if ldcompute else 0

    # Adjustment of tendencies when temperature reaches 0
    with computation(PARALLEL), interval(...):
        theta_tt = TT / exn
        if (theta_t - theta_tt) * (theta_t + theta_b - theta_tt) < 0:
            delta_t_micro = 0

        if abs(theta_a_tnd > 1e-20):
            delta_t_tmp = (theta_tt - theta_b - theta_t) / theta_a_tnd
            if delta_t_tmp > 0:
                delta_t_micro = min(delta_t_micro, delta_t_tmp)

    # Tendencies adjustment if a speci disappears
    # (c)
    with computation(PARALLEL), interval(...):
        delta_t_micro = mixing_ratio_step_limiter(
            rc_a_tnd, rc_b, rc_t, delta_t_micro, C_RTMIN, MNH_TINY
        )
    # (r)
    with computation(PARALLEL), interval(...):
        delta_t_micro = mixing_ratio_step_limiter(
            rr_a_tnd, rr_b, rr_t, delta_t_micro, R_RTMIN, MNH_TINY
        )
    # (i)
    with computation(PARALLEL), interval(...):
        delta_t_micro = mixing_ratio_step_limiter(
            ri_a_tnd, ri_b, ri_t, delta_t_micro, I_RTMIN, MNH_TINY
        )
    # (s)
    with computation(PARALLEL), interval(...):
        delta_t_micro = mixing_ratio_step_limiter(
            rs_a_tnd, rs_b, rs_t, delta_t_micro, S_RTMIN, MNH_TINY
        )
    # (g)
    with computation(PARALLEL), interval(...):
        delta_t_micro = mixing_ratio_step_limiter(
            rg_a_tnd, rg_b, rg_t, delta_t_micro, G_RTMIN, MNH_TINY
        )

    # We stop when the end of the timestep is reached
    with computation(PARALLEL), interval(...):
        ldcompute = False if t_micro + delta_t_micro > TSTEP else ldcompute

    ################## Mixing ratio limiter ###################
    # TODO : TSTEP_TS out of the loop
    with computation(PARALLEL), interval(...):
        if TSTEP_TS != 0:
            if t_micro + delta_t_micro > t_soft + delta_t_soft:
                delta_t_micro = t_soft + delta_t_soft - t_micro
                ldcompute = False

    with computation(PARALLEL), interval(...):
        # TODO: add condition on LL_ANY_ITER
        time_threshold_tmp = (
            (sign(1, rc_a_tnd) * MRSTEP + rc_0r_t - rc_t - rc_b)
            if abs(rc_a_tnd) > 1e-20
            else -1
        )

    # l363
    with computation(PARALLEL), interval(...):
        if (
            time_threshold_tmp >= 0
            and time_threshold_tmp < delta_t_micro
            and (rc_t > C_RTMIN or rc_a_tnd > 0)
        ):
            delta_t_micro = min(delta_t_micro, time_threshold_tmp)
            ldcompute = False

            # Translation note : ldcompute is LLCOMPUTE in mode_ice4_stepping.F90

    # l370
    # Translation note : l370 to l378 in mode_ice4_stepping. F90 contracted in a single stencil
    with computation(PARALLEL), interval(...):
        r_b_max = abs(rr_b)

    ################ (r) #############
    with computation(PARALLEL), interval(...):
        time_threshold_tmp = (
            (sign(1, rr_a_tnd) * MRSTEP + rr_0r_t - rr_t - rr_b)
            if abs(rr_a_tnd) > 1e-20
            else -1
        )

    with computation(PARALLEL), interval(...):
        if (
            time_threshold_tmp >= 0
            and time_threshold_tmp < delta_t_micro
            and (rr_t > R_RTMIN or rr_a_tnd > 0)
        ):
            delta_t_micro = min(delta_t_micro, time_threshold_tmp)
            ldcompute = False

    with computation(PARALLEL), interval(...):
        r_b_max = max(r_b_max, abs(rr_b))

    ################ (i) #############
    with computation(PARALLEL), interval(...):
        time_threshold_tmp = (
            (sign(1, ri_a_tnd) * MRSTEP + ri_0r_t - ri_t - ri_b)
            if abs(ri_a_tnd) > 1e-20
            else -1
        )

    with computation(PARALLEL), interval(...):
        if (
            time_threshold_tmp >= 0
            and time_threshold_tmp < delta_t_micro
            and (rc_t > I_RTMIN or ri_a_tnd > 0)
        ):
            delta_t_micro = min(delta_t_micro, time_threshold_tmp)
            ldcompute = False

    with computation(PARALLEL), interval(...):
        r_b_max = max(r_b_max, abs(ri_b))

    ################ (s) #############
    with computation(PARALLEL), interval(...):
        time_threshold_tmp = (
            (sign(1, rs_a_tnd) * MRSTEP + rs_0r_t - rs_t - rs_b)
            if abs(rs_a_tnd) > 1e-20
            else -1
        )

    with computation(PARALLEL), interval(...):
        if (
            time_threshold_tmp >= 0
            and time_threshold_tmp < delta_t_micro
            and (rs_t > S_RTMIN or rs_a_tnd > 0)
        ):
            delta_t_micro = min(delta_t_micro, time_threshold_tmp)
            ldcompute = False

    with computation(PARALLEL), interval(...):
        r_b_max = max(r_b_max, abs(rs_b))

    ################ (g) #############
    with computation(PARALLEL), interval(...):
        time_threshold_tmp = (
            (sign(1, rg_a_tnd) * MRSTEP + rg_0r_t - rg_t - rg_b)
            if abs(rg_a_tnd) > 1e-20
            else -1
        )

    with computation(PARALLEL), interval(...):
        if (
            time_threshold_tmp >= 0
            and time_threshold_tmp < delta_t_micro
            and (rg_t > G_RTMIN or rg_a_tnd > 0)
        ):
            delta_t_micro = min(delta_t_micro, time_threshold_tmp)
            ldcompute = False

    with computation(PARALLEL), interval(...):
        r_b_max = max(r_b_max, abs(rg_b))  # (g)

    # Limiter on max mixing ratio
    with computation(PARALLEL), interval(...):
        if r_b_max > MRSTEP:
            delta_t_micro = 0
            ldcompute = False
