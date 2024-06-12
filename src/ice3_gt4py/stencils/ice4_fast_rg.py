# -*- coding: utf-8 -*-
from __future__ import annotations

from gt4py.cartesian.gtscript import (
    __INLINED,
    exp,
    log,
    computation,
    interval,
    PARALLEL,
)
from gt4py.cartesian import gtscript
from ifs_physics_common.framework.stencil import stencil_collection
from ifs_physics_common.utils.f2py import ported_method

from ice3_gt4py.functions.interp_micro import (
    index_micro2d_dry_g,
    index_micro2d_dry_r,
    index_micro2d_dry_s,
)


@ported_method(from_file="PHYEX/src/common/micro/mode_ice4_fast_rg.F90")
@stencil_collection("ice4_fast_rg")
def ice4_fast_rg(
    ldcompute: gtscript.Field["bool"],
    t: gtscript.Field["float"],
    rhodref: gtscript.Field["float"],
    pres: gtscript.Field["float"],
    rv_t: gtscript.Field["float"],
    rr_t: gtscript.Field["float"],
    ri_t: gtscript.Field["float"],
    rg_t: gtscript.Field["float"],
    rc_t: gtscript.Field["float"],
    rs_t: gtscript.Field["float"],
    ci_t: gtscript.Field["float"],
    ka: gtscript.Field["float"],
    dv: gtscript.Field["float"],
    cj: gtscript.Field["float"],
    lbdar: gtscript.Field["float"],
    lbdas: gtscript.Field["float"],
    lbdag: gtscript.Field["float"],
    ricfrrg: gtscript.Field["float"],
    rrcfrig: gtscript.Field["float"],
    ricfrr: gtscript.Field["float"],
    rg_rcdry_tnd: gtscript.Field["float"],
    rg_ridry_tnd: gtscript.Field["float"],
    rg_rsdry_tnd: gtscript.Field["float"],
    rg_rrdry_tnd: gtscript.Field["float"],
    rg_riwet_tnd: gtscript.Field["float"],
    rg_rswet_tnd: gtscript.Field["float"],
    rg_freez1_tnd: gtscript.Field["float"],
    rg_freez2_tnd: gtscript.Field["float"],
    rgmltr: gtscript.Field["float"],
    ker_sdryg: gtscript.GlobalTable[float, (40, 40)],
    ker_rdryg: gtscript.GlobalTable[float, (40, 40)],
    index_floor_s: gtscript.Field["int"],
    index_floor_g: gtscript.Field["int"],
    index_floor_r: gtscript.Field["int"],
):
    """Compute fast graupel sources

    Args:
        ldcompute (gtscript.Field[int]): switch to compute microphysical processes on column
        t (gtscript.Field[float]): temperature
        rhodref (gtscript.Field[float]): reference density
        ri_t (gtscript.Field[float]): ice mixing ratio at t

        rg_t (gtscript.Field[float]): graupel m.r. at t
        rc_t (gtscript.Field[float]): cloud droplets m.r. at t
        rs_t (gtscript.Field[float]): snow m.r. at t
        ci_t (gtscript.Field[float]): _description_
        dv (gtscript.Field[float]): diffusivity of water vapor
        ka (gtscript.Field[float]): thermal conductivity of the air
        cj (gtscript.Field[float]): function to compute the ventilation coefficient
        lbdar (gtscript.Field[float]): slope parameter for rain
        lbdas (gtscript.Field[float]): slope parameter for snow
        lbdag (gtscript.Field[float]): slope parameter for graupel
        ricfrrg (gtscript.Field[float]): rain contact freezing
        rrcfrig (gtscript.Field[float]): rain contact freezing
        ricfrr (gtscript.Field[float]): rain contact freezing
        rg_rcdry_tnd (gtscript.Field[float]): Graupel wet growth
        rg_ridry_tnd (gtscript.Field[float]): Graupel wet growth
        rg_riwet_tnd (gtscript.Field[float]): Graupel wet growth
        rg_rsdry_tnd (gtscript.Field[float]): Graupel wet growth
        rg_rswet_tnd (gtscript.Field[float]): Graupel wet growth
        gdry (gtscript.Field[int]): boolean gtscript.Field
    """

    from __externals__ import (
        LDSOFT,
        LCRFLIMIT,
        ALPI,
        ALPW,
        BETAI,
        BETAW,
        BS,
        CEXVT,
        CI,
        CL,
        COLEXIG,
        COLIG,
        COLSG,
        CPV,
        CXG,
        CXS,
        DG,
        EPSILO,
        ESTT,
        EX0DEPG,
        EX1DEPG,
        EXICFRR,
        EXRCFRI,
        FCDRYG,
        FIDRYG,
        FRDRYG,
        FSDRYG,
        G_RTMIN,
        GAMI,
        GAMW,
        I_RTMIN,
        ICFRR,
        LBSDRYG1,
        LBSDRYG2,
        LBSDRYG3,
        LEVLIMIT,
        LMTT,
        LNULLWETG,
        LVTT,
        LWETGPOST,
        O0DEPG,
        O1DEPG,
        R_RTMIN,
        RCFRI,
        RV,
        S_RTMIN,
        TT,
    )

    # 6.1 rain contact freezing
    with computation(PARALLEL), interval(...):
        if ri_t > I_RTMIN and rr_t > R_RTMIN and ldcompute:
            # not LDSOFT : compute the tendencies
            if __INLINED(not LDSOFT):
                ricfrrg = ICFRR * ri_t * lbdar**EXICFRR * rhodref ** (-CEXVT)
                rrcfrig = RCFRI * ci_t * lbdar**EXRCFRI * rhodref ** (-CEXVT)

                if __INLINED(LCRFLIMIT):
                    zw0d = max(
                        0,
                        min(
                            1,
                            (ricfrrg * CI + rrcfrig * CL)
                            * (TT - t)
                            / max(1e-20, LVTT * rrcfrig),
                        ),
                    )
                    rrcfrig = zw0d * rrcfrig
                    ricfrr = (1 - zw0d) * rrcfrig
                    ricfrrg = zw0d * ricfrrg

                else:
                    ricfrr = 0

        else:
            ricfrrg = 0
            rrcfrig = 0
            ricfrr = 0

    # 6.3 compute graupel growth
    with computation(PARALLEL), interval(...):
        if rg_t > G_RTMIN and rc_t > R_RTMIN and ldcompute:
            if __INLINED(not LDSOFT):
                rg_rcdry_tnd = lbdag ** (CXG - DG - 2.0) * rhodref ** (-CEXVT)
                rg_rcdry_tnd = rg_rcdry_tnd * FCDRYG * rc_t

        else:
            rg_rcdry_tnd = 0

        if rg_t > G_RTMIN and ri_t > I_RTMIN and ldcompute:
            if __INLINED(not LDSOFT):
                rg_ridry_tnd = lbdag ** (CXG - DG - 2.0) * rhodref ** (-CEXVT)
                rg_ridry_tnd = FIDRYG * exp(COLEXIG * (t - TT)) * ri_t * rg_ridry_tnd
                rg_riwet_tnd = rg_ridry_tnd / (COLIG * exp(COLEXIG * (t - TT)))

        else:
            rg_ridry_tnd = 0
            rg_riwet_tnd = 0

    # 6.2.1 wet and dry collection of rs on graupel
    # Translation note : l171 in mode_ice4_fast_rg.F90
    with computation(PARALLEL), interval(...):
        if rs_t > S_RTMIN and rg_t > G_RTMIN and ldcompute:
            gdry = True  # GDRY is a boolean gtscript.Field in f90

        else:
            gdry = False
            rg_rsdry_tnd = 0
            rg_rswet_tnd = 0

    with computation(PARALLEL), interval(...):
        if __INLINED(not LDSOFT):
            if gdry:
                index_floor_s, index_float_s = index_micro2d_dry_s(lbdas)
                index_floor_g, index_float_g = index_micro2d_dry_g(lbdag)
                zw_tmp = index_float_g * (
                    index_float_s * ker_sdryg.A[index_floor_g + 1, index_floor_s + 1]
                    + (1 - index_float_s)
                    * ker_sdryg.A[index_floor_g + 1, index_floor_s]
                ) + (1 - index_float_g) * (
                    index_float_s * ker_sdryg.A[index_floor_g, index_floor_s + 1]
                    + (1 - index_float_s) * ker_sdryg.A[index_floor_g, index_floor_s]
                )

    with computation(PARALLEL), interval(...):
        # Translation note : #ifdef REPRO48 l192 to l198 kept
        #                                   l200 to l206 removed
        if gdry:
            rg_rswet_tnd = (
                FSDRYG
                * zw_tmp
                / COLSG
                * (lbdas * (CXS - BS))
                * (lbdag**CXG)
                * (rhodref ** (-CEXVT))
                * (
                    LBSDRYG1 / (lbdag**2)
                    + LBSDRYG2 / (lbdag * lbdas)
                    + LBSDRYG3 / (lbdas**2)
                )
            )

            rg_rsdry_tnd = rg_rswet_tnd * COLSG * exp(t - TT)

    # 6.2.6 accreation of raindrops on the graupeln
    with computation(PARALLEL), interval(...):
        if rr_t < R_RTMIN and rg_t < G_RTMIN and ldcompute:
            gdry = True
        else:
            gdry = False
            rg_rrdry_tnd = 0

    with computation(PARALLEL), interval(...):
        if __INLINED(not LDSOFT):
            index_floor_g, index_float_g = index_micro2d_dry_g(lbdag)
            index_floor_r, index_float_r = index_micro2d_dry_r(lbdar)
            zw_tmp = index_float_r * (
                index_float_g * ker_rdryg.A[index_floor_r + 1, index_floor_g + 1]
                + (1 - index_float_g) * ker_rdryg.A[index_floor_r + 1, index_floor_g]
            ) + (1 - index_float_r) * (
                index_float_g * ker_rdryg.A[index_floor_r, index_floor_g + 1]
                + (1 - index_float_g) * ker_rdryg.A[index_floor_r, index_floor_g]
            )

    # # l233
    with computation(PARALLEL), interval(...):
        if __INLINED(not LDSOFT):
            if gdry:
                rg_rrdry_tnd = (
                    FRDRYG
                    * zw_tmp
                    * (lbdar ** (-4))
                    * (lbdag**CXG)
                    * (rhodref ** (-CEXVT - 1))
                    * (
                        LBSDRYG1 / (lbdag**2)
                        + LBSDRYG2 / (lbdag * lbdar)
                        + LBSDRYG3 / (lbdar**2)
                    )
                )

    # l245
    with computation(PARALLEL), interval(...):
        rdryg_init_tmp = rg_rcdry_tnd + rg_ridry_tnd + rg_rsdry_tnd + rg_rrdry_tnd

    # Translation note l300 to l316 removed (no hail)

    # Freezing rate and growth mode
    # Translation note : l251 in mode_ice4_fast_rg.F90
    with computation(PARALLEL), interval(...):
        if rg_t > G_RTMIN and ldcompute:
            # Duplicated code with ice4_fast_rs
            if __INLINED(not LDSOFT):
                rg_freez1_tnd = rv_t * pres / (EPSILO + rv_t)
                if LEVLIMIT:
                    rg_freez1_tnd = min(
                        rg_freez1_tnd, exp(ALPI - BETAI / t - GAMI * log(t))
                    )

                rg_freez1_tnd = ka * (TT - t) + dv * (LVTT + (CPV - CL) * (t - TT)) * (
                    ESTT - rg_freez1_tnd
                ) / (RV * t)
                rg_freez1_tnd *= (
                    O0DEPG * lbdag**EX0DEPG + O1DEPG * cj * lbdag**EX1DEPG
                ) / (rhodref * (LMTT - CL * (TT - t)))
                rg_freez2_tnd = (rhodref * (LMTT + (CI - CL) * (TT - t))) / (
                    rhodref * (LMTT - CL * (TT - t))
                )

            rwetg_init_tmp = max(
                rg_riwet_tnd + rg_rswet_tnd,
                max(0, rg_freez1_tnd + rg_freez2_tnd * (rg_riwet_tnd + rg_rswet_tnd)),
            )

            # Growth mode
            # bool calculation :
            ldwetg = (
                1
                if (
                    max(0, rwetg_init_tmp - rg_riwet_tnd - rg_rswet_tnd)
                    <= max(0, rdryg_init_tmp - rg_ridry_tnd - rg_rsdry_tnd)
                )
                else 0
            )

            if __INLINED(not LNULLWETG):
                ldwetg = 1 if (ldwetg == 1 and rdryg_init_tmp > 0) else 0

            else:
                ldwetg = 1 if (ldwetg == 1 and rwetg_init_tmp > 0) else 0

            if __INLINED(not LWETGPOST):
                ldwetg = 1 if (ldwetg == 1 and t < TT) else 0

            lldryg = (
                1
                if (
                    t < TT
                    and rdryg_init_tmp > 1e-20
                    and max(0, rwetg_init_tmp - rg_riwet_tnd - rg_rswet_tnd)
                    > max(0, rg_rsdry_tnd - rg_ridry_tnd - rg_rsdry_tnd)
                )
                else 0
            )

        else:
            rg_freez1_tnd = 0
            rg_freez2_tnd = 0
            rwetg_init_tmp = 0
            ldwetg = 0
            lldryg = 0

    # l317
    with computation(PARALLEL), interval(...):
        if ldwetg == 1:
            rr_wetg = -(rg_riwet_tnd + rg_rswet_tnd + rg_rcdry_tnd - rwetg_init_tmp)
            rc_wetg = rg_rcdry_tnd
            ri_wetg = rg_riwet_tnd
            rs_wetg = rg_rswet_tnd

        else:
            rr_wetg = 0
            rc_wetg = 0
            ri_wetg = 0
            rs_wetg = 0

        if lldryg == 1:
            rc_dry = rg_rcdry_tnd
            rr_dry = rg_rrdry_tnd
            ri_dry = rg_ridry_tnd
            rs_dry = rg_rsdry_tnd

        else:
            rc_dry = 0
            rr_dry = 0
            ri_dry = 0
            rs_dry = 0

    # 6.5 Melting of the graupel
    with computation(PARALLEL), interval(...):
        if rg_t > G_RTMIN and t > TT and ldcompute:
            if __INLINED(not LDSOFT):
                rgmltr = rv_t * pres / (EPSILO + rv_t)
                if LEVLIMIT:
                    rgmltr = min(rgmltr, exp(ALPW - BETAW / t - GAMW * log(t)))

                rgmltr = ka * (TT - t) + dv * (LVTT + (CPV - CL) * (t - TT)) * (
                    ESTT - rgmltr
                ) / (RV * t)
                rgmltr = max(
                    0,
                    (
                        -rgmltr
                        * (O0DEPG * lbdag**EX0DEPG + O1DEPG * cj * lbdag**EX1DEPG)
                        - (rg_rcdry_tnd + rg_rrdry_tnd) * (rhodref * CL * (TT - t))
                    )
                    / (rhodref * LMTT),
                )

        else:
            rgmltr = 0


@ported_method(
    from_file="PHYEX/src/common/micro/mode_ice4_tendencies.F90",
    from_line=386,
    to_line=390,
)
@stencil_collection("ice4_fast_rg_pre_processing")
def ice4_fast_rg_pre_post_processing(
    rgsi: gtscript.Field["float"],
    rgsi_mr: gtscript.Field["float"],
    rvdepg: gtscript.Field["float"],
    rsmltg: gtscript.Field["float"],
    rraccsg: gtscript.Field["float"],
    rsaccrg: gtscript.Field["float"],
    rcrimsg: gtscript.Field["float"],
    rsrimcg: gtscript.Field["float"],
    rrhong_mr: gtscript.Field["float"],
    rsrimcg_mr: gtscript.Field["float"],
):

    with computation(PARALLEL), interval(...):
        rgsi = rvdepg + rsmltg + rraccsg + rsaccrg + rcrimsg + rsrimcg
        rgsi_mr = rrhong_mr + rsrimcg_mr
