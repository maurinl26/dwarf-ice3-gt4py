# -*- coding: utf-8 -*-
from __future__ import annotations

from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import (
    __externals__,
    __INLINED,
    exp,
    log,
    computation,
    interval,
    PARALLEL,
)
from ifs_physics_common.framework.stencil import stencil_collection
from ifs_physics_common.utils.f2py import ported_method

from ice3_gt4py.functions.interp_micro import (
    index_interp_micro_1d,
    index_micro2d_acc_r,
    index_micro2d_acc_s,
)
from ice3_gt4py.functions.sign import sign


@ported_method(from_file="PHYEX/src/common/micro/mode_ice4_fast_rs.F90")
@stencil_collection("ice4_fast_rs")
def ice4_fast_rs(
    ldcompute: gtscript.Field["bool"],
    rhodref: gtscript.Field["float"],
    lv_fact: gtscript.Field["float"],
    ls_fact: gtscript.Field["float"],
    pres: gtscript.Field["float"],
    dv: gtscript.Field["float"],
    ka: gtscript.Field["float"],
    cj: gtscript.Field["float"],
    lbdar: gtscript.Field["float"],
    lbdas: gtscript.Field["float"],
    t: gtscript.Field["float"],
    rv_t: gtscript.Field["float"],
    rc_t: gtscript.Field["float"],
    rr_t: gtscript.Field["float"],
    rs_t: gtscript.Field["float"],
    riaggs: gtscript.Field["float"],
    rcrimss: gtscript.Field["float"],
    rcrimsg: gtscript.Field["float"],
    rsrimcg: gtscript.Field["float"],
    rraccss: gtscript.Field["float"],
    rraccsg: gtscript.Field["float"],
    rsaccrg: gtscript.Field["float"],
    rs_mltg_tnd: gtscript.Field["float"],
    rc_mltsr_tnd: gtscript.Field["float"],
    rs_rcrims_tnd: gtscript.Field["float"],
    rs_rcrimss_tnd: gtscript.Field["float"],
    rs_rsrimcg_tnd: gtscript.Field["float"],
    rs_rraccs_tnd: gtscript.Field["float"],
    rs_rraccss_tnd: gtscript.Field["float"],
    rs_rsaccrg_tnd: gtscript.Field["float"],
    rs_freez1_tnd: gtscript.Field["float"],
    rs_freez2_tnd: gtscript.Field["float"],
    gaminc_rim1: gtscript.GlobalTable[float, (80)],
    gaminc_rim2: gtscript.GlobalTable[float, (80)],
    gaminc_rim4: gtscript.GlobalTable[float, (80)],
    ker_raccs: gtscript.GlobalTable[float, (40, 40)],
    ker_raccss: gtscript.GlobalTable[float, (40, 40)],
    ker_saccrg: gtscript.GlobalTable[float, (40, 40)],
    index_floor: gtscript.Field["int"],
    index_floor_r: gtscript.Field["int"],
    index_floor_s: gtscript.Field["int"],
):
    """Computes fast processes for snow

    Args:
        ldcompute (gtscript.Field[bool]): mask for sources computations
        rhodref (gtscript.Field[float]): dry density of air
        lv_fact (gtscript.Field[float]): latent heat of vapourisation over heat capacity
        ls_fact (gtscript.Field[float]): latent heat of sublimation over heat capacity
        pres (gtscript.Field[float]): absolute pressure
        rr_t (gtscript.Field[float]): rain m.r. at t
        rs_t (gtscript.Field[float]): snow m.r. at t
        riaggs (gtscript.Field[float]): ice aggregation to snow
        rsrimcg (gtscript.Field[float]): snow riming over graupel
        rraccss (gtscript.Field[float]): accretion of rain and aggregates
        rsaccrg (gtscript.Field[float]): accretion of graupel
        rs_mltg_tnd (gtscript.Field[float]): melting of snow
        rs_rsrimcg_tnd (gtscript.Field[float]): heavy riming of the aggregates
        rs_rraccs_tnd (gtscript.Field[float]): accretion of rain and aggregates
        rs_rraccss_tnd (gtscript.Field[float]): accretion of rain and aggregates
        rs_rsaccrg_tnd (gtscript.Field[float]): accretion of rain and aggregates
        rs_freez1_tnd (gtscript.Field[float]): snow freezing source (tendency)
        rs_freez2_tnd (gtscript.Field[float]): snow freezing source (tendency)
        gaminc_rim1 (gtscript.GlobalTable[float,): look up table for riming
        gaminc_rim2 (gtscript.GlobalTable[float,): look up table for riming
        gaminc_rim4 (gtscript.GlobalTable[float,): look up table for riming
        ker_raccs (gtscript.GlobalTable[float,): look-up table for rain accretion over snow
        ker_raccss (gtscript.GlobalTable[float,): look-up table for rain accretion over snow
        ker_saccrg (gtscript.GlobalTable[float,): look-up table for snow accretion over graupel
        index_floor (gtscript.Field[int]): integer index for riming look up tables
        index_floor_r (gtscript.Field[int]): integer index for accretion look up tables
        index_floor_s (gtscript.Field[int]): integer index for accretion look up tables
    """
    from __externals__ import (
        LDSOFT,
        ALPI,
        ALPW,
        BETAI,
        BETAW,
        BS,
        C_RTMIN,
        CEXVT,
        CI,
        CL,
        CPV,
        CRIMSG,
        CRIMSS,
        CXS,
        EPSILO,
        ESTT,
        EX0DEPS,
        EX1DEPS,
        EXCRIMSG,
        EXCRIMSS,
        EXSRIMCG,
        EXSRIMCG2,
        FRACCSS,
        FSACCRG,
        FSCVMG,
        GAMI,
        GAMW,
        LBRACCS1,
        LBRACCS2,
        LBRACCS3,
        LBSACCR1,
        LBSACCR2,
        LBSACCR3,
        LEVLIMIT,
        LMTT,
        LVTT,
        O0DEPS,
        O1DEPS,
        R_RTMIN,
        RV,
        S_RTMIN,
        SNOW_RIMING,
        SRIMCG,
        SRIMCG2,
        SRIMCG3,
        TT,
    )

    # 5.0 maximum freezing rate
    with computation(PARALLEL), interval(...):
        # Translation note l106 removed not LDSOFT
        if rs_t < S_RTMIN and ldcompute:
            rs_freez1_tnd = rv_t * pres / (EPSILO + rv_t)
            if LEVLIMIT:
                rs_freez1_tnd = min(
                    rs_freez1_tnd, exp(ALPI - BETAI / t - GAMI * log(t))
                )

            rs_freez1_tnd = ka * (TT - t) + dv * (LVTT + (CPV - CL) * (t - TT)) * (
                ESTT - rs_freez1_tnd
            ) / (RV * t)

            # Translation note l115 to l177 kept #ifdef REPRO48
            # Translation note l119 to l122 removed #else REPRO48

            rs_freez1_tnd *= (
                O0DEPS * lbdas**EX0DEPS + O1DEPS * cj * lbdas**EX1DEPS
            ) / (rhodref * (LMTT - CL * (TT - t)))
            rs_freez2_tnd = (rhodref * (LMTT + (CI - CL) * (TT - t))) / (
                rhodref * (LMTT - CL * (TT - t))
            )

            # Translation note l129 removed
            freez_rate_tmp = 0

        else:
            rs_freez1_tnd = 0
            rs_freez2_tnd = 0
            freez_rate_tmp = 0

    # 5.1 cloud droplet riming of the aggregates
    with computation(PARALLEL), interval(...):
        if rc_t > C_RTMIN and rs_t > S_RTMIN and ldcompute:
            zw_tmp = lbdas

            # Translation note : l144 kept
            #                    l146 removed

            grim_tmp = True

        else:
            grim_tmp = False
            rs_rcrims_tnd = 0
            rs_rcrimss_tnd = 0
            rs_rsrimcg_tnd = 0

    # Interpolation + Lookup Table
    with computation(PARALLEL), interval(...):
        if __INLINED(not LDSOFT):
            if grim_tmp:
                # Translation note : LDPACK is False l46 to l88 removed in interp_micro.func.h
                #                                    l90 to l123 kept
                index_floor, index_float = index_interp_micro_1d(zw_tmp)
                zw1_tmp = (
                    index_float * gaminc_rim1.A[index_floor + 1]
                    + (1 - index_float) * gaminc_rim1.A[index_floor]
                )
                zw2_tmp = (
                    index_float * gaminc_rim2.A[index_floor + 1]
                    + (1 - index_float) * gaminc_rim2.A[index_floor]
                )
                zw3_tmp = (
                    index_float * gaminc_rim4.A[index_floor + 1]
                    + (1 - index_float) * gaminc_rim4.A[index_floor]
                )

    # 5.1.4 riming of the small sized aggregates
    with computation(PARALLEL), interval(...):
        if grim_tmp:
            # Translation note : #ifdef REPRO48 l170 to l172 kept
            #                                   l174 to l178 removed
            rs_rcrimss_tnd = (
                CRIMSS * zw1_tmp * rc_t * lbdas**EXCRIMSS * rhodref ** (-CEXVT)
            )

    # 5.1.6 riming convesion of the large size aggregates
    with computation(PARALLEL), interval(...):
        if grim_tmp:
            # Translation note : #ifdef REPRO48 l189 to l191 kept
            #                                   l193 to l197 removed
            rs_rcrims_tnd = CRIMSG * rc_t * lbdas**EXCRIMSG * rhodref ** (-CEXVT)

    # if parami  csnowriming == M90
    with computation(PARALLEL), interval(...):
        # PARAMI%CSNOWRIMING == M90
        # TODO : refactor if statement out of stencil for performance
        if __INLINED(SNOW_RIMING == 0):
            if grim_tmp:
                zw_tmp = rs_rsrimcg_tnd - rs_rcrimss_tnd
                # Translation note : #ifdef REPRO48 l208 kept
                #                                   l210 and l211 removed
                rs_rsrimcg_tnd = SRIMCG * lbdas**EXSRIMCG * (1 - zw2_tmp)

                # Translation note : #ifdef REPRO48 l214 to l217 kept
                #                                   l219 to l223 removed
                rs_rsrimcg_tnd = (
                    zw_tmp
                    * rs_rsrimcg_tnd
                    / max(
                        1e-20,
                        SRIMCG3 * SRIMCG2 * lbdas**EXSRIMCG2 * (1 - zw3_tmp)
                        - SRIMCG3 * rs_rsrimcg_tnd,
                    )
                )

        else:
            rs_rsrimcg_tnd = 0

    #
    with computation(PARALLEL), interval(...):
        if grim_tmp and t < TT:
            rcrimss = min(freez_rate_tmp, rs_rcrimss_tnd)
            freez_rate_tmp = max(0, freez_rate_tmp - rcrimss)

            # proportion we are able to freeze
            zw0_tmp = min(1, freez_rate_tmp / max(1e-20, rs_rcrims_tnd - rcrimss))
            rcrimsg = zw0_tmp * max(0, rs_rcrims_tnd - rcrimss)  # rc_rimsg
            freez_rate_tmp = max(0, freez_rate_tmp - rcrimsg)
            rsrimcg = zw0_tmp * rs_rsrimcg_tnd

            rsrimcg *= max(0, -sign(1, -rcrimsg))
            rcrimsg = max(0, rcrimsg)

        else:
            rcrimss = 0
            rcrimsg = 0
            rsrimcg = 0

    # 5.2. rain accretion onto the aggregates
    with computation(PARALLEL), interval(...):
        if rr_t > R_RTMIN and rs_t > S_RTMIN and ldcompute:
            gacc_tmp = True
        else:
            gacc_tmp = False
            rs_rraccs_tnd = 0
            rs_rraccss_tnd = 0
            rs_rsaccrg_tnd = 0

    with computation(PARALLEL), interval(...):
        # Translation note : LDPACK is False l159 to l223 removed in interp_micro.func.h
        #                                    l226 to l266 kept

        if __INLINED(not LDSOFT):
            if gacc_tmp:
                rs_rraccs_tnd = 0
                rs_rraccss_tnd = 0
                rs_rsaccrg_tnd = 0

                index_floor_r, index_float_r = index_micro2d_acc_r(lbdar)
                index_floor_s, index_float_s = index_micro2d_acc_s(lbdas)
                zw1_tmp = index_float_s * (
                    index_float_r * ker_raccss.A[index_floor_s + 1, index_floor_r + 1]
                    + (1 - index_float_r)
                    * ker_raccss.A[index_floor_s + 1, index_floor_r]
                ) + (1 - index_float_s) * (
                    index_float_r * ker_raccss.A[index_floor_s, index_floor_r + 1]
                    + (1 - index_float_r) * ker_raccss.A[index_floor_s, index_floor_r]
                )

                zw2_tmp = index_float_s * (
                    index_float_r * ker_raccs.A[index_floor_s + 1, index_floor_r + 1]
                    + (1 - index_float_r)
                    * ker_raccs.A[index_floor_s + 1, index_floor_r]
                ) + (1 - index_float_s) * (
                    index_float_r * ker_raccs.A[index_floor_s, index_floor_r + 1]
                    + (1 - index_float_r) * ker_raccs.A[index_floor_s, index_floor_r]
                )

                zw3_tmp = index_float_s * (
                    index_float_r * ker_saccrg.A[index_floor_s + 1, index_floor_r + 1]
                    + (1 - index_float_r)
                    * ker_saccrg.A[index_floor_s + 1, index_floor_r]
                ) + (1 - index_float_s) * (
                    index_float_r * ker_saccrg.A[index_floor_s, index_floor_r + 1]
                    + (1 - index_float_r) * ker_saccrg.A[index_floor_s, index_floor_r]
                )

    #         # CALL INTERP_MICRO_2D

    # 5.2.4. raindrop accreation on the small sized aggregates
    with computation(PARALLEL), interval(...):
        if gacc_tmp:
            # Translation note : REPRO48 l279 to l283 kept
            #                            l285 to l289 removed

            zw_tmp = (
                FRACCSS
                * (lbdas**CXS)
                * (rhodref ** (-CEXVT))
                * (
                    LBRACCS1 / (lbdas**2)
                    + LBRACCS2 / (lbdas * lbdar)
                    + LBRACCS3 / (lbdar**2)
                )
                / lbdar**4
            )

    # 5.2.6 raindrop accretion-conversion of the large sized aggregates
    with computation(PARALLEL), interval(...):
        if gacc_tmp:
            rs_rsaccrg_tnd = (
                FSACCRG
                * zw3_tmp
                * (lbdas ** (CXS - BS))
                * (rhodref ** (-CEXVT - 1))
                * (
                    LBSACCR1 / (lbdas**2)
                    + LBSACCR2 / (lbdar * lbdas)
                    + LBSACCR3 / (lbdas**2)
                )
                / lbdar
            )

    # l324
    # More restrictive ACC mask to be used for accretion by negative temperature only
    with computation(PARALLEL), interval(...):
        if gacc_tmp and t < TT:
            rraccss = min(freez_rate_tmp, rs_rraccss_tnd)
            freez_rate_tmp = max(0, freez_rate_tmp - rraccss)

            # proportion we are able to freeze
            zw_tmp = min(1, freez_rate_tmp / max(1e-20, rs_rraccss_tnd - rraccss))
            rraccsg = zw_tmp * max(0, rs_rraccs_tnd - rraccss)
            freez_rate_tmp = max(0, freez_rate_tmp - rraccsg)
            rsaccrg = zw_tmp * rs_rsaccrg_tnd

            rsaccrg *= max(0, -sign(1, -rraccsg))
            rraccsg = max(0, rraccsg)

        else:
            rraccss = 0
            rraccsg = 0
            rsaccrg = 0

    # 5.3 Conversion-Melting of the aggregates
    with computation(PARALLEL), interval(...):
        if rs_t < S_RTMIN and t > TT and ldcompute:
            if __INLINED(not LDSOFT):
                rs_mltg_tnd = rv_t * pres / (EPSILO + rv_t)
                if LEVLIMIT:
                    rs_mltg_tnd = min(
                        rs_mltg_tnd, exp(ALPW - BETAW / t - GAMW * log(t))
                    )
                rs_mltg_tnd = ka * (TT - t) + (
                    dv
                    * (LVTT + (CPV - CL) * (t - TT))
                    * (ESTT - rs_mltg_tnd)
                    / (RV * t)
                )

                # Tranlsation note : #ifdef REPRO48 l360 to l365 kept
                #                                   l367 to l374 removed
                rs_mltg_tnd = FSCVMG * max(
                    0,
                    (
                        -rs_mltg_tnd
                        * (O0DEPS * lbdas**EX0DEPS + O1DEPS * cj * lbdas * EX1DEPS)
                        - (rs_rcrims_tnd + rs_rraccs_tnd) * (rhodref * CL * (TT - t))
                    )
                    / (rhodref * LMTT),
                )

                # note that RSCVMG = RSMLT*XFSCVMG but no heat is exchanged (at the rate RSMLT)
                # because the graupeln produced by this process are still icy###
                #
                # When T < XTT, rc is collected by snow (riming) to produce snow and graupel
                # When T > XTT, if riming was still enabled, rc would produce snow and graupel with snow becomming graupel (conversion/melting) and graupel becomming rain (melting)
                # To insure consistency when crossing T=XTT, rc collected with T>XTT must be transformed in rain.
                # rc cannot produce iced species with a positive temperature but is still collected with a good efficiency by snow

                rc_mltsr_tnd = rs_rcrims_tnd

        else:
            rs_mltg_tnd = 0
            rc_mltsr_tnd = 0
