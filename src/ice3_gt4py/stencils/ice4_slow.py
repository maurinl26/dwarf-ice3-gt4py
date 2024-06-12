# -*- coding: utf-8 -*-
from __future__ import annotations

from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import exp, __INLINED, computation, interval, PARALLEL
from ifs_physics_common.framework.stencil import stencil_collection
from ifs_physics_common.utils.f2py import ported_method


@ported_method(from_file="PHYEX/src/common/micro/mode_ice4_slow.F90")
@stencil_collection("ice4_slow")
def ice4_slow(
    ldcompute: gtscript.Field["bool"],
    rhodref: gtscript.Field["float"],
    t: gtscript.Field["float"],
    ssi: gtscript.Field["float"],
    lv_fact: gtscript.Field["float"],
    ls_fact: gtscript.Field["float"],
    rv_t: gtscript.Field["float"],
    rc_t: gtscript.Field["float"],
    ri_t: gtscript.Field["float"],
    rs_t: gtscript.Field["float"],
    rg_t: gtscript.Field["float"],
    lbdas: gtscript.Field["float"],
    lbdag: gtscript.Field["float"],
    ai: gtscript.Field["float"],
    cj: gtscript.Field["float"],
    hli_hcf: gtscript.Field["float"],
    hli_hri: gtscript.Field["float"],
    rc_honi_tnd: gtscript.Field["float"],
    rv_deps_tnd: gtscript.Field["float"],
    ri_aggs_tnd: gtscript.Field["float"],
    ri_auts_tnd: gtscript.Field["float"],
    rv_depg_tnd: gtscript.Field["float"],
):
    """Compute the slow processes

    Args:
        ldcompute (gtscript.Field[float]): switch to activate processes computation on column
        rhodref (gtscript.Field[float]): reference density
        t (gtscript.Field[float]): temperature
        ssi (gtscript.Field[float]): supersaturation over ice
        lv_fact (gtscript.Field[float]): vaporisation latent heat over heat capacity
        ls_fact (gtscript.Field[float]): sublimation latent heat over heat capacity
        rv_t (gtscript.Field[float]): vapour mixing ratio at t
        ri_t (gtscript.Field[float]): ice m.r. at t
        rs_t (gtscript.Field[float]): snow m.r. at t
        rg_t (gtscript.Field[float]): graupel m.r. at t
        lbdag (gtscript.Field[float]): slope parameter of the graupel distribution
        lbdas (gtscript.Field[float]): slope parameter of the snow distribution
        ai (gtscript.Field[float]): thermodynamical function
        cj (gtscript.Field[float]): function to compute the ventilation factor
        hli_hcf (gtscript.Field[float]): low clouds cloud fraction
        hli_hri (gtscript.Field[float]): low clouds ice mixing ratio
        rc_honi_tnd (gtscript.Field[float]): homogeneous nucelation
        rv_deps_tnd (gtscript.Field[float]): deposition on snow
        ri_aggs_tnd (gtscript.Field[float]): aggregation on snow
        ri_auts_tnd (gtscript.Field[float]): autoconversion of ice
        rv_depg_tnd (gtscript.Field[float]): deposition on graupel
    """

    from __externals__ import (
        ACRIAUTI,
        ALPHA3,
        BCRIAUTI,
        BETA3,
        C_RTMIN,
        CEXVT,
        COLEXIS,
        CRIAUTI,
        EX0DEPG,
        EX0DEPS,
        EX1DEPG,
        EX1DEPS,
        EXIAGGS,
        FIAGGS,
        G_RTMIN,
        HON,
        I_RTMIN,
        O0DEPG,
        O0DEPS,
        O1DEPG,
        O1DEPS,
        S_RTMIN,
        TEXAUTI,
        TIMAUTI,
        TT,
        V_RTMIN,
        LDSOFT,
    )

    # 3.2 compute the homogeneous nucleation source : RCHONI
    with computation(PARALLEL), interval(...):
        if t < TT - 35.0 and rc_t > C_RTMIN and ldcompute:
            if __INLINED(not LDSOFT):
                rc_honi_tnd = min(
                    1000, HON * rhodref * rc_t * exp(ALPHA3 * (t - TT) - BETA3)
                )

        else:
            rc_honi_tnd = 0

    # 3.4 compute the deposition, aggregation and autoconversion sources
    # 3.4.3 compute the deposition on r_s : RVDEPS
    with computation(PARALLEL), interval(...):
        if rv_t < V_RTMIN and rs_t < S_RTMIN and ldcompute:
            # Translation note : #ifdef REPRO48 l118 to 120 kept
            # Translation note : #else REPRO48  l121 to 126 omitted

            if __INLINED(not LDSOFT):
                rv_deps_tnd = (ssi / (rhodref * ai)) * (
                    O0DEPS * lbdas**EX0DEPS + O1DEPS * cj * lbdas**EX1DEPS
                )

        else:
            rv_deps_tnd = 0

    # 3.4.4 compute the aggregation on r_s: RIAGGS
    with computation(PARALLEL), interval(...):
        if ri_t > I_RTMIN and rs_t > S_RTMIN and ldcompute:
            # Translation note : #ifdef REPRO48 l138 to 142 kept
            # Translation note : #else REPRO48 l143 to 150 omitted
            if __INLINED(not LDSOFT):
                ri_aggs_tnd = (
                    FIAGGS
                    * exp(COLEXIS * (t - TT))
                    * ri_t
                    * lbdas**EXIAGGS
                    * rhodref ** (-CEXVT)
                )

        # Translation note : OELEC = False l151 omitted
        else:
            ri_aggs_tnd = 0

    # 3.4.5 compute the autoconversion of r_i for r_s production: RIAUTS
    with computation(PARALLEL), interval(...):
        if hli_hri > I_RTMIN and ldcompute:
            if __INLINED(not LDSOFT):
                criauti_tmp = min(CRIAUTI, 10 ** (ACRIAUTI * (t - TT) + BCRIAUTI))
                ri_auts_tnd = (
                    TIMAUTI
                    * exp(TEXAUTI * (t - TT))
                    * max(0, hli_hri - criauti_tmp * hli_hcf)
                )

        else:
            ri_auts_tnd = 0

    # 3.4.6 compute the depsoition on r_g: RVDEPG
    with computation(PARALLEL), interval(...):
        if rv_t > V_RTMIN and rg_t > G_RTMIN and ldcompute:
            if __INLINED(not LDSOFT):
                rv_depg_tnd = (ssi / (rhodref * ai)) * (
                    O0DEPG * lbdag**EX0DEPG + O1DEPG * cj * lbdag**EX1DEPG
                )

        else:
            rv_depg_tnd = 0
