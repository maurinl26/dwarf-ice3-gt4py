# -*- coding: utf-8 -*-
from __future__ import annotations

from gt4py.cartesian.gtscript import (
    __INLINED,
    PARALLEL,
    atan,
    computation,
    exp,
    floor,
    interval,
    sqrt,
)
from gt4py.cartesian import gtscript
from ifs_physics_common.framework.stencil import stencil_collection
import numpy as np

from ice3_gt4py.functions.compute_ice_frac import compute_frac_ice
from ice3_gt4py.functions.ice_adjust import (
    sublimation_latent_heat,
    vaporisation_latent_heat,
)
from ice3_gt4py.functions.temperature import update_temperature
from ice3_gt4py.functions.tiwmx import e_sat_i, e_sat_w


@stencil_collection("ice_adjust")
def ice_adjust(
    sigqsat: gtscript.Field["float"],
    exnref: gtscript.Field["float"],
    exn: gtscript.Field["float"],
    rhodref: gtscript.Field["float"],
    pabs: gtscript.Field["float"],
    sigs: gtscript.Field["float"],
    cf_mf: gtscript.Field["float"],
    rc_mf: gtscript.Field["float"],
    ri_mf: gtscript.Field["float"],
    th: gtscript.Field["float"],
    rv: gtscript.Field["float"],
    rc: gtscript.Field["float"],
    ri: gtscript.Field["float"],
    rr: gtscript.Field["float"],
    rs: gtscript.Field["float"],
    rg: gtscript.Field["float"],
    ths: gtscript.Field["float"],
    rvs: gtscript.Field["float"],
    rcs: gtscript.Field["float"],
    ris: gtscript.Field["float"],
    cldfr: gtscript.Field["float"],
    ifr: gtscript.Field["float"],
    hlc_hrc: gtscript.Field["float"],
    hlc_hcf: gtscript.Field["float"],
    hli_hri: gtscript.Field["float"],
    hli_hcf: gtscript.Field["float"],
    sigrc: gtscript.Field["float"],
    cph: gtscript.Field["float"],
    lv: gtscript.Field["float"],
    ls: gtscript.Field["float"],
    criaut: gtscript.Field["float"],
    inq1: gtscript.Field[np.int64],
    src_1d: gtscript.GlobalTable[("float", (34))],
    dt: "float",
):
    """Microphysical adjustments for specific contents due to condensation.

    Args:
        sigqsat (gtscript.Field[float]): external qsat variance contribution
        exnref (gtscript.Field[float]): reference exner pressure
        exn (gtscript.Field[float]): true exner pressure
        rhodref (gtscript.Field[float]): reference density
        pabs (gtscript.Field[float]): absolute pressure at time t
        sigs (gtscript.Field[float]): standard dev for sub-grid saturation       (from turbulence scheme)
        cf_mf (gtscript.Field[float]): convective mass flux cloud fraction       (from shallow convection)
        rc_mf (gtscript.Field[float]): convective mass flux liquid mixing ratio  (from shallow convection)
        ri_mf (gtscript.Field[float]): convective mass flux ice mixing ratio     (from shallow convection)
        th (gtscript.Field[float]): potential temperature
        rv (gtscript.Field[float]): water vapour m.r. to adjust
        rc (gtscript.Field[float]): cloud water m.r. to adjust
        ri (gtscript.Field[float]): cloud ice m.r. to adjust
        rr (gtscript.Field[float]): rain water m.r. to adjust
        rs (gtscript.Field[float]): snow m.r. to adjust
        rg (gtscript.Field[float]): graupel m.r. to adjust
        ths (gtscript.Field[float]): potential temperature source
        rvs (gtscript.Field[float]): water vapour source
        rcs (gtscript.Field[float]): cloud droplets source
        ris (gtscript.Field[float]): ice source
        cldfr (gtscript.Field[float]): cloud fraction
        ifr (gtscript.Field[float]): ratio cloud ice moist part to dry part
        hlc_hrc (gtscript.Field[float]): high liquid content droplet m.r.
        hlc_hcf (gtscript.Field[float]): high liquid content cloud fraction
        hli_hri (gtscript.Field[float]): high liquid content ice m.r.
        hli_hcf (gtscript.Field[float]): high liquid content cloud fraction
        sigrc (gtscript.Field[float]): _description_
        dt (float): time step
    """

    from __externals__ import (
        ACRIAUTI,
        BCRIAUTI,
        CI,
        CL,
        CPD,
        CPV,
        CRIAUTC,
        CRIAUTI,
        NRR,
        RD,
        RV,
        SUBG_COND,
        SUBG_MF_PDF,
        TT,
        LAMBDA3,
    )

    # 2.3 Compute the variation of mixing ratio
    with computation(PARALLEL), interval(...):
        t = th * exn
        lv = vaporisation_latent_heat(t)
        ls = sublimation_latent_heat(t)

        # Rem
        rv_tmp = rv
        ri_tmp = ri
        rc_tmp = rc

    # Translation note : in Fortran, ITERMAX = 1, DO JITER =1,ITERMAX
    # Translation note : version without iteration is kept (1 iteration)
    #                   IF jiter = 1; CALL ITERATION()
    # jiter > 0

    # numer of moist variables fixed to 6 (without hail)

    # Translation note :
    # 2.4 specific heat for moist air at t+1
    with computation(PARALLEL), interval(...):
        # Translation note : case(7) removed because hail is not taken into account
        # Translation note : l453 to l456 removed
        if __INLINED(NRR == 6):
            cph = CPD + CPV * rv_tmp + CL * (rc_tmp + rr) + CI * (ri_tmp + rs + rg)
        if __INLINED(NRR == 5):
            cph = CPD + CPV * rv_tmp + CL * (rc_tmp + rr) + CI * (ri_tmp + rs)
        if __INLINED(NRR == 4):
            cph = CPD + CPV * rv_tmp + CL * (rc_tmp + rr)
        if __INLINED(NRR == 2):
            cph = CPD + CPV * rv_tmp + CL * rc_tmp + CI * ri_tmp

    # 3. subgrid condensation scheme
    # Translation note : only the case with LSUBG_COND = True retained (l475 in ice_adjust.F90)
    # sigqsat and sigs must be provided by the user
    with computation(PARALLEL), interval(...):
        cldfr = 0
        sigrc = 0

        # local gtscript.Fields
        # Translation note : 506 -> 514 kept (ocnd2 == False) # Arome default setting
        # Translation note : 515 -> 575 skipped (ocnd2 == True)
        prifact = 1  # ocnd2 == False for AROME
        ifr = 10
        frac_tmp = 0  # l340 in Condensation .f90

        # Translation note : 252 -> 263 if(present(PLV)) skipped (ls/lv are assumed to be present)
        # Translation note : 264 -> 274 if(present(PCPH)) skipped (files are assumed to be present)

        # store total water mixing ratio (244 -> 248)
        rt = rv_tmp + rc_tmp + ri_tmp * prifact

        # Translation note : 276 -> 310 (not osigmas) skipped (osigmas = True) for Arome default version
        # Translation note : 316 -> 331 (ocnd2 == True) skipped

        # l334 to l337
        pv = min(
            e_sat_w(t),
            0.99 * pabs,
        )
        piv = min(
            e_sat_i(t),
            0.99 * pabs,
        )

        # Translation note : OUSERI = True, OCND2 = False
        frac_tmp = ri_tmp / (rc_tmp + ri_tmp) if rc_tmp + ri_tmp > 1e-20 else 0
        frac_tmp = compute_frac_ice(t)

        # Supersaturation coefficients
        qsl = (RD / RV) * pv / (pabs - pv)
        qsi = (RD / RV) * piv / (pabs - piv)

        # # dtype_interpolate bewteen liquid and solid as a function of temperature
        qsl = (1 - frac_tmp) * qsl + frac_tmp * qsi
        lvs = (1 - frac_tmp) * lv + frac_tmp * ls

        # # coefficients a et b
        ah = lvs * qsl / (RV * t**2) * (1 + RV * qsl / RD)
        a = 1 / (1 + lvs / cph * ah)
        # # b = ah * a
        sbar = a * (rt - qsl + ah * lvs * (rc_tmp + ri_tmp * prifact) / CPD)

        # l369 - l390
        # Translation note : LSTATNW = False
        # Translation note : l381 retained for sigmas formulation
        sigma = (
            sqrt((2 * sigs) ** 2 + (sigqsat * qsl * a) ** 2) if sigqsat != 0 else sigs
        )

        # Translation note : l407 - l411
        sigma = max(1e-10, sigma)
        q1 = sbar / sigma

        # Translation notes : l413 to l468 skipped (HCONDENS=="GAUS")
        # Translation notes : l469 to l504 kept (HCONDENS = "CB02")
        # 9.2.3 Fractional cloudiness and cloud condensate

        # Translation note : l470 to l479
        if q1 > 0:
            cond_tmp = (
                min(exp(-1) + 0.66 * q1 + 0.086 * q1**2, 2) if q1 <= 2 else q1
            )  # we use the MIN function for continuity
        else:
            cond_tmp = exp(1.2 * q1 - 1)
        cond_tmp *= sigma

        # Translation note : l482 to l489
        # cloud fraction
        cldfr = max(0, min(1, 0.5 + 0.36 * atan(1.55 * q1))) if cond_tmp > 1e-12 else 0

        # Translation note : l487 to l489
        cond_tmp = 0 if cldfr == 0 else cond_tmp

        inq1 = floor(
            min(10, max(-22, min(-100, 2 * floor(q1))))
        )  # inner min/max prevents sigfpe when 2*zq1 does not fit dtype_into an "int"
        inc = 2 * q1 - inq1

        sigrc = min(1, (1 - inc) * src_1d.A[inq1] + inc * src_1d.A[inq1 + 1])

        # Translation note : l496 to l503
        hlc_hcf = 0
        hlc_hrc = 0
        hli_hcf = 0
        hli_hri = 0

        # # Translation notes : 506 -> 514 (not ocnd2)
        rc_tmp = (1 - frac_tmp) * cond_tmp  # liquid condensate
        ri_tmp = frac_tmp * cond_tmp  # solid condensate
        t = update_temperature(t, rc_tmp, rc_tmp, rc, ri, lv, ls)
        rv_tmp = rt - rc_tmp - ri_tmp * prifact

        # Transaltion notes : 566 -> 578 HLAMBDA3 = CB
        if __INLINED(LAMBDA3 == 0):
            sigrc *= min(3, max(1, 1 - q1))

    # Translation note : end jiter

    # l274 in ice_adjust.F90
    ##### 5.     COMPUTE THE SOURCES AND STORES THE CLOUD FRACTION #####
    with computation(PARALLEL), interval(...):
        # 5.0 compute the variation of mixing ratio
        w1 = (rc_tmp - rc) / dt
        w2 = (ri_tmp - ri) / dt

        # 5.1 compute the sources
        w1 = max(w1, -rcs) if w1 < 0 else min(w1, rvs)
        rvs -= w1
        rcs += w1
        ths += w1 * lv / (cph * exnref)

        w2 = max(w2, -ris) if w2 < 0 else min(w2, rvs)
        rvs -= w2
        rcs += w2
        ths += w2 * ls / (cph * exnref)

        # 5.2  compute the cloud fraction cldfr
        if __INLINED(not SUBG_COND):
            cldfr = 1 if (rcs + ris > 1e-12 / dt) else 0
        # Translation note : OCOMPUTE_SRC is taken False
        # Translation note : l320 to l322 removed

        # Translation note : LSUBG_COND = TRUE for Arome
        else:
            w1 = rc_mf / dt
            w2 = ri_mf / dt

            if w1 + w2 > rvs:
                w1 *= rvs / (w1 + w2)
                w2 = rvs - w1

            cldfr = min(1, cldfr + cf_mf)
            rcs += w1
            ris += w2
            rvs -= w1 + w2
            ths += (w1 * lv + w2 * ls) / cph / exnref

            # Droplets subgrid autoconversion
            # with computation(PARALLEL), interval(...):
            # LLHLC_H is True (AROME like)
            #
            criaut = CRIAUTC / rhodref

            # ice_adjust.F90 IF LLNONE; IF CSUBG_MF_PDF is None
            if __INLINED(SUBG_MF_PDF == 0):
                if w1 * dt > cf_mf * criaut:
                    hlc_hrc += w1 * dt
                    hlc_hcf = min(1, hlc_hcf + cf_mf)

            # Translation note : if LLTRIANGLE in .F90
            elif __INLINED(SUBG_MF_PDF == 1):
                if w1 * dt > cf_mf * criaut:
                    hcf = 1 - 0.5 * (criaut * cf_mf / max(1e-20, w1 * dt)) ** 2
                    hr = w1 * dt - (criaut * cf_mf) ** 3 / (
                        3 * max(1e-20, w1 * dt) ** 2
                    )

                elif 2 * w1 * dt <= cf_mf * criaut:
                    hcf = 0
                    hr = 0

                else:
                    hcf = (2 * w1 * dt - criaut * cf_mf) ** 2 / (
                        2.0 * max(1.0e-20, w1 * dt) ** 2
                    )
                    hr = (
                        4.0 * (w1 * dt) ** 3
                        - 3.0 * w1 * dt * (criaut * cf_mf) ** 2
                        + (criaut * cf_mf) ** 3
                    ) / (3 * max(1.0e-20, w1 * dt) ** 2)

                hcf *= cf_mf
                hlc_hcf = min(1, hlc_hcf + hcf)
                hlc_hrc += hr

            # Ice subgrid autoconversion
            criaut = min(
                CRIAUTI,
                10 ** (ACRIAUTI * (t - TT) + BCRIAUTI),
            )

            if __INLINED(SUBG_MF_PDF == 0):
                if w2 * dt > cf_mf * criaut:
                    hli_hri += w2 * dt
                    hli_hcf = min(1, hli_hcf + cf_mf)

            elif __INLINED(SUBG_MF_PDF == 1):
                if w2 * dt > cf_mf * criaut:
                    hli_hcf = 1 - 0.5 * ((criaut * cf_mf) / (w2 * dt)) ** 2
                    hli_hri = w2 * dt - (criaut * cf_mf) ** 3 / (3 * (w2 * dt) ** 2)

                elif 2 * w2 * dt <= cf_mf * criaut:
                    hli_hcf = 0
                    hli_hri = 0

                else:
                    hli_hcf = (2 * w2 * dt - criaut * cf_mf) ** 2 / (
                        2.0 * (w2 * dt) ** 2
                    )
                    hli_hri = (
                        4.0 * (w2 * dt) ** 3
                        - 3.0 * w2 * dt * (criaut * cf_mf) ** 2
                        + (criaut * cf_mf) ** 3
                    ) / (3 * (w2 * dt) ** 2)

                hli_hcf *= cf_mf
                hli_hcf = min(1, hli_hcf + hli_hcf)
                hli_hri += hli_hri

    # Translation note : 402 -> 427 (removed pout_x not present )
