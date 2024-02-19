# -*- coding: utf-8 -*-
from __future__ import annotations

from gt4py.cartesian.gtscript import Field, IJ, function
from gt4py.cartesian.gtscript import exp, log, sqrt, floor, atan
from ice3_gt4py.functions.compute_ice_frac import compute_frac_ice
from ice3_gt4py.functions.src_1d import src_1d
from ice3_gt4py.functions.temperature import update_temperature


from ifs_physics_common.framework.stencil import stencil_collection

#
@stencil_collection("upstream_sedimentation")
def upstream_sedimentation(
    dt: "float",
    rhodref: Field["float"],
    oorhodz: Field["float"],  # 1 / (rho * dz)
    dz: Field["float"],
    pabst: Field["float"],
    tht: Field["float"],
    rc_in: Field["float"],
    rr_in: Field["float"],
    ri_in: Field["float"],
    rs_in: Field["float"],
    rg_in: Field["float"],
    ray_in: Field["float"],  # Optional PRAY
    lbc_in: Field["float"],  # Optional LBC
    fsedc_in: Field["float"],  # Optional PFSEDC
    conc3d_in: Field["float"],  # Optional PCONC3D
    rcs_tnd: Field["float"],
    rrs_tnd: Field["float"],
    ris_tnd: Field["float"],
    rss_tnd: Field["float"],
    rgs_tnd: Field["float"],
    dt__rho_dz_tmp: Field["float"],
    inst_rr_out: Field[IJ, "float"],
    inst_rc_out: Field[IJ, "float"],
    inst_ri_out: Field[IJ, "float"],
    inst_rs_out: Field[IJ, "float"],
    inst_rg_out: Field[IJ, "float"],
    wsed_tmp: Field["float"],  # sedimentation fluxes
    t_tmp: Field["float"],  # temperature - temporary field
    remaining_time: Field[IJ, "float"],  # remaining time until time step end
    max_tstep: Field[IJ, "float"],
    fpr_c_out: Field["float"],  # upper-air precip fluxes for cloud droplets
    fpr_r_out: Field["float"],
    fpr_s_out: Field["float"],
    fpr_i_out: Field["float"],
    fpr_g_out: Field["float"],
):
    from __externals__ import (
        lbexc,
        p00,
        t00,  # 293.15
        Rd,
        cpd,
        xcc,  # from iced
        cexvt,
        split_maxcfl,  # from parami - param_ice
        rholw,  # density of liquid water
        c_rtmin,  # min content for cloud droplets
    )

    # TODO
    # remaining time to be initialized

    # 2. Compute the fluxes
    # l219 to l262
    rcs_tnd -= rc_in / dt
    ris_tnd -= ri_in / dt
    rrs_tnd -= rr_in / dt
    rss_tnd -= rs_in / dt
    rgs_tnd -= rg_in / dt

    # in internal_sedim_split
    # For cloud droplets
    # TODO : encapsulation in do while
    # TODO: extend by functions to other species
    with computation(PARALLEL), interval(...):

        wlbdc = (lbc_in * conc3d_in / (rhodref * rc_in)) ** lbexc
        ray_tmp = ray_in / wlbdc

        # TODO : replace with exner
        t_tmp = tht * (pabst / p00) ** (Rd / cpd)
        wlbda = 6.6e-8 * (p00 / pabst) * (t_tmp / t00)
        cc_tmp = xcc * (1 + 1.26 * wlbda / ray_tmp)
        wsed_tmp = rhodref ** (-cexvt + 1) * wlbdc * cc_tmp * fsedc_in

    # l723 : main part
    with computation(PARALLEL), interval(...):

        max_tstep = maximum_time_step(
            c_rtmin, rhodref, max_tstep, rc_in, dz, wsed_tmp, remaining_time
        )

    # l733
    # Ground level
    with computation(PARALLEL), interval(0, 1):
        remaining_time[0, 0] -= max_tstep[0, 0]
        inst_rc_out[0, 0] += instant_precipitation(wsed_tmp, max_tstep, dt)

    # l738
    with computation(PARALLEL), interval(...):
        rcs_tnd = mixing_ratio_update(max_tstep, oorhodz, wsed_tmp, rcs_tnd, rc_in, dt)
        fpr_c_out += upper_air_flux(wsed_tmp, max_tstep, dt)


@function
def upper_air_flux(
    wsed_tmp: Field["float"],
    max_tstep: Field[IJ, "float"],
    dt: "float",
):

    return wsed_tmp * (max_tstep / dt)


@function
def mixing_ratio_update(
    max_tstep: Field[IJ, "float"],
    oorhodz: Field["float"],
    wsed: Field["float"],
    rs_tnd: Field["float"],
    r_in: Field["float"],
    dt: "float",
) -> Field["float"]:
    """Update mixing ratio

    Args:
        max_tstep (Field[IJ, float]): maximum time step to use
        oorhodz (Field[float]): 1 / (rho * dz)
        wsed (Field[float]): sedimentation flux
        rs_tnd (Field[float]): tendency for mixing ratio
        r_in (Field[float]): mixing ratio at time t
        dt (float): time step

    Returns:
        Field[float]: mixing ratio up to date
    """

    mrchange_tmp = max_tstep[0, 0] * oorhodz * (wsed[0, 0, 1] - wsed[0, 0, 0])
    r_in += mrchange_tmp + rs_tnd * max_tstep
    rs_tnd += mrchange_tmp / dt

    return rs_tnd


@function
def maximum_time_step(
    rtmin: "float",
    rhodref: Field["float"],
    max_tstep: Field[IJ, "float"],
    r: Field["float"],
    dz: Field["float"],
    wsed_tmp: Field["float"],
    remaining_time: Field["float"],
):
    from __externals__ import split_maxcfl

    tstep = max_tstep
    if r > rtmin and wsed_tmp > 1e-20 and remaining_time > 0:
        tstep[0, 0] = min(
            max_tstep,
            split_maxcfl
            * rhodref[0, 0, 0]
            * r[0, 0, 0]
            * dz[0, 0, 0]
            / wsed_tmp[0, 0, 0],
        )

    return tstep


@function
def instant_precipitation(
    wsed_tmp: Field["float"], max_tstep: Field["float"], dt: "float"
) -> Field["float"]:

    from __externals__ import rholw

    return wsed_tmp[0, 0, 0] / rholw * (max_tstep / dt)

