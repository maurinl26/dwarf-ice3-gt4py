# -*- coding: utf-8 -*-
from __future__ import annotations

from gt4py.cartesian import gtscript
from gt4py.cartesian.gtscript import interval, computation, PARALLEL
from ifs_physics_common.framework.stencil import stencil_collection

from ifs_physics_common.utils.f2py import ported_method


@ported_method(
    from_file="PHYEX/src/common/micro/mode_ice4_stepping.F90",
    from_line=225,
    to_line=228,
)
@stencil_collection("ice4_stepping_tsoft_init")
def ice4_stepping_init_tsoft(
    t_micro: gtscript.Field["float"], t_soft: gtscript.Field["float"]
):
    """Initialise t_soft with value of t_micro after each loop
    on LSOFT condition.

    Args:
        t_micro (gtscript.Field[float]): time for microphsyics loops
        t_soft (gtscript.Field[float]): time for lsoft blocks loops
    """

    with computation(PARALLEL), interval(...):
        t_soft = t_micro
