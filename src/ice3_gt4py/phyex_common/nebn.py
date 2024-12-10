# -*- coding: utf-8 -*-
from dataclasses import dataclass, field
from enum import Enum
from typing import Literal

from ifs_physics_common.utils.f2py import ported_class


class FracIceAdjust(Enum):
    """Enumeration for ice fraction adjustments modes

    T in case of AROME

    """

    T = 0
    O = 1
    N = 2
    S = 3


class FracIceShallow(Enum):
    """Enumeration of ice fraction for shallow mass fluxes

    T in case of AROME
    """

    T = 0
    S = 1


class Condens(Enum):
    """Enumeration for condensation variance

    HCONDENS in .F90
    CB02 for AROME
    """

    CB02 = 0
    GAUS = 1


class Lambda3(Enum):
    """LAMBDA3 in AROME

    CB by default in AROME
    """

    CB = 0


@ported_class(from_file="PHYEX/src/common/aux/modd_nebn.F90")
@dataclass
class Neb:
    """Declaration of

    Args:
        tminmix (float): minimum temperature for mixed phase
        tmaxmix (float): maximum temperature for mixed phase
        hgt_qs (float): switch for height dependant VQSIGSAT
        frac_ice_adjust (str): ice fraction for adjustments
        frac_ice_shallow (str): ice fraction for shallow_mf
        vsigqsat (float): coeff applied to qsat variance contribution
        condens (str): subgrid condensation PDF
        lambda3 (str): lambda3 choice for subgrid cloud scheme
        statnw (bool): updated full statistical cloud scheme
        sigmas (bool): switch for using sigma_s from turbulence scheme
        subg_cond (bool): switch for subgrid condensation

    """

    HPROGRAM: Literal["AROME", "MESO-NH", "LMDZ"]

    TMINMIX: float = field(default=273.16)
    TMAXMIX: float = field(default=253.16)
    LHGT_QS: bool = field(default=False)
    FRAC_ICE_ADJUST: int = field(default=FracIceAdjust.S.value)
    FRAC_ICE_SHALLOW: int = field(default=FracIceShallow.S.value)
    VSIGQSAT: float = field(default=0.02)
    CONDENS: int = field(default=Condens.CB02.value)
    LAMBDA3: int = field(default=Lambda3.CB.value)
    LSTATNW: bool = field(default=False)
    LSIGMAS: bool = field(default=True)
    LSUBG_COND: bool = field(default=False)

    def __post_init__(self):
        if self.HPROGRAM == "AROME":
            self.FRAC_ICE_ADJUST = FracIceAdjust.T.value
            self.FRAC_ICE_SHALLOW = FracIceShallow.T.value
            self.VSIGQSAT = 0.02
            self.LSIGMAS = True
            self.LSUBG_COND = True

        elif self.HPROGRAM == "LMDZ":
            self.LSUBG_COND = True
