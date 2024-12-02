# -*- coding: utf-8 -*-
import datetime
import logging
from functools import cached_property, partial
from pathlib import Path
from typing import TYPE_CHECKING, Literal, Tuple

import fmodpy
import numpy as np
from ifs_physics_common.framework.components import ComputationalGridComponent
from ifs_physics_common.framework.config import GT4PyConfig
from ifs_physics_common.framework.grid import ComputationalGrid, DimSymbol, I, J, K
from ifs_physics_common.framework.storage import allocate_data_array
from ifs_physics_common.utils.typingx import DataArray, NDArrayLikeDict

import ice3_gt4py.stencils
from ice3_gt4py.initialisation.utils import initialize_field
from ice3_gt4py.phyex_common.phyex import Phyex
from ice3_gt4py.utils.allocate import allocate


##### For tests ####
def allocate_state_ice4_rrhong(
    computational_grid: ComputationalGrid, gt4py_config: GT4PyConfig
) -> NDArrayLikeDict:
    """Allocate field to state keys following type (float, int, bool) and dimensions (2D, 3D).

    Args:
        computational_grid (ComputationalGrid): grid indexes
        gt4py_config (GT4PyConfig): gt4py configuration

    Returns:
        NDArrayLikeDict: dictionnary of field with associated keys for field name
    """

    def _allocate(
        grid_id: Tuple[DimSymbol, ...],
        units: str,
        dtype: Literal["bool", "float", "int"],
    ) -> DataArray:
        return allocate_data_array(
            computational_grid, grid_id, units, gt4py_config=gt4py_config, dtype=dtype
        )

    allocate_b_ij = partial(_allocate, grid_id=(I, J), units="", dtype="bool")
    allocate_b = partial(_allocate, grid_id=(I, J, K), units="", dtype="bool")
    allocate_f = partial(_allocate, grid_id=(I, J, K), units="", dtype="float")
    allocate_h = partial(_allocate, grid_id=(I, J, K - 1 / 2), units="", dtype="float")
    allocate_ij = partial(_allocate, grid_id=(I, J), units="", dtype="float")
    allocate_i_ij = partial(_allocate, grid_id=(I, J), units="", dtype="int")

    return {
        "ldcompute": allocate_b(),
        "exn": allocate_f(),
        "ls_fact": allocate_f(),
        "lv_fact": allocate_f(),
        "t": allocate_f(),
        "tht": allocate_f(),
        "rr_t": allocate_f(),
        "rrhong_mr": allocate_f(),
    }


class TestIce4RRHONG(ComputationalGridComponent):
    def __init__(
        self,
        computational_grid: ComputationalGrid,
        gt4py_config: GT4PyConfig,
        phyex: Phyex,
        fortran_module: str,
        fortran_subroutine: str,
        fortran_script: str,
        gt4py_stencil: str,
    ) -> None:
        super().__init__(
            computational_grid=computational_grid, gt4py_config=gt4py_config
        )

        externals = phyex.to_externals()
        self.ice4_rrhong_gt4py = self.compile_stencil(gt4py_stencil, externals)

        self.externals = {
            "xtt": externals["TT"],
            "r_rtmin": externals["R_RTMIN"],
            "lfeedbackt": externals["LFEEDBACKT"],
        }

        self.generate_gt4py_state()

        # aro_filter stands for the parts before 'call ice_adjust' in aro_adjust.f90
        # stencils_fortran_directory = "./src/ice3_gt4py/stencils_fortran"
        # fortran_script = "mode_ice4_rrhong.F90"

        project_dir = Path.cwd()
        stencils_fortran_dir = Path(
            project_dir, "src", "ice3_gt4py", "stencils_fortran"
        )

        fortran_script_path = Path(stencils_fortran_dir, fortran_script)
        self.fortran_script = fmodpy.fimport(fortran_script_path)
        self.fortran_module = fortran_module
        self.fortran_subroutine = fortran_subroutine
        fortran_module = getattr(self.fortran_script, self.fortran_module)
        self.fortran_stencil = getattr(fortran_module, self.fortran_subroutine)

    @cached_property
    def dims(self):
        return {"kproma": KPROMA, "ksize": KSIZE}

    @cached_property
    def fields_mapping(self):
        return {
            "ldcompute": "ldcompute",
            "pexn": "exn",
            "plsfact": "ls_fact",
            "plvfact": "lv_fact",
            "pt": "t",
            "ptht": "tht",
            "prrt": "rr_t",
            "prrhong_mr": "rrhong_mr",
        }

    @cached_property
    def fields_in(self):
        return {
            "ldcompute": np.ones((self.dims["kproma"]), dtype=np.int32),
            "pexn": np.array(np.random.rand(self.dims["kproma"]), "f", order="F"),
            "plvfact": np.array(np.random.rand(self.dims["kproma"]), "f", order="F"),
            "plsfact": np.array(np.random.rand(self.dims["kproma"]), "f", order="F"),
            "pt": np.array(np.random.rand(self.dims["kproma"]), "f", order="F"),
            "prrt": np.array(np.random.rand(self.dims["kproma"]), "f", order="F"),
            "ptht": np.array(np.random.rand(self.dims["kproma"]), "f", order="F"),
        }

    @cached_property
    def fields_out(self):
        return {"prrhong_mr": np.zeros((self.dims["kproma"]), "f", order="F")}

    def generate_gt4py_state(self):

        self.state_gt4py = allocate_state_ice4_rrhong(
            self.computational_grid, self.gt4py_config
        )
        fields = {**self.fields_in, **self.fields_out}
        for key_fortran, key_gt4py in self.fields_mapping.items():
            initialize_field(
                self.state_gt4py[key_gt4py], 2*fields[key_fortran][:, np.newaxis]
            )

    def test(self):
        """Call fortran stencil"""

        prrhong_mr_mean = self.fields_out['prrhong_mr'].mean()
        rrhong_mr = self.state_gt4py['rrhong_mr'].values.mean()
        logging.info(f"Input field, rrhong_mr (fortran) : {prrhong_mr_mean}")
        logging.info(f"Input field, rrhong_mr (gt4py) : {rrhong_mr}")

        plsfact_mean = self.fields_in["plsfact"].mean()
        ls_fact_mean = self.state_gt4py['ls_fact'].values.mean()
        logging.info(f"Input field, ls_fact (fortran) : {plsfact_mean}")
        logging.info(f"Input field, ls_fact (gt4py) : {ls_fact_mean}\n ")
        
        state_fortran = {
            **self.fields_in,
            **self.fields_out,
        }
        
        for field_name, array in state_fortran.items():
            logging.info(f"Fortran field name {field_name}, array shape {array.shape}, array type {type(array)}")
            print(array)

        field_fortran = self.fortran_stencil(
            **self.dims, **self.externals, **self.fields_in, **self.fields_out
        )
        

        self.ice4_rrhong_gt4py(
            **self.state_gt4py,
        )

        # logging.info(f"Mean Fortran : {prrhong_mr.mean()}")
        field_gt4py = self.state_gt4py["rrhong_mr"][...].values
        
        logging.info(f"Outputs")
        logging.info(f"Field fortran shape : {field_fortran.shape}")
        logging.info(f"GT4Py shape : {field_gt4py.shape}")
        logging.info(f"Field fortran mean : {field_fortran.mean()}")
        logging.info(f"GT4Py mean : {field_gt4py.mean()}\n")
        
        logging.info(f"Reshaping GT4Py")
        reshaped_gt4py_field = field_gt4py[:,:,:1].reshape(-1, 1).reshape(-1, 1)
        
        logging.info(f"GT4Py (Reshaped) shape : {reshaped_gt4py_field.shape}")
        logging.info(f"Fortran shape : {field_fortran.shape}")
        logging.info(f"Mean GT4Py (Reshaped) {reshaped_gt4py_field.mean()}")
        logging.info(f"Mean Fortran {field_fortran.mean()}")
        
        
        mean_absolute_error = abs(field_fortran - reshaped_gt4py_field).mean()
        logging.info(f"MAE : {mean_absolute_error}")
        


if __name__ == "__main__":

    KPROMA, KSIZE = 50, 50

    # TODO : set in env values
    backend = "gt:cpu_ifirst"
    rebuild = True
    validate_args = True

    phyex = Phyex(program="AROME")

    logging.info("Initializing grid ...")

    # Grid has only 1 dimension since fields are packed in fortran version
    grid = ComputationalGrid(50, 1, 1)
    dt = datetime.timedelta(seconds=1)

    ######## Backend and gt4py config #######
    logging.info(f"With backend {backend}")
    gt4py_config = GT4PyConfig(
        backend=backend, rebuild=rebuild, validate_args=validate_args, verbose=False
    )

    logging.info("Calling ice4_rrhong with dicts")

    test_ice4_rrhong = TestIce4RRHONG(
        computational_grid=grid,
        gt4py_config=gt4py_config,
        phyex=phyex,
        fortran_module="mode_ice4_rrhong",
        fortran_subroutine="ice4_rrhong",
        fortran_script="mode_ice4_rrhong.F90",
        gt4py_stencil="ice4_rrhong",
    ).test()