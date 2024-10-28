# -*- coding: utf-8 -*-
from abc import abstractmethod
from functools import cached_property, partial
import fmodpy
from ifs_physics_common.framework.config import GT4PyConfig
from ifs_physics_common.framework.components import ComputationalGridComponent
from ifs_physics_common.framework.grid import ComputationalGrid
import numpy as np

from ice3_gt4py.phyex_common.phyex import Phyex

from ifs_physics_common.framework.config import GT4PyConfig
from ifs_physics_common.framework.grid import ComputationalGrid

from pathlib import Path

import logging


class FieldInterface:
    
    @cached_property
    @abstractmethod
    def externals(self):
        """Dictionnary of externals
        key is fortran_name
        value is the value stored in phyex dataclasses
        """
        pass
        

    @cached_property
    @abstractmethod
    def dims(self) -> dict:
        """Compute fortran dimensions based on gt4py grid attributes

        Returns:
            dict: _description_
        """
        pass


    @cached_property
    @abstractmethod
    def array_shape(self) -> dict:
        pass

    @cached_property
    @abstractmethod
    def fields_in(self):
        pass

    @cached_property
    @abstractmethod
    def fields_out(self):
        pass

    @cached_property
    @abstractmethod
    def fields_inout(self):
        pass



class FortranComponent(ComputationalGridComponent, FieldInterface):
    def __init__(
        self,
        computational_grid: ComputationalGrid,
        gt4py_config: GT4PyConfig,
        externals: dict,
        fortran_subroutine: str,
        fortran_script: str,
        fortran_module: str,
        gt4py_stencil: str,
    ) -> None:
        super(self, ComputationalGridComponent).__init__(
            computational_grid=computational_grid, gt4py_config=gt4py_config
        )
        self.externals = externals
        self.compile_fortran_stencil(
            fortran_module=fortran_module,
            fortran_script=fortran_script,
            fortran_subroutine=fortran_subroutine,
        )


    def compile_fortran_stencil(
        self, fortran_script, fortran_module, fortran_subroutine
    ):
        current_directory = Path.cwd()
        logging.info(f"Root directory {current_directory}")
        root_directory = current_directory

        stencils_directory = Path(
            root_directory, "src", "ice3_gt4py", "stencils_fortran"
        )
        script_path = Path(stencils_directory, fortran_script)

        logging.info(f"Fortran script path {script_path}")
        self.fortran_script = fmodpy.fimport(script_path)
        fortran_module = getattr(self.fortran_script, fortran_module)
        self.fortran_stencil = getattr(fortran_module, fortran_subroutine)

    ##### Calls #####
    def call_fortran_stencil(self, fields: dict):
        """Call fortran stencil on a given field dict

        externals and dims are handled by component attributes itself

        Args:
            fields (dict): dictionnary of numpy arrays
        """
        
        ############# Preparing input ########
        field_attributes = {**self.fields_in, **self.fields_out, **self.fields_inout}
        state_fortran = dict()
        for key, array in fields.items():
            fortran_name = field_attributes[key]["fortran_name"]
            state_fortran.update({fortran_name: array})


        for field_name, array in state_fortran.items():
            logging.info(f"Field name {field_name}, array shape {array.shape}, array type {type(array)}")
            
        logging.info(f"Input, dims {self.dims}")
        logging.info(f"Input, externals {self.externals}")

        ############### Call ##################
        output_fields = self.fortran_stencil(
                **self.dims, **self.externals, **state_fortran 
            )
        
        ############## Preparing output #######
        output_fields_attributes = {**self.fields_inout, **self.fields_out}
        logging.info(f"Length of outputs {len(output_fields_attributes)}")
        if len(output_fields_attributes) == 1:
            logging.info(f"{output_fields_attributes}")
            return {next(iter(output_fields_attributes.keys())): np.array(output_fields)}
        elif len(output_fields_attributes) > 1:
            return {
                    value["fortran_name"]: output_fields[i] for i, value in enumerate(output_fields_attributes.values())
                }
            
class GT4PyComponent(ComputationalGridComponent, FieldInterface):
    def __init__(
        self,
        computational_grid: ComputationalGrid,
        gt4py_config: GT4PyConfig,
        externals: dict,
        fortran_subroutine: str,
        fortran_script: str,
        fortran_module: str,
        gt4py_stencil: str,
    ) -> None:
        super(self, ComputationalGridComponent).__init__(
            computational_grid=computational_grid, gt4py_config=gt4py_config
        )
        self.externals = externals
        self.compile_fortran_stencil(
            fortran_module=fortran_module,
            fortran_script=fortran_script,
            fortran_subroutine=fortran_subroutine,
        )

        self.compile_gt4py_stencil(gt4py_stencil, self.externals)


        



    #### Compilations ####
    def compile_gt4py_stencil(self, gt4py_stencil: str, externals: dict):
        """Compile GT4Py script given

        Args:
            gt4py_stencil (str): _description_
            externals (dict): _description_
        """
        self.gt4py_stencil = self.compile_stencil(gt4py_stencil, externals)


    def call_gt4py_stencil(self, fields: dict):
        """Call gt4py_stencil from a numpy array"""
        self.gt4py_stencil(**fields)
        return fields

