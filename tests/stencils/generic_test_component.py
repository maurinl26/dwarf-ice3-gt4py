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


class TestComponent(ComputationalGridComponent):
    def __init__(
        self,
        computational_grid: ComputationalGrid,
        gt4py_config: GT4PyConfig,
        phyex: Phyex,
        fortran_subroutine: str,
        fortran_script: str,
        fortran_module: str,
        gt4py_stencil: str,
    ) -> None:
        super().__init__(
            computational_grid=computational_grid, gt4py_config=gt4py_config
        )
        self.phyex_externals = phyex.to_externals()

        self.compile_fortran_stencil(
            fortran_module=fortran_module,
            fortran_script=fortran_script,
            fortran_subroutine=fortran_subroutine,
        )

        self.compile_gt4py_stencil(gt4py_stencil, self.phyex_externals)

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

    #### Compilations ####
    def compile_gt4py_stencil(self, gt4py_stencil: str, externals: dict):
        """Compile GT4Py script given

        Args:
            gt4py_stencil (str): _description_
            externals (dict): _description_
        """
        self.gt4py_stencil = self.compile_stencil(gt4py_stencil, externals)

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
        field_attributes = {**self.fields_in, **self.fields_out, **self.fields_inout}
        state_fortran = dict()
        for key, field in fields.items():
            fortran_name = field_attributes[key]["fortran_name"]
            state_fortran.update({fortran_name: field})

        output_fields_tuple = self.fortran_stencil(
            **state_fortran, **self.dims, **self.externals
        )
        
        output_fields = dict()
        fields_to_name = {**self.fields_inout, **self.fields_out}
        assert(len(output_fields_tuple) == len(list(fields_to_name.keys())))
        
        i = 0
        for field_name, field_attributes in fields_to_name.items():
            fortran_name = field_attributes["fortran_name"]
            output_fields.update({
                fortran_name: output_fields_tuple[i]
            })
            i+=1
            
        return output_fields

    def call_gt4py_stencil(self, fields: dict):
        """Call gt4py_stencil from a numpy array"""
        self.gt4py_stencil(**fields)
        return fields