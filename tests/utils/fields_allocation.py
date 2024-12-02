# -*- coding: utf-8 -*-
import logging

import numpy as np
from ifs_physics_common.framework.components import ComputationalGridComponent 
from ifs_physics_common.framework.grid import I, J, K
from ifs_physics_common.utils.typingx import NDArrayLikeDict
from utils.array_dict_operations import remove_y_axis, unpack
from utils.initialize_fields import initialize_field
from utils.allocate_state import allocate_state

from repro.default_config import test_grid, default_gt4py_config


####### Field allocation functions #######
def allocate_gt4py_fields(
    component: ComputationalGridComponent, 
    fields: dict
) -> NDArrayLikeDict:
    """Allocate storage for gt4py fields and
    initialize fields with given np arrays

    Args:
        component (ComputationalGridComponent): a ComputationalGridComponent with well described fields
        fields (dict): _description_

    Returns:
        NDArrayLikeDict: _description_
    """

    fields_metadata = {
        **component.fields_in,
        **component.fields_out,
        **component.fields_inout,
    }
    state_gt4py = allocate_state(test_grid, default_gt4py_config, fields_metadata)
    for key, field_array in fields.items():
        initialize_field(state_gt4py[key], field_array)
    
    return state_gt4py


def draw_fields(component: ComputationalGridComponent) -> NDArrayLikeDict:
    """Draw random fields according to component description

    Args:
        component (ComputationalGridComponent): a ComputationalGridComponent, with
        well described fields

    Returns:
        NDArrayLikeDict: dictionnary of random arrays associated with their field name
    """
    
    np.random.seed(23)
    
    if isinstance(component.array_shape, tuple):
        return {
        **{
            key: np.array(
                np.random.rand(*component.array_shape),
                dtype=float,
                order="F",
            )
            for key in component.fields_in.keys()
        },
        **{
            key: np.array(
                np.random.rand(*component.array_shape),
                float,
                order="F",
            )
            for key in component.fields_inout.keys()
        },
        **{
            key: np.array(
                np.random.rand(*component.array_shape),
                float,
                order="F",
            )
            for key in component.fields_out.keys()
        },
    }

    else:
        return TypeError("Array shape is not a tuple")
        


def compare_output(component, fortran_fields: dict, gt4py_state: dict):
    """Compare fortran and gt4py field mean on inout and out fields for a TestComponent

    Args:
        fortran_fields (dict): output fields from fortran
        gt4py_state (dict): output fields from gt4py
    """
    
    
    
    
        
    fortran_fields = unpack(fortran_fields, component.computational_grid)
    gt4py_fields = remove_y_axis(gt4py_state)
    
    absolute_differences = dict()
    fields_to_compare = {**component.fields_inout, **component.fields_out}
    for field_name, field_attributes in fields_to_compare.items():
        logging.info(f"Field to compare : {field_name}")
        logging.info(f"fortran field shape : {fortran_fields[field_name].shape}")
        logging.info(f"fortran field mean : {fortran_fields[field_name].mean()}")   
             
        # Removing njt dimension
        assert(gt4py_fields[field_name].shape == fortran_fields[field_name].shape)
        logging.info(f"gt4py field shape {gt4py_fields[field_name].shape}")
        logging.info(f"gt4py field mean : {gt4py_fields[field_name].values.mean()}")
        absolute_diff = abs(gt4py_fields[field_name] - fortran_fields[field_name]).values.mean()
        logging.info(f"{field_name}, absolute mean difference {absolute_diff}\n")
        absolute_differences.update({
            field_name: absolute_diff
        })
        
    return absolute_differences


def compare_input(component, fortran_fields: dict, gt4py_state: dict):
    """Compare fortran and gt4py field mean on inout and out fields for a TestComponent

    Args:
        fortran_fields (dict): output fields from fortran
        gt4py_state (dict): output fields from gt4py
    """
    absolute_differences = dict()
    fields_to_compare = {**component.fields_in, **component.fields_inout}
    logging.info(f"Input fields to compare {fields_to_compare.keys()}")
    logging.info(f"Fortran field keys {fortran_fields.keys()}")
    for field_name, field_attributes in fields_to_compare.items():
        logging.info(f"Input field name : {field_name}")
        fortran_name = field_attributes["fortran_name"]
        logging.info(f"(Input) fortran name {fortran_name}")
        fortran_field = fortran_fields[field_name][:, np.newaxis,:]
        logging.info(f"(Input) fortran field shape {fortran_field.shape}")
        logging.info(f"(Input) fortran field mean : {fortran_field.mean()}")
        
        # 2D fields + removing shadow level
        gt4py_field = gt4py_state[field_name]
        logging.info(f"gt4py field shape {gt4py_field.shape}")
        logging.info(f"gt4py field mean : {gt4py_field.values.mean()}")
        absolute_diff = abs(gt4py_field - fortran_field).values.mean()
        logging.info(f"{field_name}, absolute mean difference {absolute_diff}\n")
        absolute_differences.update({
            field_name: absolute_diff
        })
    
    for field, diff in absolute_differences.items():
            logging.info(f"Field name : {field}, error on input field : {diff}")
    

def run_test(component: ComputationalGridComponent):
    """Draw random arrays and call gt4py and fortran stencils side-by-side

    Args:
        component (ComputationalGridComponent): component to test
    """

    logging.info(f"\n Start test {component.__class__.__name__}")
    fields = draw_fields(component)
    state_gt4py = allocate_gt4py_fields(component, fields)
    
    logging.info(f"Compare input  fields")
    compare_input(component, fields, state_gt4py)
    
    logging.info("Calling fortran field")
    fortran_output_fields = component.call_fortran_stencil(fields)
    logging.info(f"fortran_output_fields {fortran_output_fields.keys()}")

    logging.info("Calling gt4py field")
    gt4py_output_fields = component.call_gt4py_stencil(state_gt4py)
    
    # TODO: remove shadow first level

    logging.info("Compare output fields")
    absolute_differences = compare_output(component=component, fortran_fields=fortran_output_fields, gt4py_state=gt4py_output_fields)
    logging.info(f"End test {component.__class__.__name__}\n")
    
    return absolute_differences
    
    