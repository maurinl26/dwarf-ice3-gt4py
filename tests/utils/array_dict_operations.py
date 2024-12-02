from ifs_physics_common.framework.grid import I, J, K
import numpy as np

############# Utils #############
def absolute_differences(fortran_fields, gt4py_fields, fields_to_compare):
        return {
            abs(gt4py_fields[field_name] - fortran_fields[field_name]).values.mean()
            if (gt4py_fields[field_name].shape == fortran_fields[field_name].shape) else None
            for field_name in fields_to_compare.keys() 
        }
        
def remove_y_axis(fields):
        return {
            key: np.squeeze(array, axis=1) for key, array in fields.items()
        }
    
def unpack(fields, component_grid):
        """Unpack as a 2d field

        Args:
            fortran_fields (_type_): _description_
            component_grid (_type_): _description_

        Returns:
            _type_: _description_
        """
        
        nit, njt, nkt = component_grid.grids[(I, J, K)].shape
        return {
            key: array.reshape(nit*njt, nkt)  for key, array in fields.items()
        }