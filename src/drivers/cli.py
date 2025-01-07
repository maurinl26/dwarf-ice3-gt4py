# -*- coding: utf-8 -*-
from pathlib import Path
import subprocess

import numpy as np
from ice3_gt4py.components.ice_adjust_split import IceAdjustSplit
import typer
import logging
import datetime
import time
import sys
import xarray as xr

from ifs_physics_common.framework.config import GT4PyConfig
from ifs_physics_common.framework.grid import ComputationalGrid

from drivers.compare import compare_fields
from drivers.core import write_dataset, write_performance_tracking
from ice3_gt4py.components.ice_adjust import IceAdjust
from ice3_gt4py.components.rain_ice import RainIce
from ice3_gt4py.initialisation.state_ice_adjust import (
    get_state_ice_adjust,
)
from ice3_gt4py.initialisation.state_rain_ice import get_state_rain_ice
from ice3_gt4py.phyex_common.phyex import Phyex
from ice3_gt4py.utils.reader import NetCDFReader

logging.basicConfig(level=logging.DEBUG, stream=sys.stdout)
logging.getLogger()

app = typer.Typer()

######################## GT4Py drivers #######################
@app.command()
def run_ice_adjust(
    backend: str,
    dataset: str,
    output_path: str,
    tracking_file: str,
    rebuild: bool = True,
    validate_args: bool = False,
):
    """Run ice_adjust component"""

    ##### Grid #####
    logging.info("Initializing grid ...")
    nx = 10000
    ny = 1
    nz = 15
    grid = ComputationalGrid(nx, ny, nz)
    dt = datetime.timedelta(seconds=1)

    ################## Phyex #################
    logging.info("Initializing Phyex ...")
    phyex = Phyex("AROME")

    ######## Backend and gt4py config #######
    logging.info(f"With backend {backend}")
    gt4py_config = GT4PyConfig(
        backend=backend, rebuild=rebuild, validate_args=validate_args, verbose=True
    )

    ######## Instanciation + compilation #####
    logging.info(f"Compilation for IceAdjust stencils")
    start = time.time()
    ice_adjust = IceAdjust(grid, gt4py_config, phyex)
    stop = time.time()
    elapsed_time = stop - start
    logging.info(f"Compilation duration for IceAdjust : {elapsed_time} s")

    ####### Create state for AroAdjust #######
    reader = NetCDFReader(Path(dataset))

    logging.info("Getting state for IceAdjust")
    state = get_state_ice_adjust(grid, gt4py_config=gt4py_config, netcdf_reader=reader)
    logging.info(f"Keys : {list(state.keys())}")

    ###### Launching IceAdjust ###############
    logging.info("Launching IceAdjust")

    # TODO: decorator for tracking
    start = time.time()
    tends, diags = ice_adjust(state, dt)
    stop = time.time()
    elapsed_time = stop - start
    logging.info(f"Execution duration for IceAdjust : {elapsed_time} s")

    #################### Write dataset ######################
    write_dataset(state, (nx, ny, nz), output_path)

    ############### Compute differences per field ###########
    metrics = compare_fields(dataset, output_path, "ice_adjust")

    ####################### Tracking ########################
    write_performance_tracking(gt4py_config, metrics, tracking_file)


@app.command()
def run_ice_adjust_split(
    backend: str,
    dataset: str,
    output_path: str,
    tracking_file: str,
    rebuild: bool = True,
    validate_args: bool = False,
):
    """Run ice_adjust splitted version to avoid
    interpolation problems for sigrc
    """

    ##### Grid #####
    logging.info("Initializing grid ...")
    nx = 9472
    ny = 1
    nz = 15
    grid = ComputationalGrid(nx, ny, nz)
    dt = datetime.timedelta(seconds=50)

    ################## Phyex #################
    logging.info("Initializing Phyex ...")
    phyex = Phyex("AROME")

    ######## Backend and gt4py config #######
    logging.info(f"With backend {backend}")
    gt4py_config = GT4PyConfig(
        backend=backend, rebuild=rebuild, validate_args=validate_args, verbose=True
    )

    ######## Instanciation + compilation #####
    logging.info(f"Compilation for IceAdjust stencils")
    start_compilation = time.time()
    ice_adjust = IceAdjustSplit(grid, gt4py_config, phyex)
    stop_compilation = time.time()
    elapsed_time = stop_compilation - start_compilation
    logging.info(f"Compilation duration for IceAdjust : {elapsed_time} s")

    ####### Create state for AroAdjust #######
    reader = NetCDFReader(Path(dataset))

    logging.info("Getting state for IceAdjust")
    state = get_state_ice_adjust(grid, gt4py_config=gt4py_config, netcdf_reader=reader)

    ####### Check inputs #####################
    logging.info("Fields IN")
    field = state["rhodref"].isel(z=slice(0, 15))
    logging.info(f"Field rhodref, shape {field.shape}")
    logging.info(f"Field rhodref, mean {field.mean().values}, std {field.std().values}")

    logging.info("Fields INOUT (before call)")
    logging.info(
        f"Field rvs, mean {state['rvs'].isel(z=slice(0,15)).mean().values}, std {state['rvs'].isel(z=slice(0,15)).std().values}"
    )
    logging.info(
        f"Field rcs, mean {state['rcs'].isel(z=slice(0,15)).mean().values}, std {state['rcs'].isel(z=slice(0,15)).std().values}"
    )
    logging.info(
        f"Field ris, mean {state['ris'].isel(z=slice(0,15)).mean().values}, std {state['ris'].isel(z=slice(0,15)).std().values}"
    )
    logging.info(
        f"Field ths, mean {state['ths'].isel(z=slice(0,15)).mean().values}, std {state['ths'].isel(z=slice(0,15)).std().values}"
    )

    ###### Launching IceAdjust ###############
    logging.info("Launching IceAdjust")

    # TODO: decorator for tracking
    start = time.time()
    tends, diags = ice_adjust(state, dt)
    stop = time.time()
    elapsed_time = stop - start
    logging.info(f"Execution duration for IceAdjust : {elapsed_time} s")

    logging.info("Fields INOUT (after call)")
    logging.info(
        f"Field rvs, mean {state['rvs'].isel(z=slice(0,15)).mean().values}, std {state['rvs'].isel(z=slice(0,15)).std().values}"
    )
    logging.info(
        f"Field rcs, mean {state['rcs'].isel(z=slice(0,15)).mean().values}, std {state['rcs'].isel(z=slice(0,15)).std().values}"
    )
    logging.info(
        f"Field ris, mean {state['ris'].isel(z=slice(0,15)).mean().values}, std {state['ris'].isel(z=slice(0,15)).std().values}"
    )
    logging.info(
        f"Field ths, mean {state['ths'].isel(z=slice(0,15)).mean().values}, std {state['ths'].isel(z=slice(0,15)).std().values}"
    )

    logging.info(
        f"Field hlc_hrc, mean {state['hlc_hrc'].isel(z=slice(0,15)).mean().values}, std {state['hlc_hrc'].isel(z=slice(0,15)).std().values}"
    )
    logging.info(
        f"Field hlc_hcf, mean {state['hlc_hcf'].isel(z=slice(0,15)).mean().values}, std {state['hlc_hcf'].isel(z=slice(0,15)).std().values}"
    )
    logging.info(
        f"Field hli_hcf, mean {state['hli_hcf'].isel(z=slice(0,15)).mean().values}, std {state['hli_hcf'].isel(z=slice(0,15)).std().values}"
    )
    logging.info(
        f"Field hli_hri, mean {state['hli_hri'].isel(z=slice(0,15)).mean().values}, std {state['hli_hri'].isel(z=slice(0,15)).std().values}"
    )

    #################### Write dataset ######################
    write_dataset(state, (nx, ny, nz), output_path)

    ############### Compute differences per field ###########
    # metrics = compare_fields(dataset, output_path, "ice_adjust")

    # ####################### Tracking ########################
    # write_performance_tracking(gt4py_config, metrics, tracking_file)


@app.command()
def run_rain_ice(
    backend: str,
    dataset: str,
    output_path: str,
    tracking_file: str,
    rebuild: bool = True,
    validate_args: bool = False,
):
    """Run aro_rain_ice component"""

    ##### Grid #####
    logging.info("Initializing grid ...")
    nx = 5000
    ny = 1
    nz = 15
    grid = ComputationalGrid(nx, ny, nz)
    dt = datetime.timedelta(seconds=1)

    ################## Phyex #################
    logging.info("Initializing Phyex ...")
    cprogram = "AROME"
    phyex = Phyex(cprogram)

    ######## Backend and gt4py config #######
    logging.info(f"With backend {backend}")
    gt4py_config = GT4PyConfig(
        backend=backend, rebuild=rebuild, validate_args=validate_args, verbose=True
    )

    ######## Instanciation + compilation #####
    logging.info(f"Compilation for RainIce stencils")
    start = time.time()
    rain_ice = RainIce(grid, gt4py_config, phyex)
    stop = time.time()
    elapsed_time = stop - start
    logging.info(f"Compilation duration for RainIce : {elapsed_time} s")

    ####### Create state for AroAdjust #######
    reader = NetCDFReader(Path(dataset))

    logging.info("Getting state for RainIce")
    state = get_state_rain_ice(grid, gt4py_config=gt4py_config, netcdf_reader=reader)
    logging.info(f"Keys : {list(state.keys())}")

    ###### Launching RainIce ###############
    logging.info("Launching RainIce")

    start = time.time()
    tends, diags = rain_ice(state, dt)
    stop = time.time()
    elapsed_time = stop - start
    logging.info(f"Execution duration for RainIce : {elapsed_time} s")

    logging.info(f"Extracting state data to {output_path}")
    output_fields = xr.Dataset(state)
    for key, field in state.items():
        if key not in ["time"]:
            if key not in ["sea", "town", "inprr", "inprc", "inprg", "inprs"]:
                array = xr.DataArray(
                    data=field.data[:, :, 1:],
                    dims=["I", "J", "K"],
                    coords={
                        "I": range(nx),
                        "J": range(ny),
                        "K": range(nz),
                    },
                    name=f"{key}",
                )
                output_fields[key] = array
            else:
                array = xr.DataArray(
                    data=field.data,
                    dims=["I", "J"],
                    coords={
                        "I": range(nx),
                        "J": range(ny),
                    },
                    name=f"{key}",
                )
                output_fields[key] = array
    output_fields.to_netcdf(Path(output_path))

    ################## Metrics ###################################
    metrics = compare_fields(dataset, output_path, "rain_ice")

    ################### Performance tracking ######################
    write_performance_tracking(gt4py_config, metrics, tracking_file)


##################### Fortran drivers #########################
@app.command()
def run_ice_adjust_fortran(
    testdir: str, name: str, archfile: str, checkOpt: str, extrapolation_opts: str
):
    """Call and run main_ice_adjust (Fortran)"""

    try:

        logging.info("Setting env variables")
        subprocess.run(
            ["source", f"{testdir}/{name}/build/with_fcm/arch_{archfile}/arch.env"]
        )

        logging.info("Job submit")
        subprocess.run(
            [
                "submit",
                "Output_run",
                "Stderr_run",
                f"{testdir}/{name}/build/with_fcm/arch_${archfile}/build/bin/main_ice_adjust.exe",
                f"{checkOpt}",
                f"{extrapolation_opts}",
            ]
        )

    except RuntimeError as e:
        logging.error("Fortran ice_adjust execution failed")
        logging.error(f"{e}")


@app.command()
def run_rain_ice_fortran(
    testdir: str, name: str, archfile: str, checkOpt: str, extrapolation_opts: str
):
    """Call and run main_rain_ice (Fortran)"""
    pass


if __name__ == "__main__":
    app()
