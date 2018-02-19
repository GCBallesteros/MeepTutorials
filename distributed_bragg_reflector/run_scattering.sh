#!/bin/bash
meep structure?=false scattering_dbr.ctl > out && grep flux1 out > reference_power;
meep structure?=true  scattering_dbr.ctl > out && grep flux1 out > structure_power;
h5topng scattering_dbr-eps*.h5
