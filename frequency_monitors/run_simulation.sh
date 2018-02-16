#!/bin/bash
meep ring_resonator.ctl > out
grep ldos out > ldos
h5topng ring_resonator-eps-*.h5