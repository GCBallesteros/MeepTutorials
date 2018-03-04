# Mode Profile Frequency Monitors

> This blog post is still under preparation. It may or may not work :)

This is the code for the frequency monitor blog. It shows how one can create a monitor that will show the fields at all frequencies running a single simulation.

## Running
Make *show_freq_monitor.py* and *add_attr.py* executable.

1. Run the simulation
  ```
  ./run_simulation.sh
  ```

2. To show the interactive frequency monitor:
  ```
  ./show_freq_monitor --fname ring_resonator-freq-monitor.h5 --npoints 200
  ```


You can view the different options available for the frequency monitor viewer by doing:
```
./show_freq_monitor.py --help to all the options available.
```


## Extra Content

The notebook in this folder will plot the LDOS which is output after running the simulations.

For more details plese visit *maxwellrules.com*.
