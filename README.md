# ASIC event generator
This is an event generator for an ASIC, developed for the PANDA experiment. I wrote this, because I needed some input for my simulations and check if the logic of the ASIC is doing, what it is supposed to do. Therefore there are many specifics to the chip I'm working on, but maybe you can use some parts of this.

### Concept
This script emulates the voltage input of two discriminators (time and energy). They have different threshold, time branch lower than the energy branch. The script produces voltages for the given range of time in discrete steps. Randomly, an event is generated and the according voltage increase is stored in this steps.

Finally, for all stored voltages it is checked if they are above the given threshold levels. Discriminator output `DOT = True` means voltage higher than threshold for time branch, `DOT = False` means lower than threshold. If the DOT/DOE values change, the time since the last flip (in ps) is saved to a file. Later in the simulation, these times will be used as wait information between the change of the DOT/DOE states.

## Installation / Usage
No special requirements. Python > 2.6 is needed, after that you can run the program with

    python generate_asic_data.py <time> [<threads>]

* **`time`**: Run time in seconds.
* **`threads`**: (Optional) Number of threads, 0 is number of CPU cores. *Default: 0*

The file creates a directory `asic_data_?ms` with the run time in ms. Inside the directory DOT/DOE switching information (see concept) is stored for each channel in the format `ch?_DOx.dat`.

### Configuration
Inside the python document some more parameters are given to configure the emulation process. See the `_conf` parameter and its comments for details.
