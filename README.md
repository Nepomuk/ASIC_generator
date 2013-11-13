# ASIC event generator
This is an event generator for an ASIC, developed for the PANDA experiment. I wrote this, because I needed some input for my simulations and check if the logic of the ASIC is doing, what it is supposed to do. Therefore there are many specifics to the chip I'm working on, but maybe you can use some parts of this.

### Concept
This script emulates the voltage input of two discriminators (time and energy). They have different threshold, time branch lower than the energy branch. The script produces voltages for the given range of time in discrete steps. Randomly, an event is generated and the according voltage increase is stored in this steps.

Finally, for all stored voltages it is checked if they are above the given threshold levels. Discriminator output `DOT = True` means voltage higher than threshold for time branch, `DOT = False` means lower than threshold. If the DOT/DOE values change, the time since the last flip (in ps) is saved to a file. Later in the simulation, these times will be used as wait information between the change of the DOT/DOE states.

## Installation / Usage
No special requirements. Python > 2.6 is needed, after that you can run the program with

    python generate_asic_data.py <time> [<threads>]

* **`time`**: Run time in milliseconds.
* **`threads`**: (Optional) Number of threads, 0 is number of CPU cores. *Default: 0*

The file creates a directory `asic_data_?ms` with the run time in ms. Inside the directory DOT/DOE switching information (see concept) is stored for each channel in the format `ch?_DOx.dat`, while the information about the generated pulses are stored in `ch?_trues.dat` and `ch?_dark.dat`.

### Configuration
Inside the python document, some parameters are given to configure the emulation process. See the `_conf` parameter and its comments for details. A selection is also explained here.

* **Simulation accuracy**: By changing the `stepping` parameter, you can change the step size for simulating the voltage values. Default is *1e-11 s*. Remeber: lower values cause a longer simulation time, the change is almost linearly.
* **Rates**: With `eventRate` and `darkRate` the rate of events/dark counts per second is defined. Default ist *160 kHz* of events and *1 MHz* of dark counts.
* **Thresholds**: The threshold of the discriminator is set by `tThreshold` and `eThreshold` for the time and energy branch, respectively. Default is *0.5 mV* and *10 mV*.
* **Output parameters**: There are some parameters defining the saving to files. You may change the directory (`directory`) and the filenames (`filename` list). For the directory, you should include a {time:.0f} tag that will be replaced by the run time in milliseconds. Filenames have a {channel:d} tag for the channel number, which should be included as well.

### Example output files
*Note*: Data for first channel, others obviously have the same structure.

#### ch0_trues.dat / ch0_dark.dat
Data structure:

    absolute_time pulse_length

Both values in seconds with 12 decimal places (equals to ps resolution).

Example data (6 events):

    0.000006520456 0.000000101280
    0.000030262739 0.000000104310
    0.000035721378 0.000000081840
    0.000045872770 0.000000205050
    0.000049860426 0.000000048890
    0.000056118203 0.000000061700


#### ch0_DOT.dat / ch0_DOE.dat
Data structure:

    time_since_last_DOT_change

Value in picoseconds as an integer, filled up with leading zeros to 12 digits.

Example data (6 events):

    000006520900
    000000100830
    000023641450
    000000103860
    000005354740
    000000081430
    000010070160
    000000204440
    000003782940
    000000048560
    000006209250
    000000061340

