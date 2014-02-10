#!/usr/bin/python
# -*- coding: utf-8

# ------------------------------------------------------------------------------
#
#    Generate some input data for the ASIC's simulation
#   ------------------------------------------------------
#
#  Author:  André Goerres (a.goerres@fz-juelich.de)
#  Created: 2013-11-08
#
#  Usage:  python generate_asic_data.py <time [ms]> [<number_of_threads>]
#
#    <time [ms]>            The time span to generate in milliseconds.
#    <number_of_threads>    (Optional) Number of threads to use.
#                           Default: number of CPUs
#
#  Description:  With this script you can generate yourself some input for the
#    digital simulation of the TOFPET/PASTA ASICs. The generator creates event
#    and dark count pulses with the given parameters. Using the resulting course
#    of voltage, the discriminator decision is simulated. For every change in
#    the discriminator's output, the time difference since the last change is
#    saved in a file (DOT/DOE files). Lateron in the simulation, these time
#    values are used as wait times, until the discriminator value should be
#    changed.
#
# ------------------------------------------------------------------------------

# Import stuff
import sys          # system functions (like exit)
import time         # time functions
import os           # some useful functions for files
import random       # for pseudo-random numbers
import math         # you know, the math magic
import array        # we store the data in arrays
import csv          # parse the input file


# default values
_runTime = 0.001
_nThreads = 0


# other global values
_nChannels = 64
_conf = {
    # simulation parameters
    "stepping": 5e-11,   # s

    "eventRate": {
        "trues": 160e3,  # Hz
        "dark":    1e6,  # Hz
    },
    "tThreshold":  1e-3, # V
    "eThreshold": 10e-3, # V

    "chargeDistAlpha": 2, # shape of the distribution (gamma distribution)
    "chargeDistBeta":  5, # maximum of distribution (in fC)

    "darkPulseHeight": 4e-3, # V

    # pulse shape
    "useSimulationInput": True,  # setting False uses comb. of x^2 and exp(-x)
    "simulationInputFiles": {
        "neg": "input_pulse/NegativePolarity.csv",
        "pos": "input_pulse/PositivePolarity.csv",
        "minCharge": 1,     # fC
        "stepping": 1,      # fC
    },
    "polarity": "neg",      # possible values: pos, neg

    # output parameters
    "enableSaving": True,   # disable saving the output, more for debugging
    "directory": "asic_data_{time:.0f}ms/",  # include trailing slash!
    "filename": {
        "stats": "ch{channel:d}_DOT.dat",
        "trues": "ch{channel:d}_trues.dat",
        "dark": "ch{channel:d}_dark.dat",
        "DOT": "ch{channel:d}_DOT.dat",
        "DOE": "ch{channel:d}_DOE.dat"
    },
}


# thread pool stuff
import multiprocessing as mp
_maxThreads = mp.cpu_count()

class Worker(mp.Process):
    """Thread executing tasks from a given tasks queue"""
    def __init__(self, tasks, id):
        mp.Process.__init__(self)
        self.tasks = tasks
        self.daemon = True
        self.name = "{0}".format(id)
        self.status = WorkerStatus()
        self.start()

    def run(self):
        while True:
            func, args, kargs = self.tasks.get()
            try: func(self, *args, **kargs)
            except Exception, e: print e
            self.tasks.task_done()

class WorkerStatus():
    """Manage the worker's status in shared memory values"""
    def __init__(self):
        self.channel = mp.Value('i', -1)
        self.step = mp.Value('i', 0)
        self.progress = mp.Value('i', 0)
        self.lock = mp.Lock()

    def setChannel(self, c):
        with self.lock:
            self.channel.value = c

    def setStep(self, s):
        with self.lock:
            self.step.value = s

    def setProgress(self, p):
        with self.lock:
            self.progress.value = p

    def getChannel(self):
        return self.channel.value

    def getStep(self):
        return self.step.value

    def getProgress(self):
        return self.progress.value


class ThreadPool:
    """Pool of threads consuming tasks from a queue"""
    def __init__(self, num_threads):
        self.tasks = mp.JoinableQueue(num_threads)
        self.workers = []
        for i in range(num_threads):
            self.workers.append( Worker(self.tasks, i) )

    def add_task(self, func, *args, **kargs):
        """Add a task to the queue"""
        self.tasks.put((func, args, kargs))

    def wait_completion(self):
        """Wait for completion of all the tasks in the queue"""
        self.tasks.join()
# thread pool stuff end


# The signal with shape
#
# To save memory, not every point is saved. The list of data is splitted into
# blocks, which are empty on start. Only if an event is being created, a block
# is filled with an array of data points and the values are saved.
class Signal:
    def __init__(self, nSteps):
        self._blockSize = 2048
        self._nSteps = long(nSteps)
        self._nBlocks = long(nSteps / self._blockSize) + 1
        self._data = [None]*self._nBlocks
        self._nEvents = 0

    # initialize a block with zeros
    def createPage(self, iBlock):
        self._data[iBlock] = array.array('f', (0,)*self._blockSize)

    # set a certain data point to a value
    def setValue(self, i, value):
        if 0 <= i < self._nSteps:
            iBlock = int(i / self._blockSize)
            iSlot = i % self._blockSize   # slot inside block

            if self._data[iBlock] == None:
                self.createPage(iBlock)

            self._data[iBlock][iSlot] = value
            return True
        else:
            return False

    # add a value to a certain data point
    def addValue(self, i, value):
        if 0 <= i < self._nSteps:
            oldValue = self.getValue(i)
            if oldValue is None:
                return False
            else:
                return self.setValue(i, value + oldValue)
        else:
            return False

    # get the value of a certain data point
    def getValue(self, i):
        if 0 <= i < self._nSteps:
            iBlock = int(i / self._blockSize)
            iSlot = i % self._blockSize   # slot inside block

            if self._data[iBlock] is None:
                return 0.0
            elif len(self._data[iBlock]) == self._blockSize:
                return self._data[iBlock][iSlot]
            else:
                return None
        else:
            return None

    def getSteps(self):
        return self._nSteps

    def incrEventCounter(self):
        self._nEvents += 1

    def getEvents(self):
        return self._nEvents

    # Search the next event from a given point onwards. Because the separation
    # of the data into blocks is only known internally in the class, this helps
    # the user to not have to loop in small steps over a long period of
    # nothingness.
    def searchNextFrom(self, i, thr):
        if 0 <= i < self._nSteps:
            iBlock = long(i / self._blockSize)
            iSlot = i % self._blockSize   # slot inside block

            while iBlock < self._nBlocks:
                # there is nothing here to see, go on to next block
                if self._data[iBlock] == None:
                    iBlock += 1
                    iSlot = 0
                    continue

                # block has data, lets see if there is something interesting
                if len(self._data[iBlock]) == self._blockSize:
                    if self._data[iBlock][iSlot] >= thr:
                        return iBlock * self._blockSize + iSlot

                # next patient, please
                if iSlot < self._blockSize-1:
                    iSlot += 1
                else:
                    iBlock += 1
                    iSlot = 0

            # nothing found
            return False
        else:
            return False


class PulseShape:
    def __init__(self, useSimInput, inputFile):
        self.useSimInput = useSimInput
        self.shape = {}
        self.cutoff = 0
        self.minCharge = _conf['simulationInputFiles']['minCharge']
        self.maxCharge = _conf['simulationInputFiles']['minCharge']

        if useSimInput:
            self.loadCSV(inputFile)

    def loadCSV(self, csvFileName):
        with open(csvFileName, 'rb') as csvFile:
            reader = csv.reader(csvFile, delimiter=',')
            firstrow = True
            for row in reader:
                # skip first row (header)
                if firstrow:
                    firstrow = False
                    continue

                for i in range(len(row)/2):
                    if len(row[2*i].strip()) == 0 or len(row[2*i+1].strip()) == 0:
                        continue
                    t = float(row[2*i])
                    V = float(row[2*i+1])
                    if self.cutoff == 0:
                        self.cutoff = V
                        continue
                    elif V <= self.cutoff:
                        continue
                    V_corr = V - self.cutoff

                    curCharge = _conf['simulationInputFiles']['minCharge'] + i*_conf['simulationInputFiles']['stepping']
                    if curCharge > self.maxCharge:
                        self.maxCharge = curCharge

                    if not curCharge in self.shape:
                        self.shape[curCharge] = { 't': [], 'V': [], 't-offset': t }
                        t_rel = 0.0
                    else:
                        t_rel = t - self.shape[curCharge]['t-offset']

                    self.shape[curCharge]['t'].append(t_rel)
                    self.shape[curCharge]['V'].append(V_corr)
            self.datapoints = sum([len(v['t']) for k,v in self.shape.iteritems()])

    def selectPulse(self, charge):
        charge = round(charge)
        if charge < self.minCharge:
            charge = self.minCharge
        elif charge > self.maxCharge:
            charge = self.maxCharge
        if not charge in self.shape:
            return None

        return self.shape[charge]

    # generate pulses based on simulation output files
    def addEvent(self, signal, charge, t0):
        pulse = self.selectPulse(charge)
        if pulse is None:
            return

        iBin0 = long(t0 / _conf['stepping'])
        iBin = iBin0
        simTimeNextIndex = 0
        lastEstimate = 0
        t_thr = None

        while True:
            t = float(iBin - iBin0) * _conf['stepping']    # t in s
            for i in range(simTimeNextIndex, len(pulse['t'])):
                if pulse['t'][i] > t:
                    simTimeNextIndex = i
                    simTimePrevIndex = i
                    simTimeNext = pulse['t'][i]
                    simTimePrev = pulse['t'][i-1]
                    break

            # do some linear interpolation to get a better estimation of the voltage
            ratio = (t - simTimePrev) / (simTimeNext - simTimePrev)
            voltageEstimate = (pulse['V'][simTimeNextIndex] - pulse['V'][simTimePrevIndex]) * ratio + pulse['V'][simTimePrevIndex]

            # if it couldn't be set, it is probably out of range
            if not signal.addValue(iBin, voltageEstimate):
                break

            # get the time, when the threshold has been crossed
            if t_thr == None and lastEstimate < _conf['tThreshold'] <= voltageEstimate:
                # t_last = t - _conf['stepping']
                # voltageRatio = (_conf['tThreshold'] - lastEstimate) / (voltageEstimate - lastEstimate)
                # t_thr = t0 + t_last + voltageRatio * (t - t_last)
                t_thr = t0 + t
                # print "{0:.4f} {1:.4f} {2:.4f} {3:.2f}".format(lastEstimate, voltageEstimate, _conf['tThreshold'], voltageRatio)
                # print "{0:.1e} {1:.3e} {2:.1e}".format(t_last, t_last+voltageRatio * (t - t_last), t)
                # print

            iBin += 1
            lastEstimate = voltageEstimate

            # stop, when time is larger than pulse length
            if t >= pulse['t'][-1]:
                break

        if t_thr == None:
            t_thr = t0

        return t_thr, pulse['t'][-1]


def darkPulseValue(t_ns):
    return _conf['darkPulseHeight']*1.85 * t_ns**2 * math.exp(-t_ns)


# Generate the pulse shape of one dark event. The function is an overlay of
# a quadratic factor (x^2) and an exponential factor (exp(-x)).
def addDarkPulse(signal, t0):
    iBin0 = long(t0 / _conf['stepping'])
    iBin = iBin0
    t_thr = None

    while True:
        t = float(iBin - iBin0) * _conf['stepping']   # t in s
        t_ns = t * 1e9  # t in ns
        currentValue = darkPulseValue(t_ns)

        # get the time, when the threshold is crossed (linear interpolation)
        if t_thr == None and currentValue < _conf['tThreshold']:
            t_next = t_ns + _conf['stepping']
            nextValue = darkPulseValue(t_next*1e9)
            if nextValue > _conf['tThreshold']:
                # ratio = (_conf['tThreshold'] - currentValue) / (nextValue - currentValue)
                # t_thr = t0+t + ratio * (t_next - t)
                t_thr = t0 + t_next

        # if it couldn't be set, it is probably out of range
        if not signal.addValue(iBin, currentValue):
            break

        iBin += 1
        # stop after 20 ns
        if t_ns >= 20:
            break

    if t_thr == None:
        t_thr = t0

    signal.incrEventCounter()
    return t_thr, t


# Generate the event pulses (= trues) with a random interval in between two
# pulses. The parameter of one pulse is also chosen randomly.
def generateHits(thread, signal, pulse, eventType):
    if eventType != 'trues' and eventType != 'dark':
        return

    channel = thread.status.getChannel()
    fileName = _conf['directory'] + _conf['filename'][eventType].format(channel=int(channel))

    if _conf['enableSaving']:
        f = file(fileName, 'w')

    t = float(0)
    while True:
        saveThisEvent = True
        t_wait = random.gammavariate(1.0, 1/_conf['eventRate'][eventType])

        progress = int(round(t / _runTime * 100))
        if progress > thread.status.getProgress():
            thread.status.setProgress(progress)
            time.sleep(0.001)   # give the output some time to print

        # time of the next event
        t_next = t + t_wait
        if t_next > _runTime:
            break

        # generate pulse (maximum at 5 fC (=MIP))
        if eventType == 'trues':
            fC = random.gammavariate(_conf['chargeDistAlpha'], _conf['chargeDistBeta'])
            t_thr, length = pulse.addEvent(signal, fC, t_next)
            if fC < 2:
                saveThisEvent = False
        else:
            fC = 0.01;  # calculated integral, done with WolframAlpha (http://wolfr.am/1jJtxac)
            t_thr, length = addDarkPulse(signal, t_next)

        # write to file
        if _conf['enableSaving'] and saveThisEvent:
            f.write("{0:.12f} {1:06.3f}\n".format(t_thr, fC))

        # time for the next event
        if length > t_wait:
            t += length # fix pile-up
        t += t_wait

    if _conf['enableSaving']:
        f.close()


# Take course of voltage and generate the discriminator's output from that. As
# an output, only the time difference between changes of the DOx states are
# saved. These times are used lateron as wait times.
def getDiscriminatorOutput(thread, signal):
    channel = thread.status.getChannel()
    filenameDOT = _conf['directory'] + _conf['filename']['DOT'].format(channel=int(channel))
    filenameDOE = _conf['directory'] + _conf['filename']['DOE'].format(channel=int(channel))

    if _conf['enableSaving']:
        fileDOT = file(filenameDOT, 'w')
        fileDOE = file(filenameDOE, 'w')

    i = long(0)
    t_lastToggleDOT = 0
    t_lastToggleDOE = 0
    lastDOT = False
    lastDOE = False

    while i < signal.getSteps():
        progress = int(round(float(i) / signal.getSteps() * 100))
        if progress > thread.status.getProgress():
            thread.status.setProgress(progress)
            time.sleep(0.001)   # give the output some time to print

        t_ps = long(round(i * _conf['stepping'] * 1e12))   # t in ps
        voltage = signal.getValue(i)
        DOT = voltage >= _conf['tThreshold']
        DOE = voltage >= _conf['eThreshold']

        # DOx toggled
        if lastDOT != DOT:
            if _conf['enableSaving']:
                fileDOT.write("{0:012d}\n".format(t_ps - t_lastToggleDOT))
            t_lastToggleDOT = t_ps
        if lastDOE != DOE:
            if _conf['enableSaving']:
                fileDOE.write("{0:012d}\n".format(t_ps - t_lastToggleDOE))
            t_lastToggleDOE = t_ps

        lastDOT = DOT
        lastDOE = DOE

        # event is over, go to next event
        if not DOT and not DOE:
            i = signal.searchNextFrom(i, _conf['tThreshold'])

            # nothing found, stop searching
            if not i:
                break
        else:
            i += 1

    if _conf['enableSaving']:
        fileDOT.close()
        fileDOE.close()


# the actual generation of the events
def generateChannelEvents(thread, channel, pulse):
    thread.status.setChannel(channel)

    # init the signal objects
    nSteps = long((_runTime + 1e-6)/_conf['stepping'])
    discInput = Signal(nSteps)

    # generate the true hits
    thread.status.setStep(1)
    thread.status.setProgress(0)
    generateHits(thread, discInput, pulse, 'trues')

    # generate the dark hits
    thread.status.setStep(2)
    thread.status.setProgress(0)
    generateHits(thread, discInput, pulse, 'dark')

    # get the discriminator output
    thread.status.setStep(3)
    thread.status.setProgress(0)
    getDiscriminatorOutput(thread, discInput)
    thread.status.setProgress(100)

    # print "channel {0:2d}   generated Events: {1}".format(channel, discInput.getEvents())


# Keep the status output up to date
def updateStatus(pool):
    while True:
        printStatus(pool)
        time.sleep(0.03)


# Print some status information to the command line
def printStatus(pool, lastOutput=False):
    line = ""
    for w in pool.workers:
        if w.status.getChannel() < 0:
            line += " "*13 + "|"
        else:
            line += " ch{0:02d} {1:d} {2:3d}% |".format(w.status.getChannel(), w.status.getStep(), w.status.getProgress())

    if lastOutput:
        sys.stdout.write(line[:-1])
    else:
        sys.stdout.write(line[:-1] + '\r')
    sys.stdout.flush()


def printHeader():
    # some general description
    print "Generating events for {0:.0f} ms with a stepping of {1:.1e} s.\n".format(_runTime*1000, _conf['stepping'])
    print "Status output for each task:"
    print "  1) channel number (ch##)"
    print "  2) step in generation"
    print "     1: generate true events"
    print "     2: generate dark events"
    print "     3: locate DOT/DOE toggle points"
    print "  3) percentage of process in current step\n"
    print "Event generation starts...\n"

    # initialize the status output
    statusOutputHeadline = ""
    for t in range(_nThreads):
        statusOutputHeadline += "   task {0:2d}   |".format(t)
    print statusOutputHeadline[:-3]
    print (("-"*13 + "+")*_nThreads)[:-1]


# How the program is intended to use
def printUsage():
    print "Usage: generate_asic_data.py <time [ms]> [<threads>]"
    print ""
    print "Parameter:"
    print "  time       The simulated time of events in milliseconds."
    print "  threads    Number of CPU cores to use. Default: max"
    print ""
    sys.exit()


# the main function
def main():
    global _conf, _runTime, _nThreads

    # get the input data
    if 2 <= len(sys.argv) <= 3:
        if len(sys.argv) >= 3:
            if int(sys.argv[2]) > 0:
                _nThreads = int(sys.argv[2])

        _runTime = float(sys.argv[1]) / 1000.0
    else:
        printUsage()
        sys.exit()

    # stop the time
    startTime = time.time()

    # load the input pulse from an output of a simulation
    inputFile = _conf['simulationInputFiles'][_conf['polarity']]
    pulse = PulseShape(_conf['useSimulationInput'], inputFile)

    # set output directory based on runtime
    _conf['directory'] = _conf['directory'].format(time=_runTime*1000)
    if not os.path.exists(_conf['directory']) and _conf['enableSaving']:
        os.makedirs(_conf['directory'])

    # check number of threads
    if _nThreads == 0 or _nThreads > _maxThreads:
        _nThreads = _maxThreads

    # 1) Init a Thread pool with the desired number of threads
    pool = ThreadPool(_nThreads)

    # information for the user
    printHeader()
    statusThread = mp.Process(target=updateStatus, args=(pool,))
    statusThread.daemon = True
    statusThread.start()

    for ch in range(_nChannels):
        # 2) Add the task to the queue
        pool.add_task(generateChannelEvents, ch, pulse)
        printStatus(pool)

    # 3) Wait for completion
    pool.wait_completion()
    printStatus(pool)

    # stop the timer
    stopTime = time.time()
    totalTime = stopTime - startTime
    totalTimeSec = totalTime % 60
    totalTimeMin = int((totalTime - totalTimeSec) / 60)
    print "\n\nEvent generation finished after {0:02d}:{1:04.1f} min:sec.".format(totalTimeMin, totalTimeSec)
    print "Good bye and thanks for all the patience!"


if __name__ == '__main__':
    main()
