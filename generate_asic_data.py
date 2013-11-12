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
#  Description:
#
#  Revisions:
#    1.0 Initial revision
#
# ------------------------------------------------------------------------------

# Import stuff
import sys          # system functions (like exit)
import time         # time functions
import re           # regular expressions
import os           # some useful functions for files
import random       # for pseudo-random numbers
import math         # you know, the math magic
import array        # we store the data in arrays


# default values
_runTime = 0.001
_nThreads = 0


# other global values
_nChannels = 64
_conf = {
    "eventRate": 160e3, # Hz
    "darkRate": 1e6,    # Hz
    "tThreshold": 0.5,  # mV
    "eThreshold": 10,   # mV

    "tPeakT": 1e-9,
    "alphaT": 1.0,
    "tPeakE": 4.0,
    "alphaE": 0.5,

    "maxToT": 400e-9,   # s
    "maxD1": 1.5e-9,    # s
    "maxD2": 400e-9,    # s

    "stepping": 1e-11,  # s

    "enableSaving": True,
    "directory": "asic_data_{time:.0f}ms/",  # include trailing slash!
    "filename": {
        "stats": "ch{channel:d}_DOT.dat",
        "trues": "ch{channel:d}_trues.dat",
        "dark": "ch{channel:d}_dark.dat",
        "DOT": "ch{channel:d}_DOT.dat",
        "DOE": "ch{channel:d}_DOE.dat"
    }
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
        if 0 <= i < self.nSteps:
            oldValue = self.getValue(i)
            if not oldValue:
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

            if self._data[iBlock] == None:
                return 0.0
            elif len(self._data[iBlock]) == self._blockSize:
                return self._data[iBlock][iSlot]
            else:
                return False
        else:
            return False

    def getSteps(self):
        return self._nSteps

    def incrEventCounter(self):
        self._nEvents += 1

    def getEvents(self):
        return self._nEvents

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



def addEvent(signal, t0, height):
    #           x ** (alpha - 1) * math.exp(-x / beta)
    # pdf(x) =  --------------------------------------
    #             math.gamma(alpha) * beta ** alpha
    iBin0 = long(t0 / _conf['stepping'])
    iBin = iBin0

    while True:
        t = float(iBin - iBin0) * _conf['stepping']    # t in s
        x = t * 1e9 / height
        currentValue = height * 25 * x**2 * math.exp(-x)

        # if it couldn't be set, it is probably out of range
        if not signal.setValue(iBin, currentValue):
            break

        iBin += 1
        # regular stop criterium
        if t > 1e-9 and currentValue < _conf['tThreshold']:
            break

        # pulses longer than 1 us don't make sense
        if (t - t0) > 1e-6:
            break

    signal.incrEventCounter()
    return t


def generateTrues(thread, signal):
    channel = thread.status.getChannel()
    filenameTrues = _conf['directory'] + _conf['filename']['trues'].format(channel=int(channel))

    if _conf['enableSaving']:
        fileTrues = file(filenameTrues, 'w')

    t = float(0)
    while True:
        t_wait = random.gammavariate(1.0, 1/_conf['eventRate'])

        progress = int(round(t / _runTime * 100))
        if progress > thread.status.getProgress():
            thread.status.setProgress(progress)
            time.sleep(0.001)   # give the output some time to print

        # time of the next event
        t_next = t + t_wait
        if t_next > _runTime:
            break

        # generate pulse
        height = random.uniform(1.0, 30.0)
        #print "generate event at t={0} with heigth {1}".format(t,height)
        length = addEvent(signal, t_next, height)

        # write to file
        if _conf['enableSaving']:
            fileTrues.write("{0:.12f} {1:.12f}\n".format(t, length))

        # time for the next event
        if length > t_wait:
            t += length # fix pile-up
        t += t_wait

    if _conf['enableSaving']:
        fileTrues.close()


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
def generateChannelEvents(thread, channel):
    thread.status.setChannel(channel)

    # init the signal objects
    nSteps = long((_runTime + 1e-6)/_conf['stepping'])
    discInput = Signal(nSteps)

    # generate the true hits
    thread.status.setStep(1)
    thread.status.setProgress(0)
    generateTrues(thread, discInput)

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
    print "Usage: generate_asic_data.py <time> [<cores>]"
    print ""
    print "Parameter:"
    print "  time       The simulated time of events in seconds."
    print "  cores      Number of CPU cores to use. Default: max"
    print ""
    sys.exit()


# the main function
def main():
    global _conf, _runTime, _nThreads

    if 2 <= len(sys.argv) <= 3:
        if len(sys.argv) >= 3:
            if int(sys.argv[2]) > 0:
                _nThreads = int(sys.argv[2])

        _runTime = float(sys.argv[1])
    else:
        printUsage()
        sys.exit()

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
        pool.add_task(generateChannelEvents, ch)
        # printStatus(pool)

    # 3) Wait for completion
    pool.wait_completion()
    printStatus(pool)

    print "\n\nEverything done!"


if __name__ == '__main__':
    main()
