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

    "stepping": 10e-10, # s

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
from Queue import Queue
from threading import Thread

from multiprocessing import cpu_count
_maxThreads = cpu_count()

class Worker(Thread):
    """Thread executing tasks from a given tasks queue"""
    def __init__(self, tasks):
        Thread.__init__(self)
        self.tasks = tasks
        self.daemon = True
        self.start()

    def run(self):
        while True:
            func, args, kargs = self.tasks.get()
            try: func(*args, **kargs)
            except Exception, e: print e
            self.tasks.task_done()


class ThreadPool:
    """Pool of threads consuming tasks from a queue"""
    def __init__(self, num_threads):
        self.tasks = Queue(num_threads)
        for _ in range(num_threads): Worker(self.tasks)

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
        self._nSteps = nSteps
        self._nBlocks = int(nSteps / self._blockSize) + 1
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
            iBlock = int(i / self._blockSize)
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
    alpha = 3.0
    iBin0 = int(t0 / _conf['stepping'])
    iBin = iBin0

    while True:
        t = (iBin - iBin0) * _conf['stepping']    # t in s
        x = t * 1e9 / height
        currentValue = height * 25 * x**2 * math.exp(-x)

        # if it couldn't be set, it is probably out of range
        if not signal.setValue(iBin, currentValue):
            break

        iBin += 1
        # regular stop criterium
        if (iBin - iBin0) > 10 and currentValue < _conf['tThreshold']:
            break

        # pulses longer than 1 us don't make sense
        if (t - t0) > 1e-6:
            break

    signal.incrEventCounter()
    return t


def generateTrues(signal, runTime):
    global _conf
    t = 0
    while True:
        t_wait = random.gammavariate(1.0, 1/_conf['eventRate'])

        # time of the next event
        t_next = t + t_wait
        if t_next > runTime:
            break

        # generate pulse
        height = random.uniform(1.0, 30.0)
        #print "generate event at t={0} with heigth {1}".format(t,height)
        length = addEvent(signal, t_next, height)

        # time for the next event
        if length > t_wait:
            t += length # fix pile-up
        t += t_wait


def getDiscriminatorOutput(signal, channel):
    global _conf

    filenameDOT = _conf['directory'] + _conf['filename']['DOT'].format(channel=int(channel))
    filenameDOE = _conf['directory'] + _conf['filename']['DOE'].format(channel=int(channel))

    fileDOT = file(filenameDOT, 'w')
    fileDOE = file(filenameDOE, 'w')

    i = 0
    t_lastToggleDOT = 0
    t_lastToggleDOE = 0
    lastDOT = False
    lastDOE = False

    while i < signal.getSteps():
        t_ps = int(round(i * _conf['stepping'] * 1e12))   # t in ps
        voltage = signal.getValue(i)
        DOT = voltage >= _conf['tThreshold']
        DOE = voltage >= _conf['eThreshold']

        # DOx toggled
        if lastDOT != DOT:
            fileDOT.write("{0:012d}\n".format(t_ps - t_lastToggleDOT))
            t_lastToggleDOT = t_ps
        if lastDOE != DOE:
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

    fileDOT.close()
    fileDOE.close()


# the actual generation of the events
def generateChannelEvents(runTime, channel):
    global _conf

    # init the signal objects
    nSteps = int((runTime + 1e-6)/_conf['stepping'])
    discInput = Signal(nSteps)

    # generate the true hits
    generateTrues(discInput, runTime)

    # get the discriminator output
    getDiscriminatorOutput(discInput, channel)

    print "channel {0:2d}   generated Events: {1}".format(channel, discInput.getEvents())


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
    global _conf, _runTime, _nThreads, _maxThreads, _nChannels

    if 2 <= len(sys.argv) <= 3:
        if len(sys.argv) >= 3:
            if int(sys.argv[2]) > 0:
                _nThreads = int(sys.argv[3])

        _runTime = float(sys.argv[1])
    else:
        printUsage()
        sys.exit()

    _conf['directory'] = _conf['directory'].format(time=_runTime*1000)
    if not os.path.exists(_conf['directory']):
        os.makedirs(_conf['directory'])

    # check number of threads
    if _nThreads == 0 or _nThreads > _maxThreads:
        _nThreads = _maxThreads

    # 1) Init a Thread pool with the desired number of threads
    pool = ThreadPool(_nThreads)

    for ch in range(_nChannels):
        # 2) Add the task to the queue
        pool.add_task(generateChannelEvents, _runTime, ch)

    # 3) Wait for completion
    pool.wait_completion()


if __name__ == '__main__':
    main()
