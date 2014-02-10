#!/usr/bin/python
# -*- coding: utf-8

# Import stuff
import sys          # system functions (like exit)
import pylab as pl
import argparse

# we don't have to rewrite the import of the data
import generate_asic_data as asic
_conf = asic._conf


def plotPulse(data):
  t = data['t']
  V = data['V']
  pl.plot(t, V)
  pl.show()

def main():
  parser = argparse.ArgumentParser(description='Plot an input signal for the discriminator from front-end simulations.')
  parser.add_argument('charge', metavar='fC', type=float, default=5.0, help="input charge")
  args = parser.parse_args()

  inputFile = _conf['simulationInputFiles'][_conf['polarity']]
  pulse = asic.PulseShape(_conf['useSimulationInput'], inputFile)
  pulseData = pulse.selectPulse(args.charge)
  plotPulse(pulseData)

if __name__ == '__main__':
    main()
