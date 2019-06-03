from I3Tray import *
from icecube import icetray, dataclasses, dataio
import glob
from ROOT import TH1D, TCanvas, TFile
from math import *
import numpy as np
sys.path.insert(0,"/mnt/home/neergarr/icecube/")
from optparse import OptionParser
from parser_options import *
params = RunParameters()

usage = "%usage: %prog [options]"
parser = OptionParser(usage)
parseOptions(parser, params)

Pulses='SRTTWOfflinePulsesDC'
ReadoutWindow='WaveformRange'

ExcludedDOMs = ['CalibrationErrata','BadDomsList']
from icecube import photonics_service, millipede

cascade_base = "/mnt/home/neergarr/icecube/pegleg_testing/splines/ems_mie_z20_a10.%s.fits"
muon_base = "/mnt/home/neergarr/icecube/pegleg_testing/splines/emu_%s.fits"
cascade_tables = photonics_service.I3PhotoSplineService(cascade_base % "abs", cascade_base % "prob", 0)
muon_tables = photonics_service.I3PhotoSplineService(muon_base % "abs", muon_base % "prob", 0)

muon_energy_loss = 4.5 # meters/GeV (divide track length by muon energy loss to get muon energy)
hadron_pdg = -2000001006
global llh_value
llh_value = []
global temp_parameters


#              0. 1. 2.        3.              4.   5.    6.        7.              8.               9.   
#parameters = [x,y,z,lepton azimuth,lepton zenith,time,energy, track lengeth, hadron azimuth, hadron zenith]

parameters = [59.44779725056501, -50.40980604647394, -315.12215116285233, 0.760880564153089, 1.717751242710285, 9798.736378713807, 0.5142537959869027, 215.66295627936717, 3.7698636161288634, 0.12566370614359174]


global temp_parameters
temp_parameters = list(parameters)

global history_parameters
history_parameters = []

def init(frame, j):
	if j == 0 or j == 1 or j == 2:
		temp_parameters[j] -= 100

	elif j == 3 or j == 8:
		temp_parameters[j] -= 2*np.pi

	elif j == 4 or j== 9:
		temp_parameters[j] -= .5*np.pi
	elif j==5:
		temp_parameters[j] -= 100
	else:
		temp_parameters[j] -= 50	


def mapping(frame,i,j):
	if j == 0 or j == 1 or j == 2:
                temp_parameters[j] += 100/(.5*steps)

        elif j == 3 or j == 8:
                temp_parameters[j] += (2*np.pi)/(.5*steps)

        elif j == 4 or j== 9:
                temp_parameters[j] += (.5*np.pi)/(.5*steps)
        elif j==5:
                temp_parameters[j] += 100/(.5*steps)
        else:
                temp_parameters[j] += 50/(.5*steps)


	history_parameters.append(temp_parameters[:])

def reset(frame):
	global temp_parameters
	temp_parameters = list(parameters)


event_counter  = 0
#calculates llh value
def get_llh(frame,garbage=0):
        global event_counter
        global llh_value
        llh_result = frame["pyMillipede_"+str(event_counter)+"SeedParams"]
        llh_value.append(llh_result.logl)
        event_counter += 1

#different hypotheses to choose from
def make_DC_cascade_hypothesis(frame, pdg=11):
	new_particle = dataclasses.I3Particle()
	new_particle.pos.x = temp_parameters[0]
	new_particle.pos.y = temp_parameters[1]
	new_particle.pos.z = temp_parameters[2]
	new_particle.time = temp_parameters[5]
	new_direction = dataclasses.I3Direction(temp_parameters[4], temp_parameters[3])
	new_particle.dir = new_direction
	new_particle.energy = temp_parameters[6]
	new_particle.pdg_encoding = pdg
	hypothesis = dataclasses.I3VectorI3Particle()
	hypothesis.append(new_particle)
	frame["DC_cascade_hypothesis_"+str(event_counter)] = hypothesis


def make_DC_muon_hadron_hypothesis(frame, pdg=13):
	new_muon = dataclasses.I3Particle()
	new_muon.pos.x = temp_parameters[0]
	new_muon.pos.y = temp_parameters[1]
	new_muon.pos.z = temp_parameters[2]
	new_muon.time = temp_parameters[5]
	new_direction = dataclasses.I3Direction(temp_parameters[4], temp_parameters[3])
	new_muon.dir = new_direction
	new_muon.length = temp_parameters[7]
	new_muon.energy = temp_parameters[7]/muon_energy_loss
	new_muon.pdg_encoding = pdg
	hypothesis = dataclasses.I3VectorI3Particle()
	hypothesis.append(new_muon)
	new_hadron = dataclasses.I3Particle()
	new_hadron.pos.x = temp_parameters[0]
	new_hadron.pos.y = temp_parameters[1]
	new_hadron.pos.z = temp_parameters[2]
	new_hadron.time = temp_parameters[5]
	#change temp_parameters to 3, to 4 here to say hadron and muon angles are the same
	new_direction = dataclasses.I3Direction(temp_parameters[9], temp_parameters[8])
	new_hadron.dir = new_direction
	new_hadron.energy = temp_parameters[6]
	new_hadron.pdg_encoding = hadron_pdg
	hypothesis.append(new_hadron)
	frame["DC_cascade_hypothesis_"+str(event_counter)] = hypothesis

Infile_List = glob.glob(params.Infile)

def get_DC_cascade_hypothesis(frame):
	hypothesis = frame["DC_cascade_hypothesis_"+str(event_counter)]
	return hypothesis


#Running all of the functions for j*i iterations
global steps
steps = 100

for j in xrange(0,10):
	tray = I3Tray()
	tray.AddModule("I3Reader", "reader", filenamelist=[params.GCDfile]+Infile_List)
	tray.AddModule(init, "init" + str(j), j=j)
	for i in xrange(0,steps):
		tray.AddModule(mapping,"mapping"+str(i+steps*j),i=i,j=j)
		tray.AddModule(make_DC_muon_hadron_hypothesis,"make_DC_muon_hadron_hypothesis_"+str(i+steps*j))
		tray.AddModule('PyMillipede', 'test_DC_cascade_hypothesis_'+str(i+steps*j), Hypothesis=get_DC_cascade_hypothesis, Output='pyMillipede_'+str(i+steps*j), ExcludedDOMs=ExcludedDOMs, CascadePhotonicsService=cascade_tables, MuonPhotonicsService=muon_tables, Pulses=Pulses, ReadoutWindow=ReadoutWindow)	
		tray.AddModule(get_llh,"get_llh_"+str(i+steps*j))
	tray.AddModule(reset, "reset" + str(j))
	tray.AddModule( 'TrashCan' , 'Done' )
	if (params.NEvents==-1):
		tray.Execute()
	else:
		tray.Execute(params.NEvents)
	tray.Finish()

llh_list = []
llh_list.append(history_parameters)

for i in range(1,11):
	llh_list.append((llh_value[steps*(i-1):steps*i]))

import pickle
f = open(params.Outfile, 'wb')

pickle.dump(llh_list,f)
f.close()

