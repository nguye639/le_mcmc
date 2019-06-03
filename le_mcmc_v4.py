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
#parser.add_option("-w", "--output")

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

#stepping for the MCMC
def variance(frame,pdg=11):

    rand = np.random.uniform(-1,1)
    if event_counter >=0 and event_counter < 10000:
        vertex_step = 10*np.random.uniform(-1,1)
        time_step = 50*np.random.uniform(-1,1)
        zenith_step = 10*np.random.uniform(-1,1)
        azimuth_step = 20*np.random.uniform(-1,1)
        energy_step = 5*np.random.uniform(-1,1)
        length_step = 18*np.random.uniform(-1,1)
    elif event_counter >= 10000 and event_counter < 100000:
        vertex_step = 5*np.random.uniform(-1,1)
        time_step = 25*np.random.uniform(-1,1)
        zenith_step = 5*np.random.uniform(-1,1)
        azimuth_step = 10*np.random.uniform(-1,1)
        energy_step = 2*np.random.uniform(-1,1)
        length_step = 4.5*np.random.uniform(-1,1)
    else:
        vertex_step = np.random.uniform(-1,1)
        time_step = np.random.uniform(-1,1)
        zenith_step = np.random.uniform(-1,1)
        azimuth_step = np.random.uniform(-1,1)
        energy_step = np.random.uniform(-1,1)
        length_step = np.random.uniform(-1,1)

    global parameters
    global temp_parameters
    i = 0
    for i in range(len(temp_parameters)):
        n = np.random.random()
       
        if n >= .5:
            
            if i ==0:
	        if temp_parameters[i] < 196 and temp_parameters[0] > -104:
                    temp_parameters[i] += vertex_step
                                        
            elif i ==1:
		if temp_parameters[i] < 116 and temp_parameters[1] > -184:
                    temp_parameters[i] += vertex_step
                        
            elif i ==2:
		if temp_parameters[i] < -200  and temp_parameters[2] > -500:
                    temp_parameters[i] += vertex_step
                          
            elif i == 3:
		temp_parameters[i] += azimuth_step*np.pi/180
		if temp_parameters[i] > 2*np.pi:
		    temp_parameters[i] -= 2*np.pi
		elif temp_parameters[i] < 0:
		    temp_parameters[i] += 2*np.pi
                    
            elif i == 4:
		temp_parameters[i] += zenith_step*np.pi/180
                if temp_parameters[i] > np.pi:
		    temp_parameters[i] = np.pi
		elif temp_parameters[i] < 0:
		    temp_parameters[i] = 0

	    elif i == 5:
		if temp_parameters[i]+time_step > 0:
                    temp_parameters[i] += time_step

            elif i == 6:
		if temp_parameters[i]+energy_step > 0:
                    temp_parameters[i] += energy_step

            elif i == 7:
		if temp_parameters[i]+length_step > 0:
                    temp_parameters[i] += length_step

            elif i == 8:
                temp_parameters[i] += azimuth_step*np.pi/180
                if temp_parameters[i] > 2*np.pi:
                    temp_parameters[i] -= 2*np.pi
                elif temp_parameters[i] < 0:
                    temp_parameters[i] += 2*np.pi
                     
            elif i == 9:
                temp_parameters[i] += zenith_step*np.pi/180
                if temp_parameters[i] > np.pi: 
                    temp_parameters[i] = np.pi
                elif temp_parameters[i] < 0:
                    temp_parameters[i] = 0

#choosing to update parameters                                     
def compare(frame):		
	
	global parameters
	global temp_parameters
	global llh_value
	n = np.random.random_sample()
	p = llh_value[-2] - llh_value[-1]

	if p > 34:
		x = 1
	elif p < -34:
		x = 0
	else:
		x = 10**(p)

	if x > n:
		parameters = list(temp_parameters)
		history_parameters.append(list(parameters))
	else:
		llh_value.pop()
		temp_parameters = list(parameters)			

event_counter = 0
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

def update_parameters(frame,pdg=11):
	global parameters
	global temp_parameters
	pegleg_7D = frame["IC86_Dunkman_L6_PegLeg_MultiNest7D_NumuCC"]
	parameters[0] = pegleg_7D.pos.x	
	parameters[1] = pegleg_7D.pos.y
	parameters[2] = pegleg_7D.pos.z
	parameters[3] = pegleg_7D.dir.azimuth
	parameters[4] = pegleg_7D.dir.zenith
	parameters[8] = pegleg_7D.dir.azimuth
	parameters[9] = pegleg_7D.dir.zenith
	parameters[5] = pegleg_7D.time
	parameters[6] = pegleg_7D.energy
	temp_parameters = list(parameters)


def get_truth(frame):
    #ADDED BY JESSIE
    #Returns: [x_vertex, y_vertex, z_vertex, azimuth, zenith, time, hadron energy, track length, hadron azimuth, hadron zenith]
    global truth
    truth = []

    global history_parameters
    history_parameters = []

    global parameters
    global temp_parameters

    true_nu = frame["trueNeutrino"]
    true_mu = frame["trueMuon"]
    true_casc = frame["trueCascade"]
    
    nu_energy = true_nu.energy
    nu_zen    = true_nu.dir.zenith
    nu_azi    = true_nu.dir.azimuth
    nu_x      = true_nu.pos.x
    nu_y      = true_nu.pos.y
    nu_z      = true_nu.pos.z
    nu_time   = true_nu.time
    mu_length = true_mu.length

    had_energy = true_casc.energy
    had_zen    = true_casc.dir.zenith
    had_azi    = true_casc.dir.azimuth

    truth.append([ nu_x, nu_y, nu_z, nu_azi, nu_zen, nu_time, had_energy, mu_length, had_azi, had_zen])
    print(truth)
    truth = truth[0]
    #parameters = [x_vertex, y_vertex, z_vertex, azimuth, zenith, time, hadron energy, track length, hadron azimuth, hadron zenith]

    parameters = list(truth)

    temp_parameters = list(truth)

    history_parameters.append(list(parameters))
#TO USE get_truth FUNCTION:
# truth = []
# AddModule(get_truth, "Get Truth")
# To reference things in truth: truth[0][0] = nu_x, truth[0][1] = nu_y, etc.
# NOTE: if you run this module without resetting truth = [] you can get multiple sets of truth values.
# So first time it is run truth[0][i] = truth values, second time it is run truth[1][i] is new set, etc

Infile_List = glob.glob(params.Infile)

#calculates llh value
llh_value = []
def get_llh(frame,garbage=0):
	global event_counter
	global llh_value 	
	llh_result = frame["pyMillipede_"+str(event_counter)+"SeedParams"]
	llh_value.append(llh_result.logl)
	event_counter += 1

def get_DC_cascade_hypothesis(frame):
	hypothesis = frame["DC_cascade_hypothesis_"+str(event_counter)]
	return hypothesis

#finds the fit for a single parameter
def find_fit(x,y):
    for i in range(0,len(y),1):
        if y[i] == max(y):
            j = i

    fit.append(x[j+1])

#finding best fit
def best_fit():
	global fit 
	fit = []

        x_data = []
        y_data = []
        z_data = []
        az_lepton_data = []
        ze_lepton_data = []
        time_data = []
        energy_data = []
        track_data = []
        az_hadron_data = []
        ze_hadron_data = []
        for i in range(len(history_parameters)):
                x_data.append(history_parameters[i][0])
                y_data.append(history_parameters[i][1])
                z_data.append(history_parameters[i][2])
                az_lepton_data.append(history_parameters[i][3])
                ze_lepton_data.append(history_parameters[i][4])
                time_data.append(history_parameters[i][5])
                energy_data.append(history_parameters[i][6])
                track_data.append(history_parameters[i][7])
                az_hadron_data.append(history_parameters[i][8])
                ze_hadron_data.append(history_parameters[i][9])
	x=0
        y=0

        binsize = 100
        y, x = np.histogram(x_data, bins = binsize)
        find_fit(x,y)
        y, x = np.histogram(y_data, bins = binsize)
        find_fit(x,y)
        y, x = np.histogram(z_data, bins = binsize)
        find_fit(x,y)
        y, x = np.histogram(az_lepton_data, bins = binsize)
        find_fit(x,y)
        y, x = np.histogram(ze_lepton_data, bins = binsize)
        find_fit(x,y)
        y, x = np.histogram(time_data, bins = binsize)
        find_fit(x,y)
        y, x = np.histogram(energy_data, bins = binsize)
        find_fit(x,y)
        y, x = np.histogram(track_data, bins = binsize)
        find_fit(x,y)
        y, x = np.histogram(az_hadron_data, bins = binsize)
        find_fit(x,y)
        y, x = np.histogram(ze_hadron_data, bins = binsize)
        find_fit(x,y)



#Running all of the functions for j*i iterations
for j in xrange(0,2160):
	tray = I3Tray()
	tray.AddModule("I3Reader", "reader", filenamelist=[params.GCDfile]+Infile_List)
	if j == 0:
    		tray.AddModule(get_truth, "Get Truth")

	for i in xrange(0,1000):
		if not( i == 0 and j == 0):
			tray.AddModule(variance, "variance"+str(i+1000*j), pdg = 11)
		tray.AddModule(make_DC_muon_hadron_hypothesis,"make_DC_muon_hadron_hypothesis_"+str(i+1000*j))
		tray.AddModule('PyMillipede', 'test_DC_cascade_hypothesis_'+str(i+1000*j), Hypothesis=get_DC_cascade_hypothesis, Output='pyMillipede_'+str(i+1000*j), ExcludedDOMs=ExcludedDOMs, CascadePhotonicsService=cascade_tables, MuonPhotonicsService=muon_tables, Pulses=Pulses, ReadoutWindow=ReadoutWindow)	
		tray.AddModule(get_llh,"get_llh_"+str(i+1000*j))
		if not( i == 0 and j == 0):
			tray.AddModule(compare,"compare"+str(i+1000*j))
	tray.AddModule( 'TrashCan' , 'Done' )
	if (params.NEvents==-1):
		tray.Execute()
	else:
		tray.Execute(params.NEvents)
	tray.Finish()

#Finding the best fit and error after the chain runs
best_fit()
fit = np.array(fit)
truth = np.array(truth)

print("BEST FIT: ")
print(fit)


#pickling out the list of parameters for later analysis
import pickle

pickle_out = open(params.Outfile,"wb")
pickle.dump(history_parameters,pickle_out)
pickle_out.close()
