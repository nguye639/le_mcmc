from I3Tray import *
from icecube import icetray, dataclasses, dataio
import glob
from ROOT import TH1D, TCanvas, TFile
from math import *

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

cascade_base = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.%s.fits"
muon_base = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/emu_%s.fits"
cascade_tables = photonics_service.I3PhotoSplineService(cascade_base % "abs", cascade_base % "prob", 0)
muon_tables = photonics_service.I3PhotoSplineService(muon_base % "abs", muon_base % "prob", 0)

def make_seed(frame, garbage = 0):
    truth = frame['trueNeutrino']
	
    newparticle = dataclasses.I3Particle(truth)
    newparticle.pos.x += garbage
    frame['trueNeutrino_seed_'+str(garbage)] = newparticle

Infile_List = glob.glob(params.Infile)

def print_llh(frame,garbage = 0):
	llh_result = frame["pyMillipede_"+str(garbage)+"FitParams"]
	print garbage, llh_result.logl


for i in xrange(-10,10):	
	tray = I3Tray()
	
	tray.AddModule("I3Reader", "reader", filenamelist=[params.GCDfile]+Infile_List)
	tray.AddModule(make_seed, "make_seed_"+str(i), garbage=i,streams=[icetray.I3Frame.Physics])
	tray.AddModule('PyMillipede', 'trueNeutrino_seed_'+str(i), Hypothesis=lambda fr: dataclasses.I3VectorI3Particle([fr['trueNeutrino_seed_'+str(i)]]), Output='pyMillipede_'+str(i), ExcludedDOMs=ExcludedDOMs, CascadePhotonicsService=cascade_tables, MuonPhotonicsService=muon_tables, Pulses=Pulses, ReadoutWindow=ReadoutWindow)	
	tray.AddModule(print_llh,"print_llh_"+str(i),garbage=i)
	tray.AddModule("I3Writer", "writer", filename=params.Outfile,DropOrphanStreams = [icetray.I3Frame.DAQ], streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics])
	
	tray.AddModule( 'TrashCan' , 'Done' )
	
	if (params.NEvents==-1):
		tray.Execute()
	else:
		tray.Execute(params.NEvents)
	
	tray.Finish()
