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

Infile_List = glob.glob(params.Infile)

event_counter = 0
def event_finder(frame):
	global event_counter
	if event_counter > 0:
		return False
	if not(frame.Has("IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC")):
		return False
	reco_neutrino = frame["IC86_Dunkman_L6_PegLeg_MultiNest8D_NumuCC"]
	true_neutrino = frame["trueNeutrino"]
	if true_neutrino.energy < 10 or true_neutrino.energy > 11:
		return False
	if reco_neutrino.energy < 9 or reco_neutrino.energy > 12:
		return False
	if abs(reco_neutrino.pos.x - true_neutrino.pos.x) > 10:
		return False
	if abs(reco_neutrino.pos.y - true_neutrino.pos.y) > 10:
		return False
	if abs(reco_neutrino.pos.z - true_neutrino.pos.z) > 10:
		return False
	if abs(reco_neutrino.dir.zenith - true_neutrino.dir.zenith) > 20*pi/180:
		return False
	if abs(reco_neutrino.dir.azimuth - true_neutrino.dir.azimuth) > 50*pi/180:
		return False
	event_counter += 1

tray = I3Tray()

tray.AddModule("I3Reader", "reader", filenamelist=Infile_List)
tray.AddModule(event_finder, "event_finder")
tray.AddModule("I3Writer", "writer", filename=params.Outfile,DropOrphanStreams = [icetray.I3Frame.DAQ])

tray.AddModule( 'TrashCan' , 'Done' )

if (params.NEvents==-1):
    tray.Execute()
else:
    tray.Execute(params.NEvents)

tray.Finish()
