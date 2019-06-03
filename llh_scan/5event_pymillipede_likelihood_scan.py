#!/cvmfs/icecube.opensciencegrid.org/py2-v1/icetray-start
#METAPROJECT /home/snowicki/icerec/calmdown/build/

## This script calculates a likelihood using clshim (a branch of clsim that generates
## event hypotheses) and the millipede implementation of the llh calculation. 
##
## new tactic - try using PyMillipede instead of MilliPlot, since 
## we just want to create hypotheses and calculate likelihoods
## PyMillipede is a c++ I3Module
 
from optparse import OptionParser
 
from icecube import icetray, dataclasses, dataio
from I3Tray import I3Tray
import numpy, pylab
 
import sys

usage = "%prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-i", "--inputfile", type="string", action="store", default="None", help="name of input file")
parser.add_option("-g", "--gcdfile", type="string", action="store", default="None", help="name of gcd file")
parser.add_option("-o", "--outputfile", type="string", action="store", default="None", help="name of output file")
parser.add_option("-x", "--x", type="int", action="store", default="0", help="step number being tested in x")
parser.add_option("-z", "--z", type="int", action="store", default="0", help="step number being tested in z")
parser.add_option("-a", "--xsteps", type="int", action="store", default="0", help="total number of steps in x")
parser.add_option("-b", "--zsteps", type="int", action="store", default="0", help="total number of steps in z")
parser.add_option("-c", "--xstepsize", type="int", action="store", default="0", help="size of steps in x")
parser.add_option("-d", "--zstepsize", type="int", action="store", default="0", help="size of steps in z")
parser.add_option("-r", "--randomnum", type="int", action="store", default="1337", help="seed for random number generator")
(options, args) = parser.parse_args()

infiles = [options.gcdfile, options.inputfile]

tray = I3Tray()
 
tray.AddModule('I3Reader', 'reader', filenamelist=infiles)
 
from icecube import photonics_service, millipede
 
def CLShim(tray, name):
    from icecube import icetray, dataclasses, clsim, phys_services
    import os
 
    from icecube.clsim.traysegments.common import configureOpenCLDevices, parseIceModel
 
    def getDetectorParameters(IceModel="spice_lea", DisableTilt=False, UnshadowedFraction=0.95,#1.,match the clsim setting
        UseHoleIceParameterization=True, DOMOversizeFactor=5., UnWeightedPhotons=False):
     
        DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
     
        # ice properties
        if isinstance(IceModel, str):
            mediumProperties = parseIceModel(os.path.expandvars('$I3_BUILD/clsim/resources/ice/%s' % IceModel), disableTilt=DisableTilt)
        else:
            # get ice model directly if not a string
            mediumProperties = IceModel
 
        # detector properties
        if UseHoleIceParameterization:
	    print "calculating DOM efficiency correction using hole ice"
            # the hole ice acceptance curve peaks at 0.75 instead of 1
            domEfficiencyCorrection = UnshadowedFraction*0.75*1.35 * 1.01 # DeepCore DOMs have a relative efficiency of 1.35 plus security margin of +1%
        else:
            domEfficiencyCorrection = UnshadowedFraction*1.35 * 1.01 # security margin of +1%
	print "from tray segment: domEfficiencyCorrection is  ", domEfficiencyCorrection
        domAcceptance = clsim.GetIceCubeDOMAcceptance(domRadius = DOMRadius*DOMOversizeFactor, efficiency=domEfficiencyCorrection)
	print "from tray segment: domAcceptance NumEntries is  ", domAcceptance.GetNumEntries()
        print "from tray segment: domAcceptance MinWlen is  ", domAcceptance.GetMinWlen()
        print "from tray segment: domAcceptance MaxWlen is  ", domAcceptance.GetMaxWlen()
        print "from tray segment: domAcceptance FirstWlen is  ", domAcceptance.GetFirstWavelength()
        print "from tray segment: domAcceptance WlenStepping is  ", domAcceptance.GetWavelengthStepping()
        print "from tray segment: domAcceptance 4th entry value is  ", domAcceptance.GetEntryValue(4)
        print "from tray segment: domAcceptance 4th entry wlen is  ", domAcceptance.GetEntryWavelength(4)
        domAngularSensitivity = clsim.GetIceCubeDOMAngularSensitivity(holeIce=UseHoleIceParameterization)
	print "from tray segment: domAngularSensitivity coefficients are  ", domAngularSensitivity.GetCoefficients()
	print "from tray segment: with hole ice  ", UseHoleIceParameterization
 
        if not UnWeightedPhotons:
            wavelengthGenerationBias = domAcceptance
        else:
            wavelengthGenerationBias = None
	print "from tray segment: wavelengthGenerationBias is  ", wavelengthGenerationBias
 
        wavelengthGenerators = clsim.I3CLSimRandomValuePtrSeries()
        wavelengthGenerators.append(clsim.makeCherenkovWavelengthGenerator(wavelengthGenerationBias, False, mediumProperties))
	print "from tray segment: wavelengthGenerators is  ", wavelengthGenerators#[0].NumberOfParameters() #no pybindings?
        return mediumProperties, wavelengthGenerationBias, wavelengthGenerators, domAcceptance, domAngularSensitivity, DOMOversizeFactor
 
    def makeStepConverter(UseGeant4=False):
        if UseGeant4:
        	clsim.AutoSetGeant4Environment()
        converter = clsim.I3CLSimLightSourceToStepConverterGeant4()
        ppcConverter = clsim.I3CLSimLightSourceToStepConverterPPC(photonsPerStep=200)
        parameterizationList = clsim.GetDefaultParameterizationList(ppcConverter, muonOnly=False)
        converter.SetLightSourceParameterizationSeries(parameterizationList)
         
        return converter
 
    def makeGeometry(positions, DOMOversizeFactor=5., spacing=30):
	#print positions
        import math
        string, dom = 0, 1
        lastx, lasty = float('nan'), float('nan')
     
        DOMRadius = 0.16510*icetray.I3Units.m # 13" diameter
     
        omGeos = dataclasses.I3OMGeoMap()
        moduleGeos = dataclasses.I3ModuleGeoMap()
        subdetectors = dataclasses.I3MapModuleKeyString()
        offsets = dict()
        for i, pos in enumerate(positions):
            if not math.hypot(pos.x-lastx, pos.y-lasty) < spacing:
                lastx, lasty = pos.x, pos.y
                string += 1
                dom = 1
            key = icetray.OMKey(string, dom)
            mkey = dataclasses.ModuleKey(string, dom)
            dom += 1
            omGeo = dataclasses.I3OMGeo()
            moduleGeo = dataclasses.I3ModuleGeo()
            omGeo.position = pos
            omGeo.omtype = omGeo.IceCube
            moduleGeo.pos = pos
            moduleGeo.radius = DOMRadius
            moduleGeo.module_type = moduleGeo.IceCube
            omGeos[key] = omGeo
            moduleGeos[mkey] = moduleGeo
            subdetectors[mkey] = "IceCube"
            offsets[i] = key
        frame = icetray.I3Frame()
        frame['I3OMGeoMap'] = omGeos
        frame['I3ModuleGeoMap'] = moduleGeos
        frame['Subdetectors'] = subdetectors
 
	simplegeo = clsim.I3CLSimSimpleGeometryFromI3Geometry(DOMRadius, DOMOversizeFactor, frame)

        return simplegeo
 
    def makeKernel(device, rng, detector_params_args=dict()):
        """
        Configure propagator and step generator, but do not initialize
        """
        mediumProperties, wavelengthGenerationBias, wavelengthGenerators, domAcceptance, domAngularSensitivity, DOMOversizeFactor = getDetectorParameters(**detector_params_args)
     
        propagator = clsim.I3CLSimStepToPhotonConverterOpenCL(rng)#, False)## turns "native math" off which is what is done for CPUs in the standard mode - not for GPUs, though
        propagator.SetDevice(device)
        propagator.SetMediumProperties(mediumProperties)
        propagator.SetWlenBias(wavelengthGenerationBias)
        propagator.SetWlenGenerators(wavelengthGenerators)
        propagator.SetDOMPancakeFactor(DOMOversizeFactor)
        propagator.SetMaxNumWorkitems(768*1400)

        stepGenerator = makeStepConverter()
        stepGenerator.SetMediumProperties(mediumProperties)
        stepGenerator.SetRandomService(rng)
        stepGenerator.SetWlenBias(domAcceptance)

        return stepGenerator, propagator, domAcceptance, domAngularSensitivity
     
    ##### Attempts to get CLShim to run on CPUs 
    #device = configureOpenCLDevices(UseGPUs=False, UseCPUs=True)[0]
    #device = configureOpenCLDevices(UseGPUs=False, UseCPUs=True, DoNotParallelize=True)[0]
    #device = configureOpenCLDevices(UseGPUs=False, UseCPUs=True, DoNotParallelize=False)[0]
    ##### This is if you want to use GPUs, which work just fine
    device = configureOpenCLDevices(UseGPUs=True, UseCPUs=False, DoNotParallelize=True)[0]
    print "from tray segment: is OpenCL device CPU?  ", device.IsCPU()
    print "from tray segment: is OpenCL device GPU?  ", device.IsGPU() 
    print "from tray segment: OpenCL device MaxComputeUnits  ", device.GetMaxComputeUnits()
    print "from tray segment: OpenCL device MaxWorkItemSize  ", device.GetMaxWorkItemSize()
    print "from tray segment: OpenCL device MaxWorkGroupSize  ", device.GetMaxWorkGroupSize()
    print "from tray segment: OpenCL device MaxClockFrequencyMhz  ", device.GetMaxClockFrequencyMhz()
    print "from tray segment: OpenCL device GlobalMemSize  ", device.GetGlobalMemSize()
    print "from tray segment: OpenCL device MaxConstantBufferSize  ", device.GetMaxConstantBufferSize()
    print "from tray segment: OpenCL device LocalMemSize  ", device.GetLocalMemSize()
    print "from tray segment: OpenCL device has DedicatedLocalMem?  ", device.HasDedicatedLocalMem()
    print "from tray segment: OpenCL device has ErrorCorrectionSupport?  ", device.HasErrorCorrectionSupport()
    print "from tray segment: OpenCL device is available?  ", device.IsAvailable()
    print "from tray segment: OpenCL device Vendor?  ", device.GetVendor()
    print "from tray segment: OpenCL device DriverVersion?  ", device.GetDriverVersion()
    print "from tray segment: OpenCL device DeviceVersion?  ", device.GetDeviceVersion()
    print "from tray segment: OpenCL device Extensions?  ", device.GetExtensions()

    #rng = phys_services.I3GSLRandomService(1337)
    rng = phys_services.I3GSLRandomService(options.randomnum)
    stepGenerator, propagator, domAcceptance, domAngularSensitivity = makeKernel(device, rng)
 
    tray.AddService('I3CLShimFactory', name, AngularAcceptance=domAngularSensitivity,
        WavelengthAcceptance=domAcceptance, StepConverter=stepGenerator, PhotonConverter=propagator,
        GeometryFactory=makeGeometry, OversampleFactor=500, CacheDepth=5)
 
def ProcessEvent(frame):
	if frame['I3EventHeader'].event_id not in [2209, 9993, 10404, 24550, 30171]:
        #if frame['I3EventHeader'].event_id != one_true_event_id:
                print("Skipping event %d" % (frame['I3EventHeader'].event_id) )
                return False
        else:
                print("Running on event %d" % (frame['I3EventHeader'].event_id) )
                return True

tray.AddModule(ProcessEvent, 'shouldweprocessthisevent')
 
pxs = "CLShim"
cascade_base = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.%s.fits"
muon_base = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/emu_%s.fits"
cascade_tables = photonics_service.I3PhotoSplineService(cascade_base % "abs", cascade_base % "prob", 0)
muon_tables = photonics_service.I3PhotoSplineService(muon_base % "abs", muon_base % "prob", 0)

tray.AddSegment(CLShim, pxs)
#icetray.logging.I3Logger.global_logger.set_level_for_unit('I3PhotoSplineService', icetray.logging.I3LogLevel.LOG_TRACE);

## pulse series and time range JP uses for pegleg (MSU)
## time window is used for noise generation
Pulses='SRTTWOfflinePulsesDC'
ReadoutWindow='WaveformRange'

ExcludedDOMs = ['CalibrationErrata','BadDomsList']
 
## have to read in the X, Z as input arguments? then have to look at the same space for all events
## or... instead of reading in x & z positions, read in the x & z step numbers
def make_seed(frame):
	truth = frame['trueCascade']
	trueX = truth.pos.x
	trueY = truth.pos.y
	trueZ = truth.pos.z
	#print trueX, trueZ
	xsteps = options.xsteps
	#ysteps = options.zsteps
	zsteps = options.zsteps
	#print xsteps, zsteps
	xstepsize = options.xstepsize
	#ystepsize = options.zstepsize
	zstepsize = options.zstepsize
	#print xstepsize, zstepsize
	startX = trueX - (xsteps/2.)*xstepsize
	#startY = trueY - (ysteps/2.)*ystepsize
	startZ = trueZ - (zsteps/2.)*zstepsize
	#print startX, startZ
	testX = startX + options.x*xstepsize
	#testY = startY + options.z*ystepsize
	testZ = startZ + options.z*zstepsize
	#print testX, testZ
	newparticle = dataclasses.I3Particle(truth)
	newparticle.pos.x = testX
	#newparticle.pos.y = testY
	newparticle.pos.z = testZ
	#print newparticle
	frame['trueCascade_seed'] = newparticle

tray.AddModule(make_seed, "make_seed", streams=[icetray.I3Frame.Physics])
tray.AddModule('PyMillipede', 'trueCascade_seed', Hypothesis=lambda fr: dataclasses.I3VectorI3Particle([fr['trueCascade_seed']]), Output='pyMillipede_likelihood', ExcludedDOMs=ExcludedDOMs, CascadePhotonicsService=pxs, MuonPhotonicsService=pxs, Pulses=Pulses, ReadoutWindow=ReadoutWindow)
tray.AddModule('PyMillipede', 'trueCascade_seed2', Hypothesis=lambda fr: dataclasses.I3VectorI3Particle([fr['trueCascade_seed']]), Output='pyMillipede_likelihood_tables', ExcludedDOMs=ExcludedDOMs, CascadePhotonicsService=cascade_tables, MuonPhotonicsService=muon_tables, Pulses=Pulses, ReadoutWindow=ReadoutWindow)

tray.AddModule('I3Writer', 'writer', filename=options.outputfile, DropOrphanStreams=[icetray.I3Frame.DAQ], streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics])#streams=["DAQ", "Physics"]
 
tray.AddModule('TrashCan', 'YesWeCan')
tray.Execute()
tray.Finish()
