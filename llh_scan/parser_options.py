
class RunParameters:
  def __init__(self):
    self.Infile = ''
    self.Outfile = ''
    self.GCDfile = ''
    self.NEvents = 0
    self.NSkip = 0
    self.tolerance = 0
    self.penalty = 0

def parseOptions(parser, params):
	parser.add_option(
                  "-i", "--inputfile",
                  type      = "string",
                  action    = "store",
                  default   = "None",
                  dest      = "INPUT",
                  metavar   = "<input file>",
                  help      = "Name of the input file",
                  )
	parser.add_option(
                  "-g", "--gcdfile",
                  type      = "string",
                  action    = "store",
                  default   = "None",
                  dest      = "GCD",
                  metavar   = "<geo file>",
                  help      = "Name of GCD file",
                  )
	parser.add_option(
                  "-n", "--numevents",
                  type      = "int",
                  action    = "store",
                  default   = -1,
                  dest      = "NEVENTS",
                  help      = "Number of physics events to process (default: all)",
                  )
	parser.add_option(
                  "-s", "--skip",
                  type      = "int",
                  action    = "store",
                  default   = 0,
                  dest      = "NSKIP",
                  help      = "Number of events we want to skip (default: 0)",
                  )
	parser.add_option(
                  "-o", "--outputfile",
                  type      = "string",
                  action    = "store",
                  default   = "tmp",
                  dest      = "OUTPUT",
                  metavar   = "<output file(s) name>",
                  help      = "Name of the output file(s), i.e. .root and .i3.gz names",
                  )
        parser.add_option(
		'-t', '--tolerance',
		 type = float,
		 default = .3
		 )
        parser.add_option(
		'-p', '--penalty',
		 type = float,
		 default = .1
		 )

	(options,args) = parser.parse_args()
	params.Infile = options.INPUT
	params.Outfile = options.OUTPUT
	params.GCDfile = options.GCD
	params.NEvents = options.NEVENTS
	params.NSkip = options.NSKIP
	params.tolerance = options.tolerance
	params.penalty = options.penalty

