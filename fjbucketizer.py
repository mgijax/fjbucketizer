#!/usr/local/bin/python 
#
# fjbucketizer.py
#
#--------------------------------------------

import sys
import time
import string
from tempfile import mkstemp
import os
import os.path
from TableTools import TA,TB,TD,TF,TI,TJ,TP,TS,TU,TX,FJ
from optparse import OptionParser

USAGE="%prog [options] --f1 file1.gff --f2 file2.gff\n(For help, use -h option.)"

def now():
	return time.asctime(time.localtime(time.time()))

class GUPipeline :
    #----------------------------------------------------
    def __init__(self, argv):
	# 1. initialize the option parser
	#
	self.argv = argv
	self.parser = OptionParser(USAGE)
	self.tempfiles = []
	self.types1 = []
	self.types2 = []

	self.parser.add_option("-1", "--f1", "--file1",
		dest="file1",
		metavar="FILE1",
		default=None,
		help="The first GFF file.")

	self.parser.add_option("-2", "--f2", "--file2",
		dest="file2",
		metavar="FILE2",
		default=None,
		help="The second GFF file.")

	self.parser.add_option("-k", "--minOverlap", 
		dest="k",
		metavar=" AMT",
		default="1",
		help="The minimum required overlap. (Default: 1)")

	self.parser.add_option("--t1", 
	    dest="types1",
	    action="append",
	    metavar="GFFTYPE",
	    default=[],
	    help="A GFF type to select from file 1. Repeatable. (Default = all types)" )

	self.parser.add_option("--t2", 
	    dest="types2",
	    action="append",
	    metavar="GFFTYPE",
	    default=[],
	    help="A GFF type to select from file 2. Repeatable.(Default = all types)" )

	self.parser.add_option("--nt1", 
	    dest="notTypes1",
	    action="append",
	    metavar="GFFTYPE",
	    default=[],
	    help="A GFF type to FILTER OUT of file 1. Repeatable. (Default = filters no types)" )

	self.parser.add_option("--nt2", 
	    dest="notTypes2",
	    action="append",
	    metavar="GFFTYPE",
	    default=[],
	    help="A GFF type to filter out of file 2. Repeatable.(Default = filters no types)" )

	self.parser.add_option("-i", "--ignoreStrand", 
	    dest="ignoreStrand",
	    action="store_true",
	    default = False,
	    help="Ignore strand when determining overlaps (Default = strands must match)")

	self.parser.add_option("-n", "--noAggregate", 
	    dest="noAggregate",
	    action="store_true",
	    default = False,
	    help="Input features are independent. (Default: features connected by common ID in col 9)")

	self.parser.add_option("-s", "--noSelfHits", 
	    dest="noSelfHits",
	    action="store_true",
	    default = False,
	    help="Ignore self-hit in overlap detection. (default: reports self-hits)")

	self.parser.add_option("-p", "--preserve", 
	    dest="preserve",
	    action="store_true",
	    default = False,
	    help="Preserves temporary files. (default: temp files deleted)" )

	self.parser.add_option("-d", "--od", 
	    dest="outDir",
	    metavar="DIRECTORY",
	    default=".",
	    help="Path of output directory. (Default = current directory)")

	self.parser.add_option("-t", 
	    dest="template",
	    metavar="TEMPLATE",
	    default="bucket_%s.txt",
	    help="Template string for naming output bucket files. " + \
	    "Must contain the substring '%s', which is replaced by the bucket class. " + \
	    "The classes are: '1-0', '0-1', '1-1', '1-n', 'n-1', 'n-m'. " + \
	    "(Default: bucket_%s.txt)")

	self.parser.add_option("-l", "--lf", "--logFile",
	    dest="logFile",
	    metavar="FILE",
	    default=None,
	    help="Log file. (Default = writes to stderr)")

	# 2. Parse the command line
	#
	(self.options,xxx) = self.parser.parse_args(argv)

	# 3. Validate the args
	#
	if self.options.logFile is not None:
	    sys.stderr = open(self.options.logFile, 'a')

	if self.options.file1 is None \
	or self.options.file2 is None:
	    self.parser.error("Must specify both --f1 and --f2.")

	if not os.path.isdir(self.options.outDir):
	    self.parser.error("Output directory " + self.options.outDir + \
	    	" does not exist or is not a directory.")

	self.types1 = str( self.options.types1 ).lower()
	self.types2 = str( self.options.types2 ).lower()
	self.notTypes1 = str( self.options.notTypes1 ).lower()
	self.notTypes2 = str( self.options.notTypes2 ).lower()

	self.minOverlap = self.options.k

	if not '%s' in self.options.template:
	    self.options.template += '%s'

	self.options.template = \
	        os.path.join( self.options.outDir, self.options.template)
	
    #----------------------------------------------------
    def mkTmp(self):
	tf = mkstemp(dir=self.options.outDir)
	os.close(tf[0])
	if not self.options.preserve:
	    self.tempfiles.append(tf[1])
	return tf[1]

    #----------------------------------------------------
    def cleanupTempFiles(self):
	for tf in self.tempfiles:
	    os.remove(tf)

    #----------------------------------------------------
    def debug(self,s,ts=False):
	if ts:
	    sys.stderr.write(now()+': ')
	sys.stderr.write(s)

    #----------------------------------------------------
    def execStep(self, tool, args):
	self.debug( tool.__name__ + " " + string.join(args, " ") + "\n", True)
    	t = tool(args)
	t.go()
	t.closeFiles()
	return t.nOutputRows

    #----------------------------------------------------
    def go_noAggregate(self):

	### Select rows from first file.

	args = [
	    "--file1="+self.options.file1,
	    ]

	if len(self.options.types1) > 0:
	    args.append( "?string.lower(IN[3]) in %s" % self.types1 )

	if len(self.options.notTypes1) > 0:
	    args.append( "?string.lower(IN[3]) not in %s" % self.notTypes1 )

	if len(args) > 1:
	    f1 = self.mkTmp()
	    args.append( "--out-file=" + f1 )
	    self.execStep(TF, args)
	else:
	    f1 = self.options.file1

	### Select rows from second file.

	args = [
	    "--file1="+self.options.file2,
	    ]

	if len(self.options.types2) > 0:
	    args.append( "?string.lower(IN[3]) in %s" % self.types2 )

	if len(self.options.notTypes2) > 0:
	    args.append( "?string.lower(IN[3]) not in %s" % self.notTypes2 )

	if len(args) > 1:
	    f2 = self.mkTmp()
	    args.append( "--out-file=" + f2 )
	    self.execStep(TF, args)
	else:
	    f2 = self.options.file2

	### find overlapping features. 

	overlaps=self.mkTmp()
	args = [
	    "-1", f1,
	    "-2", f2,
	    "-s", "both",
	    "-k", self.minOverlap,
	    "-o", overlaps,
	    ]
	if self.options.ignoreStrand:
	    args = args + ["--columns1", "1,4,5"]
	novl = self.execStep(FJ, args)

	if self.options.noSelfHits:
	    xxx = overlaps
	    overlaps = self.mkTmp()
	    args = [ 
	    	"?IN[10]!=IN[19]",
		"-1", xxx,
		"-o", overlaps,
		]
	    novl = self.execStep(TF, args)
	    #self.debug("FJ out: " + xxx)
	    #self.debug("TF out: " + overlaps)
	    #os.system("diff %s %s" %(xxx,overlaps))


	if novl == 0:
	    self.debug("No overlapping features detected.\n")
	    #self.debug("No overlapping features detected. Pipeline terminated.\n")
	    #sys.exit(-1)

	### bucketize the pairs.

	bucketized=self.mkTmp()
	self.execStep(TB, [
	    "--file1=" + overlaps,
	    "--k1=10",
	    "--k2=19",
	    "-t"+bucketized,
	    "IN[1:]",
	    "int(IN[11]!=IN[20])",	## compute column: 0==same strands 1==diff
	    "int(string.lower(IN[7])!='gene')",	## compute column: 0==all genes 1==nongene
	    ])

	sorted = self.mkTmp()
	self.execStep(TS, [
	    "--file1=" + bucketized,
	    "-k 3",
	    "-k 1",
	    "--out-file=" + sorted,
	    ])

	self.execStep(TP, [
	    "--file1=" + sorted,
	    "-o" + self.options.template,
	    "-p 3",
	    ])

	### Bucketization did not generate 1-0 and 0-1 buckets
	### (because we only fed it overlapping pairs).
	### Generate these buckets by diff'ing the inputs
	### against the fjoin output.

	self.execStep(TD, [
	    "--file1=" + f1,
	    "--k1=9",
	    "--file2=" + bucketized,
	    "--k2=13",
	    "--out-file=" + (self.options.template%"1-0") ])

	self.execStep(TD, [
	    "--file1=" + f2,
	    "--k1=9",
	    "--file2=" + bucketized,
	    "--k2=22",
	    "--out-file=" + (self.options.template%"0-1") ])

    #----------------------------------------------------
    def go_aggregate(self):
	
	patt = """r'gene_?id *[=:]? *"?([^"; ]+)"?'"""
	expr = '"GeneID:"+re.search(' + patt + ', IN[9], re.I).group(1)'
	
	# Select rows from first file, and extract gene ids.
	f1 = self.mkTmp()
	self.execStep(TF, [
	    "--file1="+self.options.file1,
	    "?len(%s)==0 or string.lower(IN[3]) in %s" % (self.types1,self.types1),
	    "IN[1:9]",
	    expr,
	    "--out-file=" + f1
	    ])

	# Select rows from second file, and extract gene ids.
	f2 = self.mkTmp()
	self.execStep(TF, [
	    "--file1="+self.options.file2,
	    "?len(%s)==0 or string.lower(IN[3]) in %s" % (self.types2,self.types2),
	    "IN[1:9]",
	    expr,
	    "--out-file=" + f2
	    ])

	# Find the unique genes in file1 and count the exons
	genes1=os.path.join(self.options.outDir, "genes1.txt")
	self.execStep(TA, [
	    "--file1="+f1,
	    "-g9",
	    "-acount",
	    "-afirst:1",
	    "-amin:4",
	    "-amax:5",
	    "-afirst:7",
	    "--out-file=" + genes1
	    ])

	# Find the unique genes in file2 and count the exons
	genes2=os.path.join(self.options.outDir, "genes2.txt")
	self.execStep(TA, [
	    "--file1="+f2,
	    "-g9",
	    "-acount",
	    "-afirst:1",
	    "-amin:4",
	    "-amax:5",
	    "-afirst:7",
	    "--out-file=" + genes2
	    ])

	# Find all overlapping feature pairs.
	ovlExons=self.mkTmp()
	args = [
	    "-1", f1,
	    "-2", f2,
	    "-o", ovlExons,
	    "-s", "both",
	    "-k", self.minOverlap,
	    ]
	if self.options.ignoreStrand:
	    #args.append( "--ignore-strand" )
	    args = args + ["--columns1", "1,4,5"]
	novl = self.execStep(FJ, args)

	if self.options.noSelfHits:
	    xxx = ovlExons
	    ovlExons = self.mkTmp()
	    args = [ 
	    	"?IN[10] != IN[19]",
		"-1", xxx,
		"-o", ovlExons,
		]
	    novl = self.execStep(TF, args)

	if novl == 0:
	    self.debug("No overlapping features detected. Pipeline terminated.\n")
	    sys.exit(-1)

	#  Aggregate overlapping exon pairs into
	#    overlapping gene pairs. Count the exons involved.
	ovlGenes = self.mkTmp()
	self.execStep(TA, [
	    "-g10,19",
	    "-acount:5",
	    "-acount:14",
	    "--file1="+ovlExons,
	    "--out-file=" + ovlGenes
	    ])

	#  Join with genes1 to pull in total exon counts.
	#    Do an outer join so that every gene in genes1 is
	#    represented.
	tmp1 = self.mkTmp()
	self.execStep(TJ, [
	    "--file1=" + ovlGenes,
	    "--file2=" + genes1,
	    "--k1=1",
	    "--k2=1",
	    "--right",
	    "-n.",
	    "--out-file=" + tmp1
	    ])

	#  Join with genes2 to pull in total exon counts.
	#    Do a bidi-outer join so that every gene in genes2 is
	#    represented, as is every gene in genes1..
	tmp2 = self.mkTmp()
	self.execStep(TJ, [
	    "--file1=" + tmp1,
	    "--file2=" + genes2,
	    "--k1=2",
	    "--k2=1",
	    "--left",
	    "--right",
	    "-n.",
	    "--out-file=" + tmp2
	    ])

	#  Filter for final output formatting.
	tmp3 = self.mkTmp()
	self.execStep(TF, [
	    "--file1=" + tmp2,
	    "--out-file=" + tmp3,
	    "IN[7]=='.' and IN[13] or IN[7]" ,
	    #"IN[10]=='.' and IN[16] or IN[10]" ,
	    "IN[10]=='.' and IN[16] or IN[16]=='.' and IN[10] or IN[16]==IN[10] and IN[16] or '???'",
	    "IN[5]", "IN[3]", "IN[6]", "IN[8]", "IN[9]", "IN[10]",
	    "IN[11]", "IN[4]", "IN[12]", "IN[14]", "IN[15]", "IN[16]",
	    ])

	#  Bucketize the overlapping genes. Output separate file for
	# each bucket.
	self.execStep(TB, [
	    "--file1=" + tmp3,
	    "--k1=3",
	    "--k2=9",
	    "-n.",
	    "-o" + self.options.outDir,
	    "-t" + self.options.template,
	    "IN[1:3]", "IN[4:]",   ## remove the bucket id column
	    ])

    def go(self):
	self.debug("======================================\n")
	self.debug("Starting GU pipeline\n", True)
	self.debug("Command line:\n%s\n" % (" ".join(self.argv)))

	if(self.options.noAggregate):
	    self.go_noAggregate()
	else:
	    self.go_aggregate()

	self.debug("Pipeline completed.\n", True)
	self.debug("Cleaning up...\n")

	#  Delete the temp files
	self.cleanupTempFiles()
	self.debug("Goodbye.\n\n")

#--------------------------------
try:
    GUPipeline(sys.argv).go()
except:
    sys.exit(-1)

