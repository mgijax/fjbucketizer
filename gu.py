#
# gu - Gene Unification Pipeline
#
# This pipeline compares two GFF input files for coordinate overlaps, and
# bucketizes the overlapping feature pairs. The comparison step takes two lists 
# of features [a1, a2, ...], [b1, b2, ...], outputs all feature pairs [ai, bj] such
# that ai overlaps bj.
#
# These feature pairs form the edges in a bipartite graph of overlapping
# a's and b's. Bucketization finds the connected components of the graph
# and segregates those components (or clusters) into buckets, according to 
# how many a's and b's there are in the cluster (one-to-one, one-to-many, etc).
# The buckets are output to separate files (bucket_1-1.txt, bucket_1-n.txt, ...)
#
# There is one line output per overlapping pair [ai,bj]. This line contains
# all the columns of ai and bj plus three additional columns about the cluster
# to which this pair belongs. In detail:
#       1. The id of the cluster (connected component) to which this pair belongs.
#       2. The bucket label (e.g. "n-1").
#       3. The actual counts of a's and b's in the cluster, e.g. "3-1" for 3 a's and 1 b.
#       4- All of the columns of ai followed by all of the columns of bj.
#
# GU operates in one of two modes. The above description corresponds to 'no-aggregate'
# mode: the overlapping features are bucketized directly. Sometimes, however, you want
# the results aggregated to a higher level. For example, you may have files of exons,
# but what you want is buckets of overlapping genes. For this, there is 'aggregate' mode,
# in which lists of overlapping exons (for example) are turned into lists of overlapping genes
# before being bucketized. Because aggregation loses the details of the underlying features,
# the output format is somewhat different from no-aggregate mode. 
#
# Sample invocation:
#     python3 gu.py -a -C --f1 /data/research/yz/New_GU/gff3/MGI.exome.gff3 --f2 /data/research/yz/ENSEMBL/mouse_99/gff3/ENSEMBL_all.gff3 --are1="MGI:[0-9]+" --are2="ENSEMBL:ENSMUSG[0-9]+"
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
import re

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

        self.parser.add_option("--ire", 
            dest="ire",
            metavar="REGEX",
            default=r'(ID=)?(?P<id>[^;]+);?',
            help="Id regular expression. A regular expression that will be used to extract ids from column 9. " +
              "(Also: what you'll feel as you try to get this parameter right...) " +
              "Only applicable if you specify -s (no self hits). " +
              "If the same feature can appear in both input files, it will naturally overlap itself. " +
              "To filter out such instances, GU needs to know how to find the feature ids in column 9. " +
              "This is done with regular expression patterns. " +
              "With -ire, you specify a pattern that will be used for both input files. " +
              "You can specify different patterns for " +
              "files 1 and 2 by specifying --ire1 and --ire2 instead. " +
              "By default, GU looks for patterns like 'ID=blah'. If your feature ids look like this, you " +
              "do not need to do anything. If not, read on. " +
              "Example: to aggregate both inputs by MGI id: --ire 'MGI:[0-9]+'. " +
              "Example: to aggregate input1 by MGI id and input 2 by VEGA id: " +
              "--ire1 'MGI:[0-9]+'  --ire2 'OTTMUSG[0-9]+'. " +
              "Advanced usage: Sometimes you need to define a regular expression where the actual " +
              "id you wish to extract is a sub-part of the whole pattern. This is called a 'capture'. " +
              "To capture the id, surround that part of the pattern with the magic symbols " +
              "'(?P<id>' and ')'. For example, suppose you want to capture the MGI id only when it is " +
              "part of a dbxref attribute like this: 'Dbxref=MGI:MGI:012345;'. You could use the following " +
              "regular expression: 'Dbxref=MGI:(?P<id>MGI:[0-9]+);'."
            )

        self.parser.add_option("--ire1", 
            dest="ire1",
            metavar="REGEX",
            default=None,
            help="Specify regex for input 1 only. See --ire."
            )

        self.parser.add_option("--ire2", 
            dest="ire2",
            metavar="REGEX",
            default=None,
            help="Specify regex for input 2 only. See --ire."
            )

        self.parser.add_option("--are", 
            dest="are",
            metavar="REGEX",
            default=r'Parent=(?P<id>[^"; ]+)',
            help="Aggregation regular expression. "
            )

        self.parser.add_option("--are1", 
            dest="are1",
            metavar="REGEX",
            default=None,
            help="Specify aggregation regex for input 1 only. See --are."
            )

        self.parser.add_option("--are2", 
            dest="are2",
            metavar="REGEX",
            default=None,
            help="Specify aggregation regex for input 2 only. See --are."
            )

        self.parser.add_option("-i", "--ignoreStrand", 
            dest="ignoreStrand",
            action="store_true",
            default = False,
            help="Ignore strand when determining overlaps (Default = strands must match)")

        self.parser.add_option("-a", "--aggregate", 
            dest="aggregate",
            action="store_true",
            default = False,
            help="Aggregate overlapping features by id. Opposite of -n. See also --are, --are1, and --are2.")

        self.parser.add_option("-n", "--noAggregate", 
            dest="noAggregate",
            action="store_true",
            default = False,
            help="Do not aggregate features. Opposite of -a.")

        self.parser.add_option("-C", "--chrMatchLoose", 
            dest="chrMatchLoose",
            action="store_true",
            default = False,
            help="If specified, chromosome matching is 'loose'. Otherwise it is exact. " + \
                "In loose matching, leading 'chr' is removed from chromosome field, " + \
                "so that '19' matches 'Chr19'. ")

        self.parser.add_option("-s", "--noSelfHits", 
            dest="noSelfHits",
            action="store_true",
            default = False,
            help="Ignore self-hit in overlap detection. (default: reports self-hits)")

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

        if not (self.options.aggregate or self.options.noAggregate):
            self.parser.error("Must specify either -n or -a.")
        
        if self.options.aggregate and self.options.noAggregate:
            self.parser.error("Cannot specify both -n and -a.")
        
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
        
        if self.options.chrMatchLoose:
            self.options.chrMatchExpr = "[re.sub('^(?i)chr(om(osome)?)?', '', IN[1])]+IN[2:]"
        else:
            self.options.chrMatchExpr = "IN[1:]"
        
        self.options.guDir = os.path.split(__file__)[0]
        self.options.guUtilFile = os.path.join(self.options.guDir,'guUtil.py')


        # local function
        def groomRe(name, default):
            r = getattr(self.options, name)
            if r is None:
                r = default
            if "(?P<id>" not in r:
                r = "(?P<id>%s)" % r
            r = "r'%s'" % r
            setattr(self.options, name, r )
            self.debug('%s=%s\n'%(name,r))

        groomRe('ire1',self.options.ire)
        groomRe('ire2',self.options.ire)

        groomRe('are1',self.options.are)
        groomRe('are2',self.options.are)

        
    #----------------------------------------------------
    def mkTmp(self, preserve=False):
        tf = mkstemp(dir=self.options.outDir)
        os.close(tf[0])
        if not preserve:
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
        self.debug( tool.__name__ + " " + " ".join(args) + "\n", True)
        t = tool(args)
        t.go()
        t.closeFiles()
        return t.nOutputRows

    #----------------------------------------------------
    def go_noAggregate(self):

        ### Select rows from first file.

        args = [
            "--file1="+self.options.file1,
            self.options.chrMatchExpr,
            ]

        if len(self.options.types1) > 0:
            args.append( "?str.lower(IN[3]) in %s" % self.types1 )

        if len(self.options.notTypes1) > 0:
            args.append( "?str.lower(IN[3]) not in %s" % self.notTypes1 )

        if len(args) > 2:
            f1 = self.mkTmp()
            args.append( "--out-file=" + f1 )
            self.execStep(TF, args)
        else:
            f1 = self.options.file1

        ### Select rows from second file.

        args = [
            "--file1="+self.options.file2,
            self.options.chrMatchExpr,
            ]

        if len(self.options.types2) > 0:
            args.append( "?str.lower(IN[3]) in %s" % self.types2 )

        if len(self.options.notTypes2) > 0:
            args.append( "?str.lower(IN[3]) not in %s" % self.notTypes2 )

        if len(args) > 2:
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

        if novl == 0:
            self.debug("No overlapping features detected.\n")

        ### bucketize the pairs.

        bucketized=self.mkTmp()
        self.execStep(TB, [
            "--file1=" + overlaps,
            "--k1=10",
            "--k2=19",
            "-t"+bucketized,
            "IN[1:]",
            "int(IN[11]!=IN[20])",      ## compute column: 0==same strands 1==diff
            "int(str.lower(IN[7])!='gene')", ## compute column: 0==all genes 1==nongene
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
        
        # Select rows from first file, and extract feature ids.
        f1 = self.mkTmp()
        args = [
            "--file1="+self.options.file1,
            "--exec-file="+self.options.guUtilFile,
            ]
        if len(self.options.types1) > 0:
            args.append( "?str.lower(IN[3]) in %s" % self.types1 )
        if len(self.options.notTypes1) > 0:
            args.append( "?str.lower(IN[3]) not in %s" % self.notTypes1 )
        args += [
            "IN[1:9]",
            'extractID(IN[9],%s)'%self.options.are1,
            "--out-file=" + f1
            ]

        self.execStep(TF, args)

        # Select rows from second file, and extract feature ids.
        f2 = self.mkTmp()
        args = [
            "--file1="+self.options.file2,
            "--exec-file="+self.options.guUtilFile,
            ]
        if len(self.options.types2) > 0:
            args.append( "?str.lower(IN[3]) in %s" % self.types2 )
        if len(self.options.notTypes2) > 0:
            args.append( "?str.lower(IN[3]) not in %s" % self.notTypes2 )
        args += [
            "IN[1:9]",
            'extractID(IN[9],%s)'%self.options.are2,
            "--out-file=" + f2
            ]
        self.execStep(TF, args)

        # Find the distinct higher-level features in file1 and count the base features
        genes1=os.path.join(self.options.outDir, "features1.txt")
        self.execStep(TA, [
            "--file1="+f1,
            "-g9",              # id
            "-acount",          # num. lines w/ this id
            "-afirst:1",        # first chr
            "-amin:4",          # min start val
            "-amax:5",          # max end val
            "-afirst:7",        # first strand
            "--out-file=" + genes1
            ])

        # Find the unique genes in file2 and count the exons
        genes2=os.path.join(self.options.outDir, "features2.txt")
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
            self.debug("No overlapping features detected.\n")

        #  Aggregate overlapping feature pairs into higher-level overlaps.
        #  Count the base features involved.
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
            "-t" + self.options.template,
            "IN[1:3]", "IN[4:]",   ## remove the bucket id column
            ])

    def go(self):
        self.debug("======================================\n")
        self.debug("Starting GU pipeline\n", True)
        self.debug("Command line:\n%s\n" % (" ".join(self.argv)))

        if(self.options.aggregate):
            self.go_aggregate()
        else:
            self.go_noAggregate()

        self.debug("Pipeline completed.\n", True)
        self.debug("Cleaning up...\n")

        #  Delete the temp files
        self.cleanupTempFiles()
        self.debug("Goodbye.\n\n")

#--------------------------------
GUPipeline(sys.argv).go()

