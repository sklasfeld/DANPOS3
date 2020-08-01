#!/usr/bin/env python
#2.2.2  version

import os,sys
from time import time
from functions import danpos, div
from math import log10
from copy import deepcopy
from wig import Wig
from wigs import Wigs
import argparse,sys,os
from lib import positionSelectorByGreatTSS,positionDicMinMax,vioplot,positionDistance,positionSelectorByValue,positionSelectorByGeneStructure,batchOccAroundPoints,batchOccInRegions,occAroundPoints,plot,batchOccPSD,retrieve_positions_by_value,batchPositionDistanceDistribution,batchPositionValDistribution
from wiq import rawsort,refquantile,changevalue,qnorwig,wiq,wig2wiq
from rpy2.robjects import r

#sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # This allow DANPOS to print each message on screen immediately.
class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def writelines(self, datas):
       self.stream.writelines(datas)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)

sys.stdout = Unbuffered(sys.stdout)

def printHelp():
    print('\ndanpos version 3.0.0')
    print('For help information for each function, try:\npython danpos.py <function> -h')
    print('\nFunctions:')
    print('\tdpos:\n\t\tanalyze each protein-binding position (~100\n\t\tto ~200bp wide) across the whole genome,\n\t\te.g. nucleosome positions.')
    print('\tdpeak:\n\t\tanalyze each protein-binding peak (~1 to ~1kb\n\t\twide) across the whole genome, e.g. protein\n\t\tthat binds accruately to some specific motifs.')
    print('\tdregion:\n\t\tanalyze each protein-binding region (~1 to\n\t\t~10kb wide) across the whole genome, e.g.\n\t\tsome histone modifications.')
    print('\tdtriple:\n\t\tDo analysis at all three levels including each\n\t\tregion, peak, and position. Would be useful\n\t\twhen little is known about the potential binding\n\t\tpattern.')
    print('\tprofile:\n\t\tanalyze wiggle format occupancy or differential\n\t\tsignal profile relative to gene structures or\n\t\tbed format genomic regions.')
    print('\twiq:\n\t\tnormalize wiggle files to have the same quantiles (Quantile normalization).')
    print('\twig2wiq:\n\t\tconvert wiggle format file to wiq format.')
    print('\tstat:\n\t\tsome statistics for positions, peaks, or regions.')
    print('\tselector:\n\t\tselect a subset of positions, peaks, or regions\n\t\tby value ranges or gene structures neighboring.')
    print('\tvaluesAtRanks:\n\t\tretrieve position, peak, or region values by ranks.')
    print('\nKaifu Chen, et al. chenkaifu@gmail.com, Li lab, Biostatistics department, Dan L. Duncan cancer center, Baylor College of Medicine.')
    print('')

def runDANPOS(command=''):
    '''
    Description:
        this function provide an entrance to the package DANPOS.
        It parse input parameters, print some help messages if required, else it call and pass all parameters values to the main function danpos().
    
    parameters:
        none  
    '''
    if command=='dpos':tname='position'
    elif command=='dtriple':tname='position, peak, and region'
    else:tname=command
    if (len(sys.argv)<3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least one parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\n\npython danpos.py",command, "<path> [optional arguments]\n\nfor more help, please try: python danpos.py "+command+" -h\n")
        return 0
    
    # parse all input parameters
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage="\n\npython danpos.py  <command>  <path> [optional arguments]\n\n",\
                                     description='',epilog="Kaifu Chen, et al. chenkaifu@gmail.com,  \
                                     Li lab, Biostatistics department, Dan L. Duncan cancer center, \
                                     Baylor College of Medicine.")
    #input/output
    parser.add_argument('command',default=None,\
                        help="set as '"+command+"' to run analysis for each "+tname+'.')
    parser.add_argument('path',default=None,\
                        help="Pairs of paths to sequencing data sets. \
                        The two paths in each pair must be seperated by ':', \
                        a:b means a minus b,  \
                        different pairs must be seperated by ',' e.g. file1.bed:dir2/,dir3/, \
                        each path could point to a file or a directory containing multiple files, \
                        files under each directory represent multiple replicates for the same group. \
                        suggest to use .sam or .bam format for input files, please read the documentation for \
                        details about the other supported input formats.")
    parser.add_argument('--------------------------        ',dest="separator1",metavar='',default=None,help="")
    parser.add_argument('--- general parameters ---        ',dest="separator1",metavar='',default=None,help="")
    parser.add_argument('--------------------------         ',dest="separator1",metavar='',default=None,help="")
    parser.add_argument("-m","--paired", dest="paired",metavar='',default=0,type=int,\
                        help="set to 1 if the input data is mate-pair (paired-end) reads. Ignore this when the input is wiggle format occupancy data ")
    if command=='dtriple':
        parser.add_argument('-p', '--pheight',dest="pheight",metavar='',default='1e-10',\
                            help="occupancy/intensity P value cutoff for calling invidual binding region/peak/position")
        parser.add_argument('-q', '--height',dest="height",metavar='',default=0,type=int,\
                            help="occupancy/intensity cutoff for calling invidual binding region/peak/position")
        parser.add_argument('-t', '--testcut',dest="testcut",metavar='',default='0',\
                            help="P value cutoff for calling differential region/peak/position between samples (e.g. 1e-10). \
                            Set as 0 when don't need to define peaks based on differential P value, \
                            so region/peak/position will then be defined only on occupancy cutoff, \
                            and differential P value will be calculated for each of them.")
    if command=='dregion':
        parser.add_argument('-p', '--pheight',dest="pheight",metavar='',default='1e-10',\
                            help="occupancy/intensity P value cutoff for calling invidual binding region")
        parser.add_argument('-q', '--height',dest="height",metavar='',default=0,type=int,\
                            help="occupancy/intensity cutoff for calling invidual binding region")
        parser.add_argument('-t', '--testcut',dest="testcut",metavar='',default='0',\
                            help="P value cutoff for calling differential region between samples (e.g. 1e-10). \
                            Set as 0 when don't need to define peaks based on differential P value, \
                            so regions will then be defined only on occupancy cutoff, \
                            and differential P value will be calculated for each of them.")
    if command=='dpeak':
        parser.add_argument('-p', '--pheight',dest="pheight",metavar='',default='1e-10',\
                            help="occupancy/intensity P value cutoff for calling invidual binding peak")
        parser.add_argument('-q', '--height',dest="height",metavar='',default=0,type=int,\
                            help="occupancy/intensity cutoff for calling invidual binding peak")
        parser.add_argument('-t', '--testcut',dest="testcut",metavar='',default='0',\
                            help="P value cutoff for calling differential peak between samples (e.g. 1e-10). \
                            Set as 0 when don't need to define peaks based on differential P value, \
                            so peaks will then be defined only on occupancy cutoff, \
                            and differential P value will be calculated for each of them.")
    if command=='dpos':
        parser.add_argument('-p', '--pheight',dest="pheight",metavar='',default='0',\
                            help="occupancy/intensity P value cutoff for calling invidual binding position")
        parser.add_argument('-q', '--height',dest="height",metavar='',default=5,type=int,\
                            help="occupancy/intensity cutoff for calling invidual binding position")
        parser.add_argument('-t', '--testcut',dest="testcut",metavar='',default='0',\
                            help="P value cutoff for calling differential position between samples (e.g. 1e-10). \
                            Set as 0 when don't need to define positions based on differential P value, \
                            so binding positions will then be defined only on occupancy cutoff, \
                            and differential P value will be calculated for each of them.")
    parser.add_argument('-o','--out', dest="name",metavar='',default='result',\
                        help="a name for the output directory")
    parser.add_argument('-f', '--fdr',dest="fdr",metavar='',default=1,type=int,\
                        help="set to 0 if need not to calculate FDR values (slow process).")
    parser.add_argument('-s', '--save',dest="save",metavar='',default=0,type=int,\
                        help="save middle stage files? set to 0 if don't save, otherwise set to 1")
    parser.add_argument('-b','--bg',dest='bg',metavar='',default=None,\
                        help="pairs of paths, each pair secify a genomic background data set for a MNase-/ChIP-Seq data set, \
                        a:b means b is the background data set for a, \
                        put a word 'None' when a MNase-/ChIP-Seq data set has no background data set,\
                        e.g. file1.bed:bgdir1,dir2/:None,dir3/:bg3.bed, \
                        this function is not recommended for MNase-Seq data set.\
                        ")
    if command=='dregion' or command=='dtriple':
        #region calling
        parser.add_argument('--------------------------',dest="separator1",metavar='',default=None,help="")
        parser.add_argument('---   region calling   ---',dest="separator1",metavar='',default=None,help="")
        parser.add_argument('-------------------------- ',dest="separator1",metavar='',default=None,help="")  
        if command=='dtriple':
                parser.add_argument('-r', '--call_region',dest="call_region",metavar='',default=1,type=int,\
                                help="Set to 0 if need not to analyze each invidual binding region.")
        parser.add_argument('-rd', '--region_dis',dest="region_dis",metavar='',default=3000,type=int,\
                            help="minimal tail-to-head distance between neighboring regions (bp), neighboring regions closer than -D will be merged as one single region")
        parser.add_argument('-rw', '--region_width',dest="region_width",metavar='',default=40,type=int,\
                            help="minimal width of each region")
        parser.add_argument('-rf', '--region_reference',dest="region_reference",metavar='',default=None,\
                            help="Don't call regions, but retrive values for a set of reference regions provided in a region file by this parameter.")
    if command=='dpeak' or command=='dtriple':
        #region calling
        parser.add_argument('--------------------------              ',dest="separator1",metavar='',default=None,help="")
        parser.add_argument('---    peak calling    ---              ',dest="separator1",metavar='',default=None,help="")
        parser.add_argument('--------------------------            ',dest="separator1",metavar='',default=None,help="")
        if command=='dtriple':
            parser.add_argument('-k', '--call_peak',dest="call_peak",metavar='',default=1,type=int,\
                                help="Set to 0 if need not to analyze each invidual peak.")
        parser.add_argument('-kd', '--peak_dis',dest="peak_dis",metavar='',default=40,type=int,\
                            help="minimal tail-to-head distance (bp) between neighboring peaks, neighboring peaks closer than -D will be merged as one single peak")
        parser.add_argument('-kw', '--peak_width',dest="peak_width",metavar='',default=40,type=int,\
                            help="minimal width of each peak")
        parser.add_argument('-kf', '--peak_reference',dest="peak_reference",metavar='',default=None,\
                            help="Don't call peaks, but retrive values for a set of reference peaks provided in the peak file by this parameter.")
    if command=='dpos' or command=='dtriple': 
        # position calling
        parser.add_argument('--------------------------  ',dest="separator1",metavar='',default=None,help="")
        parser.add_argument('---  Position calling  ---',dest="separator1",metavar='',default=None,help="")
        parser.add_argument('--------------------------   ',dest="separator1",metavar='',default=None,help="")  
        if command=='dtriple':
                parser.add_argument('-j', '--call_pos',dest="call_pos",metavar='',default=1,type=int,\
                                help="Set to 0 if need not to analyze each invidual binding position.")
        parser.add_argument('-jw', '--width',dest="width",metavar='',default=40,type=int,\
                            help="the window size used for scanning for the summit of each potential position.")
        parser.add_argument('-jd', '--distance',dest="distance",metavar='',default=100,type=int,\
                            help="minimal center-to-center distance between positions, positions closer than d will be merged as one single position")
        parser.add_argument('-jf', '--position_reference',dest="position_reference",metavar='',default=None,\
                            help="map each defined position to a reference position provided in the position file by this parameter.")
        parser.add_argument('-R', '--ratio',dest="ratio",metavar='',default=0.90,\
                            help="the ratio between the minimal occupancy flanking a position and the maximal occupancy in a position, only positions with values lower than this ratio will be used for defining position shift events.")
        parser.add_argument('-e', '--edge',dest="edge",metavar='',default=0,type=int,\
                            help="set to 1 if need to detect edges for each position,else set to 0")
        parser.add_argument('-g', '--gapfill',dest="gapfill",metavar='',default=0,type=int,\
                            help="do gap filling? fill the gap between two neighboring positions with an additional position \
                            if the gap size is close to the a position size, set to 0 if don't fill, otherwise set to 1")
    
    # occupancy processing
    parser.add_argument('--------------------------    ',dest="separator1",metavar='',default=None,help="")
    parser.add_argument('---occupancy processing---',dest="separator1",metavar='',default=None,help="")
    parser.add_argument('--------------------------     ',dest="separator1",metavar='',default=None,help="")
    parser.add_argument('-c','--count',dest='count',metavar='',default=None,\
                        help="specify the count of reads to be normalized to, e.g. 10000000. \
                        Or specify the count for each group, e.g. \
                        file1.bed:10000000,dir2/:20000000,dir3/:15000000, \
                        Do this only when you are clear about what you are doing, e.g. \
                        when you have spike-ins to measure the real reads count in each replicate" )
    parser.add_argument('-a', "--span", dest="span",metavar='',default=10,type=int,\
                        help="the span or step size  in the generated wiggle data")
    parser.add_argument('-z', '--smooth_width', dest="smooth_width",metavar='',default=20,type=int,\
                        help="the smooth width before position calling, set to 0 if need not to smooth")
    parser.add_argument("-L", "--exclude_low_percent",dest="exclude_low_percent",metavar='',default=0,type=float,\
                        help="the percent of extremely low occupancy positions to be excluded in determining normalization factors,\
                        may be helpful to avoid the influence of background noise on normalization")
    parser.add_argument("-H",  "--exclude_high_percent",dest="exclude_high_percent",metavar='',default=0,type=float,\
                        help="the percent of extremely high occupancy positions to be excluded in determining normalization factors,\
                        may be helpful to avoid the influence of some highly clonal or repeat regions on normalization.")
    parser.add_argument('-l', '--lmd',dest="lmd",metavar='',default=300,type=int,\
                        help="lambda width for smoothing background data before background subtraction,\
                        ignore this when the parameter -b is not specified.")
    parser.add_argument('-n', '--nor',dest="nor",metavar='',default='F',\
                        help="data normalization method, could be 'F','S' or 'N', \
                        representing normalization by fold change, normalize by sampling, or no normalization")
    parser.add_argument('-N', '--nor_region_file',dest="nor_region_file",metavar='',default=None,\
                        help="A '.wig' format file to denote the regions that could be used to calculate normalization factors.\
                        Regions to be used and not used should be assigned a value 1 and 0 in this .wig file, respectively.")
    
    parser.add_argument('--nonzero',dest="nonzero",metavar='',default=0,type=int,\
                        help="set to 1 if want to normalize basepairs with non-zero values to have the same average value between different data sets.\
                        This function will be useful when some data sets has severious clonal effect, \
                        E.g. one data set has non-zero value at 10 percent of base pairs and each non-zero vase pair has 8 fold clonal effects, \
                        whereas another data set has non-zero value at 40 percent of base pairs and each non-zero base pair has 2 fold clonal effects.")
    '''
    parser.add_argument('-s', '--statis',dest="statis",metavar='',default='P',\
                        help="the statistics method for differential test, could be 'C','P', 'F', or 'S', \
                        representing Chi-square test, Possion test, Fold change, or direct subtraction")
    '''
    
    
    # reads processing
    parser.add_argument('--------------------------      ',dest="separator1",metavar='',default=None,help="")
    parser.add_argument('---  reads processing  ---',dest="separator1",metavar='',default=None,help="")
    parser.add_argument('--------------------------       ',dest="separator1",metavar='',default=None,help="")
    parser.add_argument('-u',"--clonalcut", dest="clonalcut",metavar='',default=0,\
                        help="the cutoff for adjusting clonal signal, \
                        set as a P value larger than 0 and smaller than 1, e.g 1e-10,\
                        or set as a interger reads count,\
                        set as 0 if don't need to adjust clonal signal.", type=float)
    parser.add_argument("--frsz", dest="fs",metavar='',default=None,type=int,\
                        help="specify the average size of \
                        DNA fragments in the seuqnecing experiment. By default it is automatically \
                        detected by DANPOS. Ignore this when the input is wiggle format occupancy data")
    parser.add_argument("--mifrsz", dest="mifrsz",metavar='',default=50,type=int,\
                        help="minimal size of the DNA fragments, DANPOS will select a most \
                        probable frsz value within the range between --mifrsz value and --mafrsz value. \
                        ignore this when '--frsz' has been specified. Ignore this when the input is wiggle format occupancy data")
    parser.add_argument("--mafrsz", dest="mafrsz",metavar='',default=300,type=int,\
                        help="maximal size of the DNA fragments, DANPOS will select a most \
                        probable frsz value within the range between --mifrsz value and --mafrsz value. \
                        ignore this when '--frsz' has been specified. Ignore this when the input is wiggle format occupancy data")
    parser.add_argument("--extend", dest="extend",metavar='',default=80,type=int,\
                        help="specify the size theat each fragment will be adjusted to,\
                        the size of each fragment will be adjusted to this size when reads data\
                        is converted to occupancy data. Ignore this when the input is wiggle format occupancy data")
    
    
    
    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print("\ndanpos 2.2.2  version\n")
        parser.print_help()
        print("\n")
        return 0
    elif len(sys.argv)>=3: # at least one parameter need to be specified: a path to input file or files
        try:
            args=parser.parse_args()  #all paramter values are now saved in args
        except:
            print("\nfor more help, please try: python danpos.py",command,"-h\n")
            return 0

    if args.exclude_low_percent!=0 and args.nor_region_file!=None:
        print("\nplease don't define both -L (--exclude_low_percent) and -N (--nor_region_file) in the command line\n")
        return
    if args.exclude_high_percent!=0 and args.nor_region_file!=None:
        print("\nPlease don't define both -H (--exclude_high_percent) and -N (--nor_region_file) in the command line\n")
        return
    
    print("\ndanpos 3.0.0  version\n")
    print('command:\npython'," ".join(sys.argv)) # print the command line, this let the user to keep a log and remember what parameters they specified
    print('\n',args) # print all parameter values, this provide a eacy way for the user to double check the parameter values used by DANPOS
    
    from math import log10
    if True:#args.statis=='C' or args.statis=='P':
        temptestcut=args.testcut.split('-')
        if len(temptestcut)==2:testcut=float(temptestcut[1])-log10(float(temptestcut[0][:-1]))
        elif float(temptestcut[0])>0:testcut=0-log10(float(temptestcut[0]))
        else:testcut=0
    
    if True:#args.statis=='C' or args.statis=='P':
        temppheight=args.pheight.split('-')
        if len(temppheight)==2:pheight=float(temppheight[1])-log10(float(temppheight[0][:-1]))
        elif float(temppheight[0])>0:pheight=0-log10(float(temppheight[0]))
        else:pheight=0
    nonzero=0
    if args.nonzero==0:nonzero=0
    elif args.nonzero==1:nonzero=1
    else:
        print("'--nonzero' must be ste to either 0 or 1")
        return
    #print testcut,pheight
    
    if command=='dtriple':
        danpos(\
               #input output args
               tpath=args.path,paired=args.paired,opath=args.name,save=args.save,tbg=args.bg,fdr=args.fdr,\
               #region calling
               call_region=args.call_region,region_width=args.region_width,region_distance=args.region_dis,\
               call_peak=args.call_peak,peak_width=args.peak_width,peak_distance=args.peak_dis,\
               #position calling args
               call_position=args.call_pos,height=args.height,pheight=pheight,logp=testcut,width=args.width,distance=args.distance,edge=args.edge,fill_gap=args.gapfill,ratio=args.ratio,\
               #occupancy processing args
               ref_region=args.position_reference,ref_peak=args.peak_reference,ref_position=args.position_reference,\
               nor=args.nor,nonzero=nonzero,amount=args.count,step=args.span,smooth_width=args.smooth_width,lmd=args.lmd,\
               #reads processing
               cut=args.clonalcut,fs=args.fs,mifrsz=args.mifrsz,mafrsz=args.mafrsz,extend=args.extend,\
               # whether to bo position calling for both gain and loss of occupancy
               #exclude_low_percent=args.exclude_low_percent,exclude_high_percent=args.exclude_high_percent,\
               exclude_low_percent=args.exclude_low_percent,exclude_high_percent=args.exclude_high_percent,region_file=args.nor_region_file,\
               both=True,\
               # P value cutoff for calling occupancy positions, set to 1 will disable this parameter
               #pheight=1,\
               #the output format of wiggle format file
               wgfmt='fixed',\
               #the default value for gap filling in position calling
               fill_value=1,\
               #the differential test method,'P':Poisson, 'C':Chi-Square
               test='P',\
               #whether or not to do position calling for each replicate, set to 1 if need to do, else set to 0.
               pcfer=0)#args.pcfer)#0)
    if command=='dregion':
        danpos(\
               #input output args
               tpath=args.path,paired=args.paired,opath=args.name,save=args.save,tbg=args.bg,fdr=args.fdr,\
               #region calling
               call_region=1,region_width=args.region_width,region_distance=args.region_dis,\
               call_peak=0,\
               #peak_width=args.peak_width,peak_distance=args.peak_dis,\
               #position calling args
               call_position=0,\
               ref_region=args.region_reference,\
               #width=args.width,distance=args.distance,edge=args.edge,fill_gap=args.gapfill
               height=args.height,pheight=pheight,logp=testcut,\
               #occupancy processing args
               nor=args.nor,nonzero=nonzero,amount=args.count,step=args.span,smooth_width=args.smooth_width,lmd=args.lmd,\
               #reads processing
               cut=args.clonalcut,fs=args.fs,mifrsz=args.mifrsz,mafrsz=args.mafrsz,extend=args.extend,\
               # whether to bo position calling for both gain and loss of occupancy
               #exclude_low_percent=args.exclude_low_percent,exclude_high_percent=args.exclude_high_percent,\
               exclude_low_percent=args.exclude_low_percent,exclude_high_percent=args.exclude_high_percent,region_file=args.nor_region_file,\
               both=True,\
               # P value cutoff for calling occupancy positions, set to 1 will disable this parameter
               #pheight=1,\
               #the output format of wiggle format file
               wgfmt='fixed',\
               #the default value for gap filling in position calling
               fill_value=1,\
               #the differential test method,'P':Poisson, 'C':Chi-Square
               test='P',\
               #whether or not to do position calling for each replicate, set to 1 if need to do, else set to 0.
               pcfer=0)#args.pcfer)#0)
    if command=='dpeak':
        danpos(\
               #input output args
               tpath=args.path,paired=args.paired,opath=args.name,save=args.save,tbg=args.bg,fdr=args.fdr,\
               #region calling
               call_region=0,\
               #peak calling
               call_peak=1,peak_width=args.peak_width,peak_distance=args.peak_dis,\
               #position calling args
               call_position=0,\
               ref_peak=args.peak_reference,\
               #width=args.width,distance=args.distance,edge=args.edge,fill_gap=args.gapfill,\
               height=args.height,pheight=pheight,logp=testcut,\
               #occupancy processing args
               nor=args.nor,nonzero=nonzero,amount=args.count,step=args.span,smooth_width=args.smooth_width,lmd=args.lmd,\
               #reads processing
               cut=args.clonalcut,fs=args.fs,mifrsz=args.mifrsz,mafrsz=args.mafrsz,extend=args.extend,\
               # whether to bo position calling for both gain and loss of occupancy
               #exclude_low_percent=args.exclude_low_percent,exclude_high_percent=args.exclude_high_percent,\
               exclude_low_percent=args.exclude_low_percent,exclude_high_percent=args.exclude_high_percent,region_file=args.nor_region_file,\
               both=True,\
               # P value cutoff for calling occupancy positions, set to 1 will disable this parameter
               #pheight=1,\
               #the output format of wiggle format file
               wgfmt='fixed',\
               #the default value for gap filling in position calling
               fill_value=1,\
               #the differential test method,'P':Poisson, 'C':Chi-Square
               test='P',\
               #whether or not to do position calling for each replicate, set to 1 if need to do, else set to 0.
               pcfer=0)#args.pcfer)#0)    
    if command=='dpos':
        danpos(\
               #input output args
               tpath=args.path,paired=args.paired,opath=args.name,save=args.save,tbg=args.bg,fdr=args.fdr,\
               #region calling
               call_region=0,\
               #peak calling
               call_peak=0,\
               #position calling args
               call_position=1,width=args.width,distance=args.distance,edge=args.edge,fill_gap=args.gapfill,ratio=args.ratio,\
               ref_position=args.position_reference,\
               height=args.height,pheight=pheight,logp=testcut,\
               #occupancy processing args
               nor=args.nor,nonzero=nonzero,amount=args.count,step=args.span,smooth_width=args.smooth_width,lmd=args.lmd,\
               #reads processing
               cut=args.clonalcut,fs=args.fs,mifrsz=args.mifrsz,mafrsz=args.mafrsz,extend=args.extend,\
               # whether to bo position calling for both gain and loss of occupancy
               #exclude_low_percent=args.exclude_low_percent,exclude_high_percent=args.exclude_high_percent,\
               exclude_low_percent=args.exclude_low_percent,exclude_high_percent=args.exclude_high_percent,region_file=args.nor_region_file,\
               both=True,\
               # P value cutoff for calling occupancy positions, set to 1 will disable this parameter
               #pheight=1,\
               #the output format of wiggle format file
               wgfmt='fixed',\
               #the default value for gap filling in position calling
               fill_value=1,\
               #the differential test method,'P':Poisson, 'C':Chi-Square
               test='P',\
               #whether or not to do position calling for each replicate, set to 1 if need to do, else set to 0.
               pcfer=0)#args.pcfer)#0)                  

def wigProfileToGene(args,wgs,pos_neg='', outname='',outmode='w'):
    rcode=''
    #rcode='par(mfrow=c('+str(nrow)+','+str(ncol)+'))\n'#+rcode
    gfiles=args.genefile_paths.split(',')
    wigaliases=args.wigfile_aliases.split(',')
    colors=args.plot_colors.split(',')
    if args.genefile_aliases==None:gfnames=args.genefile_paths.split(',')
    else:gfnames=args.genefile_aliases.split(',')
    sites=args.genomic_sites.split(',')
    if 'TSS' in sites:
        print('\nprofiling for Transcript Start Sites (TSS)')
        d={}
        if len(list(wgs.keys()))>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                if i==0:outmode='w'
                else :outmode='a'
                dic=batchOccAroundPoints(wgs,outname=outname+'_TSS',groupname=gfname+pos_neg+'.tss',outmode=outmode,chrColID=1,nameColID=0,posColIDpos=3,posColIDneg=4,straColID=2,sep='\t',second_sep=None,step=args.bin_size,lines=glines,heatMap=args.heatmap,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                if len(colors)>=len(wigaliases):rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg+'.tss',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                '''
                fo=open(wfname+pos_neg+'.tss.xls','w')
                fo.write('pos\t'+'\t'.join(gfiles)+'\n')
                poses=dic[gfiles[0]].keys()
                poses.sort()
                for i in poses:
                    oline=str(i)
                    for k in gfiles:oline+='\t'+str(dic[k][i])
                    fo.write(oline+'\n')
                '''
                if len(colors)>=len(gfnames):rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg+'.tss',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))

    if 'TTS' in sites:
        print('\nprofiling for Transcription Terminal Sites (TTS)')
        d={}
        if len(list(wgs.keys()))>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                if i==0:outmode='w'
                else :outmode='a'
                dic=batchOccAroundPoints(wgs,outname=outname+'_TTS',groupname=gfname+pos_neg+'.tts',outmode=outmode,chrColID=1,nameColID=0,posColIDpos=4,posColIDneg=3,straColID=2,sep='\t',second_sep=None,step=args.bin_size,lines=glines,heatMap=args.heatmap,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                if len(colors)>=len(wigaliases):rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg+'.tts',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                '''
                fo=open(wfname+pos_neg+'.tts.xls','w')
                fo.write('pos\t'+'\t'.join(gfiles)+'\n')
                poses=dic[gfiles[0]].keys()
                poses.sort()
                for i in poses:
                    oline=str(i)
                    for k in gfiles:oline+='\t'+str(dic[k][i])
                    fo.write(oline+'\n')
                '''
                if len(colors)>=len(gfnames):rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg+'.tts',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
        
    if 'CSS' in sites:
        print('\nprofiling for Coding Start Sites (CSS)')
        d={}
        if len(list(wgs.keys()))>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                if i==0:outmode='w'
                else :outmode='a'
                dic=batchOccAroundPoints(wgs,outname=outname+'_CSS',groupname=gfname+pos_neg+'.css',outmode=outmode,chrColID=1,nameColID=0,posColIDpos=5,posColIDneg=6,straColID=2,sep='\t',second_sep=None,step=args.bin_size,lines=glines,heatMap=args.heatmap,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                if len(colors)>=len(wigaliases):rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg+'.css',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                '''
                fo=open(wfname+pos_neg+'.css.xls','w')
                fo.write('pos\t'+'\t'.join(gfiles)+'\n')
                poses=dic[gfiles[0]].keys()
                poses.sort()
                for i in poses:
                    oline=str(i)
                    for k in gfiles:oline+='\t'+str(dic[k][i])
                    fo.write(oline+'\n')
                '''
                if len(colors)>=len(gfnames):rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg+'.css',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
    
    if 'CTS' in sites:
        print('\nprofiling for Coding Terminal Sites (CTS)')
        d={}
        if len(list(wgs.keys()))>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                if i==0:outmode='w'
                else :outmode='a'
                dic=batchOccAroundPoints(wgs,outname=outname+'_CTS',groupname=gfname+pos_neg+'.cts',outmode=outmode,chrColID=1,nameColID=0,posColIDpos=6,posColIDneg=5,straColID=2,sep='\t',second_sep=None,step=args.bin_size,lines=glines,heatMap=args.heatmap,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                if len(colors)>=len(wigaliases):rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg+'.cts',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                '''
                fo=open(wfname+pos_neg+'.cts.xls','w')
                fo.write('pos\t'+'\t'.join(gfiles)+'\n')
                poses=dic[gfiles[0]].keys()
                poses.sort()
                for i in poses:
                    oline=str(i)
                    for k in gfiles:oline+='\t'+str(dic[k][i])
                    fo.write(oline+'\n')
                '''
                if len(colors)>=len(gfnames):rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg+'.cts',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
    
    if 'ESS' in sites:
        print('\nprofiling for Exon Start Sites (ESS)')
        d={}
        if len(list(wgs.keys()))>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                if i==0:outmode='w'
                else :outmode='a'
                dic=batchOccAroundPoints(wgs,outname=outname+'_ESS',groupname=gfname+pos_neg+'.ess',outmode=outmode,chrColID=1,nameColID=0,posColIDpos=8,posColIDneg=9,straColID=2,sep='\t',second_sep=',',step=args.bin_size,lines=glines,heatMap=args.heatmap,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                if len(colors)>=len(wigaliases):rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg+'.ess',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                '''
                fo=open(wfname+pos_neg+'.ess.xls','w')
                fo.write('pos\t'+'\t'.join(gfiles)+'\n')
                poses=dic[gfiles[0]].keys()
                poses.sort()
                for i in poses:
                    oline=str(i)
                    for k in gfiles:oline+='\t'+str(dic[k][i])
                    fo.write(oline+'\n')
                '''
                if len(colors)>=len(gfnames):rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg+'.ess',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))

    if 'ETS' in sites:
        print('\nprofiling for Exon Terminal Sites (ETS)')
        d={}
        if len(list(wgs.keys()))>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                if i==0:outmode='w'
                else :outmode='a'
                dic=batchOccAroundPoints(wgs,outname=outname+'_ETS',groupname=gfname+pos_neg+'.ets',outmode=outmode,chrColID=1,nameColID=0,posColIDpos=9,posColIDneg=8,straColID=2,sep='\t',second_sep=',',step=args.bin_size,\
                                         lines=glines,heatMap=args.heatmap,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                if len(colors)>=len(wigaliases):rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg+'.ets',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                '''
                fo=open(wfname+pos_neg+'.ets.xls','w')
                fo.write('pos\t'+'\t'.join(gfiles)+'\n')
                poses=dic[gfiles[0]].keys()
                poses.sort()
                for i in poses:
                    oline=str(i)
                    for k in gfiles:oline+='\t'+str(dic[k][i])
                    fo.write(oline+'\n')
                '''
                if len(colors)>=len(gfnames):rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg+'.ets',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
    if 'gene' in sites:
        print('\nprofiling for gene body')
        d={}
        if len(list(wgs.keys()))>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                if i==0:outmode='w'
                else :outmode='a'
                dic=batchOccInRegions(wgs,outname=outname+'_gene',groupname=gfname+pos_neg+'.gene',outmode=outmode,\
                                      chrColID=1,nameColID=0,startColIDpos=3,startColIDneg=4,endColIDpos=4,endColIDneg=3,straColID=2,sep='\t',second_sep=None,step=args.bin_size,\
                                      lines=glines,heatMap=args.heatmap,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP,region_size=args.region_size)
                if len(colors)>=len(wigaliases):rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg+'.gene',region_size=args.region_size,nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                '''
                fo=open(wfname+pos_neg+'.ets.xls','w')
                fo.write('pos\t'+'\t'.join(gfiles)+'\n')
                poses=dic[gfiles[0]].keys()
                poses.sort()
                for i in poses:
                    oline=str(i)
                    for k in gfiles:oline+='\t'+str(dic[k][i])
                    fo.write(oline+'\n')
                '''
                if len(colors)>=len(gfnames):rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg+'.gene',region_size=args.region_size,nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))

    return rcode


def wigProfileToBed(args,wgs,pos_neg='', outname='',outmode='w'):
    rcode=''
    #rcode='par(mfrow=c('+str(nrow)+','+str(ncol)+'))\n'#+rcode
    gfiles=args.bed3file_paths.split(',')
    wigaliases=args.wigfile_aliases.split(',')
    colors=args.plot_colors.split(',')
    if args.bed3file_aliases==None:gfnames=args.bed3file_paths.split(',')
    else:gfnames=args.bed3file_aliases.split(',')
    sites=args.genomic_sites.split(',')
    if 'center' in sites:
        print('\nprofiling for center points in Bed files')
        d={}
        if len(list(wgs.keys()))>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                lcount=len(glines)
                for j in range(lcount):
                    col=glines[j].split()
                    glines[j]='\t'.join([str(j),col[0],str(div((int(col[1])+int(col[2])),2)),'+'])
                    #print glines[j]
                if i==0:outmode='w'
                else:outmode='a'
                dic=batchOccAroundPoints(wgs,outname=outname+'_center',groupname=gfname+pos_neg,outmode=outmode,chrColID=1,nameColID=0,posColIDpos=2,posColIDneg=2,straColID=3,sep='\t',second_sep=None,step=args.bin_size,lines=glines,heatMap=args.heatmap,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                #else :dic=batchOccAroundPoints(wgs,outname=outname,groupname=gfname+pos_neg,outmode='a',chrColID=1,nameColID=0,posColIDpos=2,posColIDneg=2,straColID=3,sep='\t',second_sep=None,step=args.bin_size,lines=glines,heatMap=args.heatmap,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                if len(colors)>=len(wigaliases):rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg,nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                if len(colors)>=len(gfnames):rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg,nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
    if 'region' in sites:
        print('\nprofiling for regions in Bed files')
        d={}
        if len(list(wgs.keys()))>=1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for i in range(len(gfiles)):
                gfile,gfname=gfiles[i],gfnames[i]
                glines=open(gfile).readlines()[1:]
                lcount=len(glines)
                for j in range(lcount):
                    col=glines[j].split()
                    glines[j]='\t'.join([str(j),col[0],col[1],col[2],'+'])
                    #print glines[j]
                if i==0:outmode='w'
                else:outmode='a'
                dic=batchOccInRegions(wgs,outname=outname+'_region',groupname=gfname+pos_neg+'.gene',outmode=outmode,\
                                      chrColID=1,nameColID=0,startColIDpos=2,startColIDneg=3,endColIDpos=3,endColIDneg=2,straColID=4,sep='\t',second_sep=None,step=args.bin_size,\
                                      lines=glines,heatMap=args.heatmap,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP,region_size=args.region_size)
                #if i==0:dic=batchOccAroundPoints(wgs,outname=outname,groupname=gfname+pos_neg,outmode=outmode,chrColID=1,nameColID=0,posColIDpos=2,posColIDneg=2,straColID=3,sep='\t',second_sep=None,step=args.bin_size,lines=glines,heatMap=args.heatmap,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                #else :dic=batchOccAroundPoints(wgs,outname=outname,groupname=gfname+pos_neg,outmode='a',chrColID=1,nameColID=0,posColIDpos=2,posColIDneg=2,straColID=3,sep='\t',second_sep=None,step=args.bin_size,lines=glines,heatMap=args.heatmap,flankup=args.flank_up,flankdn=args.flank_dn,vcal=args.vcal,excludeP=args.excludeP)
                if len(colors)>=len(wigaliases):rcode+=plot(dic=dic,names=wigaliases,outname='',main=gfname+pos_neg,nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
                d[gfname]=dic
        if len(gfiles)>1:
            rcode+='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'
            for wfname in wigaliases:
                dic={}
                for gfname in gfnames:dic[gfname]=d[gfname][wfname]
                if len(colors)>=len(gfnames):rcode+=plot(dic=dic,names=gfnames,outname='',main=wfname+pos_neg,nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab=args.plot_xlab,ylab=args.plot_ylab,colors=args.plot_colors.split(','))
    
    return rcode

def profile(command='profile'):
    '''
    Description:
        This function parses input parameters, calls and passes all parameters values to the main function wigProfileToGene(), or print some help messages if required.
    parameters:
        none  
    '''
    
    if (len(sys.argv)<3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython danpos.py profile <wiggle_file_paths> [optional arguments]\n\nfor more help, please try: python danpos.py profile -h\n")
        return 0
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage="\npython danpos.py <command> <wiggle_file_paths> [optional arguments]\n\n",\
                                     description='',epilog="Kaifu Chen, et al. chenkaifu@gmail.com,  \
                                     Li lab, Biostatistics department, Dan L. Duncan cancer center, \
                                     Baylor College of Medicine.")
    parser.add_argument('command',default=None,\
                        help="set as 'profile' to analyze occupancy profile relative to gene structures or bed format genomic regions.")
    parser.add_argument('wigfile_paths',default=None,\
                        help="Paths to wiggle format data sets. \
                        each path directs to a wiggle format file, \
                        paths must be separated by the comma ','. e.g. file_a.wig,file_b.wig,file_c.wig")
    parser.add_argument('--bed3file_paths',dest='bed3file_paths',metavar='',default=None,\
                        help="Paths to 3-columns bed (bed3) format files, paths must be separated by the comma ','. Only the first 3 columns (chr,start,end) in each file will be used. \
                        Occupancy profile will be calculated around the middle point between each pair of start and end positions. ")
    parser.add_argument('--genefile_paths',dest='genefile_paths',metavar='',default=None,\
                        help="Paths to files each contains a gene sets. \
                        paths must be separated by the comma ','.\
                        Each gene set file must contain at least the following columns ordered as: name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds, \
                        we suggest to download gene set from the UCSC tables at http://genome.ucsc.edu/cgi-bin/hgTables?command=start ")
    parser.add_argument('--genomic_sites',dest='genomic_sites',metavar='',default='TSS,TTS,CSS,CTS,ESS,ETS,gene,center,region',\
                        help="The genomic site names to be analyzed, choose one or several from transcription start site (TSS), \
                        transcription terminal site (TTS), coding start site (CSS), coding terminal site (CTS), exon start site (ESS), \
                        exon terminal site (ETS), region centers in bed file (center), or regions in bed file (region). Names must be separated by the comma ',' ")
    parser.add_argument('--heatmap',dest='heatmap',metavar='',type=int,default=0,\
                        help="Set to 1 to calcute heatmap values, else set to 0 to cancel this analysis." )
    parser.add_argument('--periodicity',dest='periodicity',metavar='',type=int,default=0,\
                        help="Set to 1 to do power spectrum density (PSD) analysis, else set to 0 to cancel this analysis. \
                        This function checks the strength of periodicity in each wiggle data set.")
    parser.add_argument('--wigfile_aliases',dest='wigfile_aliases',metavar='',default=None,\
                        help="A set of short aliases for the wiggle files. \
                        Aliases must be separated by the comma ',' and arranged in the same order as the wiggle files. The data set in each wiggle file \
                        will be represent by the alias in the output files, or by the file name when no alias is provided")
    parser.add_argument('--bed3file_aliases',dest='bed3file_aliases',metavar='',default=None,\
                        help="A set of short aliases for the bed files. \
                        Aliases must be separated by the comma ',' and be arranged in the same order as the bed files. Each bed \
                        file will be represented by its alias in the output files, or by the file name when no alias is provided.")
    parser.add_argument('--genefile_aliases',dest='genefile_aliases',metavar='',default=None,\
                        help="A set of short aliases for the gene files. \
                        Aliases must be separated by the comma ',' and be arranged in the same order as the gene files. Each gene \
                        file will be represented by its alias in the output files, or by the file name when no alias is provided.")
    parser.add_argument('--name', dest="name",metavar='',default='profile',\
                        help="A name for the experiment.")
    parser.add_argument('--vcal', dest="vcal",metavar='',default='mean',\
                        help="The method to calculate plot value at each position relative to a genomic site (e.g. TSS) in a gene group, could be 'median' or 'mean'")
    parser.add_argument('--excludeP',dest="excludeP",metavar='',default=0,type=float,\
                        help="Exclude the extremely outgroup (low or high) values by this percentage.")
    parser.add_argument('--pos_neg',dest="pos_neg",metavar='',default=0,type=int,\
                        help="Set to 0, 1, -1, 2, or 3 if want to do analysis for \
                        positive and negative values together, positive values only, \
                        negative values only, positive and negative values seperately, \
                        or all kinds of analysis (together, positive, negative),\
                        this will be useful when we are analyze the differential signal between two samples.")
    parser.add_argument('--bin_size',dest="bin_size",metavar='',default=10,type=int,\
                        help="Bin size to be used for plotting, it is suggested to be the same as the step or span size in the input wiggle files.")
    parser.add_argument('--flank_up',dest="flank_up",metavar='',default=1500,type=int,\
                        help="How far to calculate from the up stream of each genomic site (e.g., TSS).")
    parser.add_argument('--flank_dn',dest="flank_dn",metavar='',default=1500,type=int,\
                        help="How far to calculate to the down stream of each genomic site (e.g., TSS).")
    parser.add_argument('--region_size',dest="region_size",metavar='',default=1500,type=int,\
                        help="normalize the length of each region (not include the flanking regions) to this size")
    parser.add_argument('--plot_row',dest="plot_row",metavar='',default=2,type=int,\
                        help="Number of rows to plot on each page.")
    parser.add_argument('--plot_column',dest="plot_column",metavar='',default=2,type=int,\
                        help="Number of columns to polt on each page.")
    parser.add_argument('--plot_xmin',dest="plot_xmin",metavar='',default=None,type=float,\
                        help="Minimal value on the x axis of each plot.")
    parser.add_argument('--plot_xmax',dest="plot_xmax",metavar='',default=None,type=float,\
                        help="Maximal value on the x axis of each plot.")
    parser.add_argument('--plot_ymin',dest="plot_ymin",metavar='',default=None,type=float,\
                        help="Minimal value on the y axis of each plot.")
    parser.add_argument('--plot_ymax',dest="plot_ymax",metavar='',default=None,type=float,\
                        help="Maximal value on the y axis of each plot.")
    parser.add_argument('--plot_xlab', dest="plot_xlab",metavar='',default='Relative distance',\
                        help="The label on the x axis.")
    parser.add_argument('--plot_ylab', dest="plot_ylab",metavar='',default='Average signal value',\
                        help="The label on the y axis.")
    parser.add_argument('--plot_colors', dest="plot_colors",metavar='',default='black,gray,red,blue,orange,purple,skyblue,cyan,green,blue4,darkgoldenrod',\
                        help="The colors to be used in the plot.")

    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print('\n')
        parser.print_help()
        print('\n')
        return 0
    elif len(sys.argv)>=4: # at least two parameter need to be specified
        try:
            args=parser.parse_args()  #all paramter values are now saved in args
        except:
            print("\nfor more help, please try: python danpos.py profile -h\n")
            return 0
    else:
        print("\nfor help, please try: python danpos.py profile -h\n")
        return 0 
    
    print('\ncommand:\npython'," ".join(sys.argv)) # print the command line, this let the user to keep a log and remember what parameters they specified
    #if args.genefile_paths==None and args.bed3file_paths==None:
    #    print 'Wrong: both genefile_paths and bed3file_paths are provided, please p'
    #print args
    if args.heatmap==0:args.heatmap=False
    elif args.heatmap==1:args.heatmap=True
    else:
        print('the argument --heatmap can only be 0 or 1')
    args.excludeP=args.excludeP*0.01 #add by kaifu on 12/16/2013
    
    print('\n\nparsing wiggle format data ...')
    if args.wigfile_aliases==None:args.wigfile_aliases=args.wigfile_paths
    wgs=Wigs()
    wigpaths=args.wigfile_paths.split(',')
    wigaliases=args.wigfile_aliases.split(',')
    if len(wigpaths)!=len(wigaliases):
        print('Error: wigfile alias count not equal to path count!\n')
        #parser.print_help()
        print('')
        return
    else:
        for i in range(len(wigpaths)):
            wg=Wig(wigpaths[i],suppress=True)
            if wg.step!=args.bin_size:
                if not wg.changeStep(args.bin_size):return
            wgs.set(wigaliases[i],wg)
            print(wigaliases[i],wgs.data[wigaliases[i]].sum())

    print('\n\nprofiling ...')
    rcode=''
    if args.pos_neg==0 or args.pos_neg==3:
        print('\n\nprofiling  for positive and negative values together...')
        outmode='w'
        if args.genefile_paths!=None:
            rcode+=wigProfileToGene(args,wgs,pos_neg='',outname=args.name,outmode=outmode)
            outmode='a'
        if args.bed3file_paths!=None:rcode+=wigProfileToBed(args,wgs,pos_neg='',outname=args.name,outmode=outmode)
    
    if args.pos_neg==1 or args.pos_neg==2 or args.pos_neg==3:
        print('\n\nprofiling  for positive values ...')
        pwgs=deepcopy(wgs)
        for k in list(pwgs.keys()):pwgs.get(k).rvNeg()
        outmode='w'
        if args.pos_neg==1 or args.pos_neg==2:
            if args.genefile_paths!=None:
                rcode+=wigProfileToGene(args,pwgs,pos_neg='.pos',outname=args.name,outmode=outmode)
                outmode='a'
            if args.bed3file_paths!=None:wigProfileToBed(args,pwgs,pos_neg='.pos',outname=args.name,outmode=outmode)
        else:
            if args.genefile_paths!=None:rcode+=wigProfileToGene(args,pwgs,pos_neg='.pos',outname=args.name,outmode='a')
            if args.bed3file_paths!=None:rcode+=wigProfileToBed(args,pwgs,pos_neg='.pos',outname=args.name,outmode='a')
        
    if args.pos_neg==-1 or args.pos_neg==2 or args.pos_neg==3:
        print('\n\nprofiling  for negative values ...')
        nwgs=deepcopy(wgs)
        for k in list(nwgs.keys()):
            nwgs.get(k).foldChange(-1.0)
            nwgs.get(k).rvNeg()
            nwgs.get(k).foldChange(-1.0)
        if args.pos_neg==2 or args.pos_neg==3:
            if args.genefile_paths!=None:rcode+=wigProfileToGene(args,nwgs,pos_neg='.neg',outname=args.name,outmode='a')
            if args.bed3file_paths!=None:rcode+=wigProfileToBed(args,nwgs,pos_neg='.neg',outname=args.name,outmode='a')
        else:
            outmode='w'
            if args.genefile_paths!=None:
                rcode+=wigProfileToGene(args,nwgs,pos_neg='.neg',outname=args.name,outmode=outmode)
                outmode='a'
            if args.bed3file_paths!=None:rcode+=wigProfileToBed(args,nwgs,pos_neg='.neg',outname=args.name,outmode=outmode)
    
    
    if args.periodicity!=0:
        print('\n\ncalculating periodicity strength...')
        dic=batchOccPSD(wgs,outname=args.name+'.signalPeriodicity')
        rcode+=plot(dic=dic,outname='',main='signalPeriodicity',nrow=args.plot_row,ncol=args.plot_column,xmin=args.plot_xmin,xmax=args.plot_xmax,ymin=args.plot_ymin,ymax=args.plot_ymax,xlab='Periodicity Length (bp)',ylab='Strength',colors=args.plot_colors.split(','))
    
    rcode='pdf("'+args.name+'.pdf")\n'+rcode
    rcode+='dev.off()\n'
    fo=open(args.name+'.R','w')
    fo.write(rcode)
    fo.close()
    r(rcode)
    print('all job done!')
    
def runPositionStatistics(command='runPositionStatistics'):
    '''
    Description:
        This function parses input parameters, calls and passes all parameters values to the main functions related to positions analysis, or print some help messages if required.
    parameters:
        none  
    '''
    
    if (len(sys.argv)<3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython danpos.py stat <file_paths> <file_aliases> [optional arguments]\n\nfor more help, please try: python danpos stat -h\n")
        return 0
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage="\npython danpos.py <command> <file_paths> <file_aliases>[optional arguments]\n\n",\
                                     description='',epilog="Kaifu Chen, et al. chenkaifu@gmail.com,  \
                                     Li lab, Biostatistics department, Dan L. Duncan cancer center, \
                                     Baylor College of Medicine.")
    parser.add_argument('command',default=None,\
                        help="set as 'stat' to run statistical analysis for positions/peaks/regions.")
    parser.add_argument('file_paths',default=None,\
                        help="Paths to the position/peak/region files. Paths must be separated by the comma ',', \
                        each position file should be in the default output format of DANPOS, e.g. the *.positions.xls, *.peaks.xls, or *.regions.xls files. \
                        each file contains positions/peaks/regions from a single sample.\
                        See DANPOS documentation for more information.\
                        e.g. a-b.allpositions.xls,c.positions.xls,d-e.positions.xls")
    parser.add_argument('file_aliases',default=None,\
                        help="A set of aliases for the position/peak/region files. \
                        Alias for different position files must be separated by the comma ',', \
                        and must be arranged in an order corresponding to the file_paths.\
                        e.g. a,b,c")
    parser.add_argument('--name', dest="name",metavar='',default='stat',\
                        help="A name for the experiment.")
    
    parser.add_argument('--plot_row',dest="plot_row",metavar='',default=2,type=int,\
                        help="Number of rows to plot on each page.")
    parser.add_argument('--plot_column',dest="plot_column",metavar='',default=2,type=int,\
                        help="Number of columns to polt on each page.")
    parser.add_argument('--plot_colors', dest="plot_colors",metavar='',default='black,gray,red,blue,orange,purple,skyblue,cyan,green,blue4,darkgoldenrod',\
                        help="The colors to be used in the plot.")

    
    parser.add_argument('-------------------------------        ',dest="separator1",metavar='',default=None,help="")
    parser.add_argument('----- position statistics -----        ',dest="separator1",metavar='',default=None,help="")
    parser.add_argument('-------------------------------         ',dest="separator1",metavar='',default=None,help="")
    
    parser.add_argument('--dis_min',dest="dis_min",metavar='',default=None,type=int,\
                        help="Minimal distance between positions.")
    parser.add_argument('--dis_max',dest="dis_max",metavar='',default=None,type=int,\
                        help="Maximal distance between positions.")
    parser.add_argument('--dis_step',dest="dis_step",metavar='',default=10,type=int,\
                        help="Bin size or step used to calculate the distribution of distance between neighboring positions.")
    
    parser.add_argument('--occ_min',dest="occ_min",metavar='',default=None,type=int,\
                        help="Minimal occupancy value of positions.")
    parser.add_argument('--occ_max',dest="occ_max",metavar='',default=None,type=int,\
                        help="Maximal occupancy value of positions.")
    parser.add_argument('--occ_step',dest="occ_step",metavar='',default=10,type=int,\
                        help="Step or bin size used to calculate occupancy values distribution.")
    
    parser.add_argument('--fuz_min',dest="fuz_min",metavar='',default=None,type=int,\
                        help="Minimal fuzziness score of positions.")
    parser.add_argument('--fuz_max',dest="fuz_max",metavar='',default=None,type=int,\
                        help="Maximal fuzziness score of positions.")
    parser.add_argument('--fuz_step',dest="fuz_step",metavar='',default=None,type=int,\
                        help="Step or bin size used to calculate fuzziness scores distribution.")
    
    parser.add_argument('-------------------------------          ',dest="separator1",metavar='',default=None,help="")
    parser.add_argument('---- peak/region statistics ---          ',dest="separator1",metavar='',default=None,help="")
    parser.add_argument('-------------------------------           ',dest="separator1",metavar='',default=None,help="")
    
    parser.add_argument('--width_min',dest="width_min",metavar='',default=None,type=int,\
                        help="Minimal width of peak or region.")
    parser.add_argument('--width_max',dest="width_max",metavar='',default=None,type=int,\
                        help="Maximal width of peak or region.")
    parser.add_argument('--width_step',dest="width_step",metavar='',default=None,type=int,\
                        help="Step or bin size used to calculate width distribution.")
    
    parser.add_argument('--total_signal_min',dest="total_signal_min",metavar='',default=None,type=int,\
                        help="Minimal total signal value in each peak or region.")
    parser.add_argument('--total_signal_max',dest="total_signal_max",metavar='',default=None,type=int,\
                        help="Maximal total signal value in each peak or region.")
    parser.add_argument('--total_signal_step',dest="total_signal_step",metavar='',default=None,type=int,\
                        help="Step or bin size used to calculate distribution of total signal value in each peak or region.")
    
    parser.add_argument('--height_min',dest="height_min",metavar='',default=None,type=int,\
                        help="Minimal peak or region height.")
    parser.add_argument('--height_max',dest="height_max",metavar='',default=None,type=int,\
                        help="Maximal peak or region height.")
    parser.add_argument('--height_step',dest="height_step",metavar='',default=None,type=int,\
                        help="Step or bin size used to calculate height distribution.")
    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print('\n')
        parser.print_help()
        print('\n')
        return 0
    elif len(sys.argv)>=3: # at least two parameter need to be specified
        try:
            args=parser.parse_args()  #all paramter values are now saved in args
        except:
            print("\nfor more help, please try: python danpos stat -h\n")
            return 0
    else:
        print("\nfor help, please try: python danpos stat -h\n")
        return 0 
    
    print('\ncommand:\npython'," ".join(sys.argv)) # print the command line, this let the user to keep a log and remember what parameters they specified

    occ,dis,fuz,width,auc,height={},{},{},{},{},{}
    positionfiles=args.file_paths.split(',')
    positionaliases=args.file_aliases.split(',')
    if len(positionfiles)!=len(positionaliases):
        print('\nError: position files and aliases counts are not the same!')
        print(len(positionfiles),'position files:',positionfiles)
        print(len(positionaliases),'alias:',positionaliases)
        print('please note that aliases for different files would be separated by the comma \',\'\n')
        return 0
    for i in range(len(positionfiles)):
        positionfile,filealias=positionfiles[i],positionaliases[i]
        title=open(positionfile).readline().split()
        '''
        aliaspair=filealias.split(':')
        if len(aliaspair)>1:
            if 'smt_pos' in title:
                print '\nWrong:\nMore than one alias are provided for',positionfile,':',aliaspair
                print 'Please note that each position file pooled/*positions.xls (generated by DANPOS) would have only one alias.\n'
                return 0
            elif not ( ('fuzziness_diff_FDR' in title) and ('0-log10fuzziness_diff_pval' in title) ):
                print '\nWrong:\nThe input file',positionfile,'is not generated by DANPOS-2.1.1 or later version.'
                print 'If your input files are generated by DANPOS-2.1.0 or earlier version, please use the files *positions.xls under the directory \"pooled/\".\n'
                return 0
            occ[aliaspair[1]]=retrieve_positions_by_value(in_file=positionfile,out_file=None,cr_col_name='chr',pos_col_name='control_smt_loca',val_col_name='control_smt_val',direction_by=[],top_value=None,bottom_value=None)
            occ[aliaspair[0]]=retrieve_positions_by_value(in_file=positionfile,out_file=None,cr_col_name='chr',pos_col_name='treat_smt_loca',val_col_name='treat_smt_val',direction_by=[],top_value=None,bottom_value=None)
            fuz[aliaspair[1]]=retrieve_positions_by_value(in_file=positionfile,out_file=None,cr_col_name='chr',pos_col_name='control_smt_loca',val_col_name='control_fuzziness_score',direction_by=[],top_value=None,bottom_value=None)
            fuz[aliaspair[0]]=retrieve_positions_by_value(in_file=positionfile,out_file=None,cr_col_name='chr',pos_col_name='treat_smt_loca',val_col_name='treat_fuzziness_score',direction_by=[],top_value=None,bottom_value=None)
        else:
            if not 'fuzziness_score' in title:
                print '\nWrong:\nThe input file',positionfile,"doesn't contain the column 'fuzziness_score'."
                print "If your input file is generated by DANPOS-2.1.1, please use the file *allpositions.xls and specify a pair of alias for the control sample and treat sample.\n"
                return 0
        '''
        if 'fuzziness_score' in title:
            occ[filealias]=retrieve_positions_by_value(in_file=positionfile,out_file=None,cr_col_name='chr',pos_col_name='smt_pos',val_col_name='smt_value',direction_by=[],top_value=None,bottom_value=None)
            dis[filealias]=positionDistance(occ[filealias])
            fuz[filealias]=retrieve_positions_by_value(in_file=positionfile,out_file=None,cr_col_name='chr',pos_col_name='smt_pos',val_col_name='fuzziness_score',direction_by=[],top_value=None,bottom_value=None)
        elif 'width_above_cutoff' in title:
            width[filealias]=retrieve_positions_by_value(in_file=positionfile,out_file=None,cr_col_name='chr',pos_col_name='center',val_col_name='width_above_cutoff',direction_by=[],top_value=None,bottom_value=None,log10trans=True)
            auc[filealias]=retrieve_positions_by_value(in_file=positionfile,out_file=None,cr_col_name='chr',pos_col_name='center',val_col_name='total_signal',direction_by=[],top_value=None,bottom_value=None,log10trans=True)
            height[filealias]=retrieve_positions_by_value(in_file=positionfile,out_file=None,cr_col_name='chr',pos_col_name='center',val_col_name='height',direction_by=[],top_value=None,bottom_value=None,log10trans=True)

    print('')
    rcode=''
    if len(occ)>0:
        minmax=positionDicMinMax(dis,lowPercent=0,highPercent=100)
        if args.dis_min==None:args.dis_min=100
        if args.dis_max==None:args.dis_max=250
        if args.dis_step==None:args.dis_step=10
        
        disdic=batchPositionDistanceDistribution(occ,outname=args.name+'.distance',min=args.dis_min,max=args.dis_max,step=args.dis_step)
        rcode=plot(dic=disdic,outname='',main='distance_distribution',nrow=args.plot_row,ncol=args.plot_column,xmin=args.dis_min,xmax=args.dis_max,ymin=None,ymax=None,xlab='Distance',ylab='Percent',colors=args.plot_colors.split(','))
        #rcode+=vioplot(dic=minmax[2],outname='',main='distance_distribution',nrow=args.plot_row,ncol=args.plot_column,ymin=args.dis_min,ymax=args.dis_max)
        
    if len(occ)>0:
        minmax=positionDicMinMax(occ,lowPercent=0,highPercent=99)
        if args.occ_max==None or args.occ_min==None:
            if args.occ_min==None:args.occ_min=minmax[0]
            if args.occ_max==None:args.occ_max=minmax[1]
            if args.occ_step==None:args.occ_step=(args.occ_max-args.occ_min)/100.0
        occdic=batchPositionValDistribution(occ,outname=args.name+'.occupancy',min=args.occ_min,max=args.occ_max,step=args.occ_step)
        rcode+=plot(dic=occdic,outname='',main='occupancy_distribution',nrow=args.plot_row,ncol=args.plot_column,xmin=args.occ_min,xmax=args.occ_max,ymin=None,ymax=None,xlab='Occupancy',ylab='Percent',colors=args.plot_colors.split(','))
        #rcode+=vioplot(dic=minmax[2],outname='',main='occupancy_distribution',nrow=args.plot_row,ncol=args.plot_column,ymin=args.occ_min,ymax=args.occ_max)

    if len(fuz)>0:
        minmax=positionDicMinMax(fuz,lowPercent=1,highPercent=99)
        if args.fuz_max==None or args.fuz_min==None:
            if args.fuz_min==None:args.fuz_min=minmax[0]
            if args.fuz_max==None:args.fuz_max=minmax[1]
            if args.fuz_step==None:args.fuz_step=(args.fuz_max-args.fuz_min)/100.0
        fuzdic=batchPositionValDistribution(fuz,outname=args.name+'.fuzziness',min=args.fuz_min,max=args.fuz_max,step=args.fuz_step)
        rcode+=plot(dic=fuzdic,outname='',main='fuzziness_distribution',nrow=args.plot_row,ncol=args.plot_column,xmin=args.fuz_min,xmax=args.fuz_max,ymin=None,ymax=None,xlab='Fuzziness',ylab='Percent',colors=args.plot_colors.split(','))
        #rcode+=vioplot(dic=minmax[2],outname='',main='fuzziness_distribution',nrow=args.plot_row,ncol=args.plot_column,ymin=args.fuz_min,ymax=args.fuz_max)

    if len(width)>0:
        minmax=positionDicMinMax(width,lowPercent=0,highPercent=100)
        
        if args.width_max==None or args.width_min==None:
            if args.width_min==None:args.width_min=max(minmax[0],1)
            if args.width_max==None:args.width_max=minmax[1]
            if args.width_step==None:args.width_step=(args.width_max-args.width_min)/100.0
        widthdic=batchPositionValDistribution(width,outname=args.name+'.width',min=args.width_min,max=args.width_max,step=args.width_step)
        rcode+=plot(dic=widthdic,outname='',main='log10width distribution',nrow=args.plot_row,ncol=args.plot_column,xmin=args.width_min,xmax=args.width_max,ymin=None,ymax=None,xlab='log10 width',ylab='Percent',colors=args.plot_colors.split(','))
        
        #rcode+=vioplot(dic=minmax[2],outname='',main='log10width distribution',nrow=args.plot_row,ncol=args.plot_column,ymin=args.width_min,ymax=args.width_max)

    if len(auc)>0:
        minmax=positionDicMinMax(auc,lowPercent=0,highPercent=100)
        
        if args.total_signal_max==None or args.total_signal_min==None:
            if args.total_signal_min==None:args.total_signal_min=max(minmax[0],1)
            if args.total_signal_max==None:args.total_signal_max=minmax[1]
            if args.total_signal_step==None:args.total_signal_step=(args.total_signal_max-args.total_signal_min)/100.0
        #print 2,minmax[:2],args.total_signal_max,args.total_signal_min,args.total_signal_step,(args.total_signal_max-args.total_signal_min)/100.0
        aucdic=batchPositionValDistribution(auc,outname=args.name+'.total_signal',min=args.total_signal_min,max=args.total_signal_max,step=args.total_signal_step)
        rcode+=plot(dic=aucdic,outname='',main='log10total_signal distribution',nrow=args.plot_row,ncol=args.plot_column,xmin=args.total_signal_min,xmax=args.total_signal_max,ymin=None,ymax=None,xlab='log10 total signal',ylab='Percent',colors=args.plot_colors.split(','))
        
        #rcode+=vioplot(dic=minmax[2],outname='',main='log10total_signal distribution',nrow=args.plot_row,ncol=args.plot_column,ymin=args.total_signal_min,ymax=args.total_signal_max)

    if len(height)>0:
        minmax=positionDicMinMax(height,lowPercent=0,highPercent=99.8)
        
        if args.height_max==None or args.height_min==None:
            if args.height_min==None:args.height_min=max(minmax[0],1)
            if args.height_max==None:args.height_max=minmax[1]
            if args.height_step==None:args.height_step=(args.height_max-args.height_min)/100.0
        heightdic=batchPositionValDistribution(width,outname=args.name+'.height',min=args.height_min,max=args.height_max,step=args.height_step)
        rcode+=plot(dic=heightdic,outname='',main='log10height distribution',nrow=args.plot_row,ncol=args.plot_column,xmin=args.height_min,xmax=args.height_max,ymin=None,ymax=None,xlab='log10 height',ylab='Percent',colors=args.plot_colors.split(','))
        #rcode+=vioplot(dic=minmax[2],outname='',main='log10height distribution',nrow=args.plot_row,ncol=args.plot_column,ymin=args.height_min,ymax=args.height_max)



    rcode='par(mfrow=c('+str(args.plot_row)+','+str(args.plot_column)+'))\n'+rcode
    rcode='pdf("'+args.name+'.pdf")\n'+rcode
    rcode+='dev.off()\n'
    fo=open(args.name+'.R','w')
    fo.write(rcode)
    fo.close()
    r(rcode)
    print('')
    


def runPositionSelector(command='runPositionSelector'):
    '''
    Description:
        This function parses input parameters, calls and passes all parameters values to the function positionSelector, or print some help messages if required.
    parameters:
        none  
    '''
    
    if (len(sys.argv)<2) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython danpos.py selector  <file_path> [optional arguments]\n\nfor more help, please try: python danpos selector -h\n")
        return 0
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage="\npython danpos.py <command> <file_path>[optional arguments]\n\n",\
                                     description='',epilog="Kaifu Chen, et al. chenkaifu@gmail.com,  \
                                     Li lab, Biostatistics department, Dan L. Duncan cancer center, \
                                     Baylor College of Medicine.")

    parser.add_argument('command',default=None,\
                        help="set as 'selector' to select positions/peaks/regions by value ranges or gene structures.")

    parser.add_argument('file_path',default=None,\
                        help="A path to the position/peak/region file, should be in the default output format of DANPOS (the *.allpositions.xls, *.positions.xls, *.allpeaks.xls, *.peaks.xls, *.allregions.xls, *.regions.xls files). \
                        See DANPOS documentation for more information.")

    parser.add_argument('--valueSelector', dest="valueSelector",metavar='',default=None,\
                        help="Selection requirements defined in the format: 'name1:50:1000,name2::1000,and', \
                        'name1' and 'name2' are column names in the position/peak/region file, \
                        e.g. 'treat_smt_val:50:1000' means selecting treat_smt_val value between 50 than 1000, \
                        'control_smt_val::1000' means selecting control_smt_val value less than 1000, \
                        the 'and' at the end means all selections are required, relace 'and' by 'or' \
                        if only and at least one or selections is required.")
    parser.add_argument('--genicSelector', dest="genicSelector",metavar='',default=None,\
                        help="One or several names selected among transcription start site (TSS), \
                        transcription terminal site (TTS), coding start site (CSS), coding terminal site (CTS), exon start site (ESS), \
                        or exon terminal site (ETS).  Each name\
                        comes with a pair of head and end locations in the format: 'site1:head:end,site2:head:end,and', \
                        e.g. 'TSS:-350:50' means selecting positions/peaks/regions that are in the TSS flanking region \
                        from 350bp upstream to 50bp downstream. The 'and' at the end means all selections are required, relace 'and' by 'or'\
                        if only and at least one selections is required.")
    
    parser.add_argument('--GREATSelector', dest="GREATSelector",metavar='',default=None,\
                        help="Do selection based on regulatory domain defined by the GREAT algorithm around TSS. \
                        E.g. GREAT:-5000:1000:1000000 define a basal regulatory domain from 5kb upstream to 1kb downstream of each TSS, \
                        and further extend each domain up to 1Mb untill reach to basal domains of neighbouring TSSes. For mor information,\
                        see http://bejerano.stanford.edu/great/public/html/")
    parser.add_argument('--out', dest="position_out",metavar='',default='selector.out.xls',\
                        help="An output file name for saving the selected positions/peaks/regions.")

    parser.add_argument('--gene_file', dest="gene_file",metavar='',default=None,\
                        help="A reference gene set, required when need to select nucleosomes by genicSelector. \
                        Gene set file must contain at least the following columns ordered as : name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, \
                        we suggest to download gene file from the UCSC tables at http://genome.ucsc.edu/cgi-bin/hgTables?command=start ")
    
    parser.add_argument('--gene_out', dest="gene_out",metavar='',default=None,\
                        help="A name of the file for saving the genes associated with selected positions/peaks/regions, as defined by genicSelector")
    
    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print('\n')
        parser.print_help()
        print('\n')
        return 0
    
    elif len(sys.argv)>=3: # at least two parameter need to be specified
        try:
            args=parser.parse_args()  #all paramter values are now saved in args
        except:
            print("\nfor more help, please try: python danpos selector -h\n")
            return 0
    else:
        print("\nfor help, please try: python danpos selector -h\n")
        return 0 
    
    print('\ncommand:\npython'," ".join(sys.argv)) # print the command line, this let the user to keep a log and remember what parameters they specified
    if args.genicSelector!=None and args.GREATSelector!=None:
        print("Wrong: please don't use genicSelector and GREATSelector at the same time!")
        return False
    if args.GREATSelector!=None:
        if args.GREATSelector[:6]!='GREAT:':
            print("Wrong: the GREATSelector must starts with 'GREAT:'")
            return False
    
    retr=open(args.file_path).readlines()
    total=max(len(retr)-1,0)
    print(total,'in total')
    if args.valueSelector!=None:retr=positionSelectorByValue(positionLines=retr,selection=args.valueSelector)
    valsel=max(len(retr)-1,0)
    print(valsel,'(',valsel*100.0/max(total,1),'percent) selected by values,')
    if args.genicSelector!=None:
        retr=positionSelectorByGeneStructure(positionLines=retr,selection=args.genicSelector,geneFile=args.gene_file,chrbinsize=1)
        genesel=max(len(retr)-1,0)
        print(genesel,'(',genesel*100.0/max(total,1),'percent) selected by genic structure and value, accounting for',genesel*100.0/max(valsel,1),'percent of value selected.')
    if args.GREATSelector!=None:
        sel=args.GREATSelector[6:]
        retr=positionSelectorByGreatTSS(positionLines=retr,selection=sel,geneFile=args.gene_file,chrbinsize=None)
        genesel=max(len(retr)-1,0)
        print(genesel,'(',genesel*100.0/max(total,1),'percent) selected by the GREAT selector, accounting for',genesel*100.0/max(valsel,1),'percent of value selected.')
    if len(retr)>1:
        fo=open(args.position_out,'w')
        if args.position_out[-3:]=='bed':
            n=0
            for line in retr[1:]:
                n+=1
                col=line.split()
                if col[1]=='-' and col[2]=='-':col[1],col[2]=str(int(col[3])-74),str(int(col[3])+74)
                elif col[1]=='-':col[1]=str(int(col[2])-74)
                elif col[2]=='-':col[2]=str(int(col[1])+74)
                if int(col[1])>int(col[2]):
                    temp=col[1]
                    col[1]=col[2]
                    col[2]=temp
                fo.write('\t'.join(col[:3]+[col[0]+':'+col[1]+'-'+col[2],'0','+'])+'\n')
        else:
            for line in retr:fo.write(line)
        fo.close()
        if args.gene_out!=None:
            gfo=open(args.gene_out,'w')
            gfo.write('gene\tselections\n')
            d={}
            for line in retr[1:]:
                col=line.split()
                tcol=col[-1].split('|')
                for t in tcol:
                    ttcol=t.split(':')
                    tttcol=ttcol[-1].split(',')
                    for ttt in tttcol:
                        ttttcol=ttt.split('/')
                        if ttttcol[0] not in d:d[ttttcol[0]]={}
                        if ttcol[0] not in d[ttttcol[0]]:d[ttttcol[0]][ttcol[0]]=[]
                        d[ttttcol[0]][ttcol[0]].append('/'.join(col[:4]))
            for gene in d:
                line=gene+'\t'
                for site in d[gene]:line+=site+':'+','.join(d[gene][site])+'|'
                gfo.write(line[:-1]+'\n')
            gfo.close()
        print('')
    else:
        print('No positions selected!\n')
    
def retrievePositionValuesAtRanks(command='retrievePositionValuesAtRanks'):
    '''
    Description:
        None
    parameters:
        none  
    '''
    
    if (len(sys.argv)<2) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython danpos.py valuesAtRanks <file_path>  [optional arguments] \n\nfor more help, please try: python danpos valuesAtRanks -h\n")
        return 0
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage="\npython danpos.py <command> <file_path>  [optional arguments] \n\n",\
                                     description='',epilog="Kaifu Chen, et al. chenkaifu@gmail.com,  \
                                     Li lab, Biostatistics department, Dan L. Duncan cancer center, \
                                     Baylor College of Medicine.")

    parser.add_argument('command',default=None,\
                        help="set as 'valuesAtRanks' to retrive position/peak/region values by ranks")

    parser.add_argument('file_path',default=None,\
                        help="a path to the position/peak/region file, should be in the default output format of DANPOS (the *.allpositions.xls, *.positions.xls, *.allpeaks.xls, *.peaks.xls, *.allregions.xls, *.regions.xls files). \
                        See DANPOS documentation for more information.")

    parser.add_argument('--valueRanks', dest="valueRanks",metavar='',default=None,\
                        help="Specify the column names and ranks of values which are to be returned, \
                        e.g. control_smt_val:50:100:-200,treat_smt_val:500:-100 will return the \
                        control_smt_val value ranked at 50 and 100 by increasing order and 200 by decreasing order, as well as\
                        the treat_smt_val value ranked at 500 by increasing order and 100 by decreasing order.\
                        The default is to return value at each quarter for each column name.")

    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print('\n')
        parser.print_help()
        print('\n')
        return 0
    
    elif len(sys.argv)>=3: # at least two parameter need to be specified
        try:
            args=parser.parse_args()  #all paramter values are now saved in args
        except:
            print("\nfor more help, please try: python danpos valuesAtRanks -h\n")
            return 0
    else:
        print("\nfor help, please try: python danpos valuesAtRanks -h\n")
        return 0 
    print('')
    positionLines=open(args.file_path).readlines()
    tcol=positionLines[0].split()
    if args.valueRanks==None:
        #print 'Please specify ranks for one or several of the follow column names:'
        #print '\n'.join(tcol[1:])
        #print ''
        #return
        tlen=len(positionLines)-1
        args.valueRanks=''
        for i in range(4,len(tcol)):args.valueRanks+=tcol[i]+':1:'+str(tlen/4)+':'+str(tlen/2)+':'+str(tlen*3/4)+':'+str(tlen)+','
        args.valueRanks=args.valueRanks[:-1]
    sels=args.valueRanks.split(',')
    for tsel in sels:
        sel=tsel.split(':')
        colid=0
        if not sel[0] in tcol[1:]:
            print("error:", sel[0],'is not a column name in the position file')
            return 0
        else:
            for i in range(len(tcol)):
                if tcol[i]==sel[0]:colid=i
            if colid==0:
                print('error: can not rank by', tcol[colid])
                return 0
        
        print(sel[0])
        values=[]
        for line in positionLines[1:]:
            col=line.split()
            try:values.append(float(col[colid]))
            except:print(col)
        values.sort()
        for r in sel[1:]:
            r=int(r)
            if r>0:print('\trank',r,':',values[r-1])
            else:print('\trank',r,':',values[r])
    print('')
    
if __name__ == "__main__":
    if len(sys.argv)>1:
        if sys.argv[1]=='dpos':runDANPOS(command='dpos')
        elif sys.argv[1]=='dpeak':runDANPOS(command='dpeak')
        elif sys.argv[1]=='dregion':runDANPOS(command='dregion')
        elif sys.argv[1]=='dtriple':runDANPOS(command='dtriple')
        elif sys.argv[1]=='stat':runPositionStatistics()
        elif sys.argv[1]=='profile':profile()
        elif sys.argv[1]=='selector':runPositionSelector()
        elif sys.argv[1]=='valuesAtRanks':retrievePositionValuesAtRanks()
        elif sys.argv[1]=='wiq':wiq()
        elif sys.argv[1]=='wig2wiq':wig2wiq()
        else:printHelp()
    else:
        print('\ndanpos version 2.2.2')
        print('For a list of functions in danpos, please try:\npython danpos.py -h')
        print('')



