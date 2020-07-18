#!/usr/bin/env python
import os,sys,re,argparse
from wig import Wig
from glob import glob
from time import time
from functions import div

def rawsort(ifile,sort_ofile,gfile,format,step=10,suppress=False,buffer=None):
    tm=time()
    if format=='wig':
        print('\nconverting',ifile,'...')
        raw_ofile=sort_ofile[:-3]+'raw.wiq'
        wg=Wig(file=ifile,gfile=gfile,step=step,suppress=suppress)
        wg.ajust_size(gfile=gfile)
        wg.save(file=raw_ofile,format="wiq",step=step,suppress=suppress)
        print('time cost:',time()-tm)
        tm=time()
    else:raw_ofile=ifile
    
    print('\nsorting',raw_ofile,'...')
    temp=ifile[:-3]+'temp'
    while os.path.isdir(temp):temp=temp+'.temp'
    os.mkdir(temp)
    if buffer!=None:cmd='sort -r -n -s -k1 -k2 -o '+sort_ofile+' --buffer-size '+str(buffer)+' --temporary-directory '+str(temp)+' '+raw_ofile
    else:cmd='sort -r -n -s -k1 -k2 -o '+sort_ofile+' --temporary-directory '+str(temp)+' '+raw_ofile
    os.system(cmd)
    if format=='wig':
        print('Removing ',raw_ofile,'...')
        os.system('rm '+raw_ofile)
    print('removing',temp)
    os.system('rm '+str(temp)+' -r')
    print('time cost:',time()-tm)
    
def refquantile(paths,ofile,gfile):
    print('\nPreparing reference ...')
    tm=time()
    files=paths.split(':')
    fi={}
    fo=open(ofile,'w')
    for file in files:fi[file]=open(file)
    nfile=len(files)
    for line in fi[files[0]]:
        col=line.split()
        v=float(col[0])
        for file in files[1:]:
            add_line=fi[file].readline()
            v+=float(add_line.split()[0])
        fo.write(str(div(v,nfile))+'\t-\t-\t-\n')
    fo.close()
    print('time cost:',time()-tm)
    
def changevalue(ifile,ref,ofile,gfile,step=10,suppress=False,buffer=None):
    from random import randint
    if ifile!=ref:print('\nnormalizing',ifile,'...')
    else:print('\nsaving reference ...')
    tm=time()
    fi,fr=open(ifile),open(ref)
    wg=Wig(step=step,gfile=gfile)
    for line in fi:
        col=line.split()
        rcol=fr.readline().split()
        if len(rcol)==0:rcol=[0.0]
        cr,pos,vl=col[2],div(int(col[3]),step),float(rcol[0])
        wg.data[cr][pos]=vl
    n=0
    for line in fr:n+=1
    if n>0:print('Warning: the input genome size is smaller than the reference genome size by',n,'wiggle steps!')
    wg.save(ofile,suppress=suppress)
    print('time cost:',time()-tm)

def qnorwig(paths,gfile,out_dir='wiq_result',ref=None,iformat='wig',rformat='wig',isorted=False,rsorted=False,step=10,suppress=False,buffer=1000000,format_only=False):
    #buffer unit is 1024
    buffer=int(buffer*1000000)
    if isorted==1:isorted=True
    else:isorted=False
    if rsorted==1:rsorted=True
    else:rsorted=False
    alltime=time()
    if not os.path.isdir(out_dir):os.mkdir(out_dir)
    files,sfiles={},{}
    for path in paths.split(':'):
        if os.path.isfile(path):files[path]=1#.append(path)
        elif os.path.isdir(path):
            if iformat=='wig':
                for file in glob(os.path.join(path,'*wig')):files[file]=1
            if iformat=='wiq':
                for file in glob(os.path.join(path,'*wiq')):files[file]=1
    for ifile in files:
        opath=re.sub('/','_',ifile)#'_'.join(elems)
        while opath[0]=='.' or opath[0]=='_' or opath[0]=='/':opath=opath[1:] 
        opath=os.path.join(out_dir,opath)
        if (iformat=='wig') or (isorted==False):
            sfiles[opath[:-3]+'sort.wiq']=opath[:-3]+'qnor.wig'
            rawsort(ifile=ifile,sort_ofile=opath[:-3]+'sort.wiq',gfile=gfile,format=iformat,step=step,suppress=suppress,buffer=buffer)
        else:sfiles[ifile]=opath[:-3]+'qnor.wig'
    if format_only:
        print('job done, total time cost:',time()-alltime,'\n')
        return
    
    if ref==None:
        tnb=0
        for c in paths:tnb+=ord(c)
        refpath=os.path.join(out_dir,'reference.'+str(tnb)+'.sort.wiq')
        refquantile(paths=':'.join(list(sfiles.keys())),ofile=refpath,gfile=gfile)
    else:
        elems=os.path.split(ref)
        opath=os.path.join(out_dir,elems[-1])
        if rformat=='wig' or rsorted==False:
            refpath=opath[:-3]+'sort.wiq'
            rawsort(ifile=ref,sort_ofile=opath[:-3]+'sort.wiq',gfile=gfile,format=rformat,step=step,suppress=suppress,buffer=buffer)
        else:refpath=ref
    for sfile in sfiles:
        changevalue(ifile=sfile,ref=refpath,ofile=sfiles[sfile],gfile=gfile,step=step,suppress=suppress)
    #if ref==None:changevalue(ifile=refpath,ref=refpath,ofile=refpath[:-8]+'qnor.wig',gfile=gfile,step=step,suppress=suppress)
    print('job done, total time cost:',time()-alltime,'\n')
def wig2wiq(command='wiq'):
    '''
    Description:
        This function parses input parameters, calls and passes all parameters values to the function qnorwig, or print some help messages if required.
    parameters:
        none  
    '''
    
    if (len(sys.argv)<2) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython danpos.py wig2wiq <chr_file> <file_path> [optional arguments]\n\nfor more help, please try: python danpos wig2wiq -h\n")
        return 0
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage="\npython danpos.py <command> <chr_file> <file_path>[optional arguments]\n\n",\
                                     description='',epilog="Kaifu Chen, et al. chenkaifu@gmail.com,  \
                                     Li lab, Biostatistics department, Dan L. Duncan cancer center, \
                                     Baylor College of Medicine.")

    parser.add_argument('command',default=None,\
                        help="set as 'wig2wiq' to convert wiggle format to wiq format.")
    
    parser.add_argument('chr_file',default=None,\
                        help="a text file with each line describing a chrosome length, e.g. 'chr1 10043'")
    
    parser.add_argument('input_paths',default=None,\
                        help="paths to the wiggle format input files.\
                        Each path could point to a file or a directory that contains the files. use ':'\
                        to seperate paths, e.g. 'Path_A:path_B")
    parser.add_argument('--out_dir', dest="out_dir",metavar='',default='wig2wiq',\
                        help="a path to the output directory")
    '''
    parser.add_argument('--count', dest="count",metavar='',default=None,type=int,\
                        help="specify the total reads count, e.g. 10000000, so the reads count \
                        in each data set will be normalized to this number")
    '''
    parser.add_argument('--step', dest="step",metavar='',default=10,type=int,\
                        help="the step size in wiggle format data.")
    
    parser.add_argument('--buffer_size', dest="buffer",metavar='',default=10,type=float,\
                        help="maximal memory size that can be used to sort file, e.g. set to 1 when need 1G memory.")
    
    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print('\n')
        parser.print_help()
        print('\n')
        return 0
    elif len(sys.argv)>=3: # at least two parameter need to be specified
        try:
            args=parser.parse_args()  #all paramter values are now saved in args
        except:
            print("\nfor more help, please try: python danpos wig2wiq -h\n")
            return 0
    else:
        print("\nfor help, please try: python danpos wig2wiq -h\n")
        return 0 

    print('\ncommand:\npython'," ".join(sys.argv)) # print the command line, this let the user to keep a log and remember what parameters they specified
    qnorwig(paths=args.input_paths,\
            gfile=args.chr_file,\
            out_dir=args.out_dir,\
            #ref=None,\
            #iformat=args.iformat,\
            #rformat=args.rformat,\
            #isorted=args.isorted,\
            #rsorted=args.rsorted,\
            step=args.step,\
            suppress=False,\
            buffer=args.buffer,\
            format_only=True)

def wiq(command='wiq'):
    '''
    Description:
        This function parses input parameters, calls and passes all parameters values to the function qnorwig, or print some help messages if required.
    parameters:
        none  
    '''
    
    if (len(sys.argv)<2) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least two parameter need to be specified, will print help message if no parameter is specified
        print("\nusage:\npython danpos.py wiq <chr_file> <file_path> [optional arguments]\n\nfor more help, please try: python danpos wiq -h\n")
        return 0
    
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,\
                                     usage="\npython danpos.py <command>  <chr_file> <file_path>[optional arguments]\n\n",\
                                     description='',epilog="Kaifu Chen, et al. chenkaifu@gmail.com,  \
                                     Li lab, Biostatistics department, Dan L. Duncan cancer center, \
                                     Baylor College of Medicine.")

    parser.add_argument('command',default=None,\
                        help="set as 'wiq' to do quantile normalization for wiggle format data.")
    
    parser.add_argument('chr_file',default=None,\
                        help="a text file with each line describing a chrosome length, e.g. 'chr1 10043'")
    
    parser.add_argument('input_paths',default=None,\
                        help="paths to the input files, each file shoule be in '.wig' or '.wiq' format.\
                        Each path could point to a file or a directory that contains the files. use ':'\
                        to seperate paths, e.g. 'Path_A:path_B")
    parser.add_argument('--out_dir', dest="out_dir",metavar='',default='wiq_result',\
                        help="a path to the output directory")
    '''
    parser.add_argument('--count', dest="count",metavar='',default=0,type=int,\
                        help="specify the total reads count, e.g. 10000000, so the reads count \
                        in each data set will be normalized to this number")
    '''
    parser.add_argument('--reference', dest="ref",metavar='',default=None,\
                        help="a path to the file containing reference data. When reference data is available, \
                        data in each input file will be normalized to have the same quantiles in reference data")
    parser.add_argument('--iformat', dest="iformat",metavar='',default='wig',\
                        help="the data format in input files, can be 'wig' or 'wiq'")
    parser.add_argument('--rformat', dest="rformat",metavar='',default='wig',\
                        help="the data format in reference file, can be 'wig' or 'wiq'")
    parser.add_argument('--isorted', dest="isorted",metavar='',default=0,type=int,\
                        help="set to 1 if the input file in wiq format and sorted.")
    parser.add_argument('--rsorted', dest="rsorted",metavar='',default=0,type=int,\
                        help="set to 1 if the reference file is in wiq format and sorted.")
    
    parser.add_argument('--step', dest="step",metavar='',default=10,type=int,\
                        help="the step size in wiggle format data.")
    
    parser.add_argument('--buffer_size', dest="buffer",metavar='',default=10,type=float,\
                        help="maximal memory size that can be used to sort file, e.g. set to 1 when need 1G memory.")
    
    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print('\n')
        parser.print_help()
        print('\n')
        return 0
    elif len(sys.argv)>=3: # at least two parameter need to be specified
        try:
            args=parser.parse_args()  #all paramter values are now saved in args
        except:
            print("\nfor more help, please try: python danpos wiq -h\n")
            return 0
    else:
        print("\nfor help, please try: python danpos wiq -h\n")
        return 0 

    print('\ncommand:\npython'," ".join(sys.argv)) # print the command line, this let the user to keep a log and remember what parameters they specified
    qnorwig(paths=args.input_paths,\
            gfile=args.chr_file,\
            out_dir=args.out_dir,\
            ref=args.ref,\
            iformat=args.iformat,\
            rformat=args.rformat,\
            isorted=args.isorted,\
            rsorted=args.rsorted,\
            step=args.step,\
            suppress=False,\
            buffer=args.buffer)

if __name__ == "__main__":
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # This allows to print each message on screen immediately.

