#!/usr/bin/env python
import os,glob,numpy
from wig import Wig
from rpy2.robjects.packages import importr
from rpy2.robjects import r,FloatVector
from copy import deepcopy
from time import time
import sys
import functions

'''
import sys,os,argparse,glob
'''
class Wigs:
    def __init__(self,path=None,step=0,suppress=False):
        '''
        Parameter:
            file: a pathe to the directory that contain the Wiggle format files
            step: each chrosome will present as an vector, each value in the vector represent a short region of the chrosome, the step value is the size of the short region.
        '''
        self.data = {}  #in the format self.data[name]=Wig class instance
        self.step=step
        if path !=None:
            self.load(path=path,suppress=suppress)
    def ensureSameChrsByRemove(self):
        '''
        Description:
            make sure each Wig instance contain the same set of chrosomes, remove the chrosomes that are not contained by some instance
        Parameter:
            None
        Value:
            None
        '''
        wigs=self.data
        chrs={}
        for wig in wigs:
            for chr in wigs[wig].data:
               if chr in chrs:chrs[chr]+=1
               else:chrs[chr]=1
        crs=list(chrs.keys())
        wnum=len(list(wigs.keys()))
        for cr in crs:
            if chrs[cr]<wnum:chrs.pop(cr)
        for wg in wigs:
            crs=list(wigs[wg].data.keys())
            for cr in crs:
                if cr not in chrs:wigs[wg].pop(cr)
        return wigs
    def foldNormalize(self,scalepairs=None,sampling_total=None,nonzero=False):
        '''
        Description:
            Normalize between Wig class instances by fold change
        Parameter:
            None
        Value:
            None
        '''
        ss=time()
        wigs=self.data
        names=list(wigs.keys())
        names.sort()
        
        if sampling_total==None:
            wsum={}
            for wig in wigs:wsum[wig]=wigs[wig].sum()
            asum=functions.div(sum(wsum.values()),len(list(wsum.keys())))
            for wig in names:
                sys.stdout.write(wig + ' from ' + str(wigs[wig].sum()) + ' to ')
                if scalepairs==None:wigs[wig].foldChange(functions.div(asum*1.0,wsum[wig]))
                else:wigs[wig].foldChange(functions.div(scalepairs[wig]*1.0,wsum[wig]))
                sys.stdout.write(str(wigs[wig].sum()) + "\n")
        else:
            average_total=functions.div(sum(sampling_total.values()),len(list(sampling_total.keys())))
            for name in names:
                sys.stdout.write(name + ' from ' + str(wigs[name].sum()) + ' to ')
                if scalepairs==None:wigs[name].foldChange(functions.div(average_total*1.0,sampling_total[name]))
                else:
                    wigs[name].foldChange(functions.div(scalepairs[name],sampling_total[name]))
                sys.stdout.write(str(wigs[name].sum()) + "\n")
        if nonzero:
            sys.stdout.write('further correction based on count of non-zero base pairs\n')
            gsizes,non0sizes={},{}
            for wig in wigs:
                gsizes[wig]=wigs[wig].gsize()
                non0sizes[wig]=wigs[wig].non0size()
            agsize=functions.div(sum(gsizes.values())*1.0,len(list(gsizes.keys())))
            for wig in wigs:
                sys.stdout.write(wig + ' from ' + str(wigs[wig].sum()) + ' to ')
                wigs[wig].foldChange(functions.div(non0sizes[wig],agsize))
                sys.stdout.write(str(wigs[wig].sum()) +
                    'based on non0size' + non0sizes[wig] + 'and genome size' 
                    + agsize + "\n")
        return 1

    def samplingTotal(self,region_file=None,region_out_file=None,exclude_low_percent=1,exclude_high_percent=1,bnum=100000,nonzero=False):
        '''
        Description:
            caculate the sum of each wig's values after excluding the low and high percentile
        Parameter:
            None
        Value:
            None
        '''
        
        #if exclude_low_percent==0 and exclude_high_percent==0:return None
        #else:
        #print 'calculating normalization factors by sampling ...'
        names=list(self.data.keys())
        if exclude_low_percent==0 and exclude_high_percent==0 and region_file==None: return None
        sampling_total={}
        if region_file==None:
            sys.stdout.write('calculate total signal in each sample after excluding the top ' +
                str(exclude_high_percent) + ' and bottom ' + str(exclude_low_percent) + 
                'percents of genomic regions with extremely high and low signal values\n')
            wsums={}
            for name in names:wsums[name]=self.data[name].sum()
            wavg=functions.div(sum(wsums.values()),len(list(wsums.values())))
            
            rfwig=deepcopy(self.data[names[0]])
            rfwig.foldChange(functions.div(wavg*1.0,wsums[names[0]]))
            for name in names[1:]:
                self.data[name].foldChange(functions.div(wavg*1.0,wsums[name]))
                rfwig.add(self.data[name])
                self.data[name].foldChange(functions.div(wsums[name]*1.0,wavg))
                
            rfwig.foldChange(functions.div(1.0,len(names)))
            lowcut,highcut=rfwig.percentile(p=[exclude_low_percent,100-exclude_high_percent],bnum=bnum,nonzero_end=nonzero)
            rg=rfwig.regionWithinValueRange(lowcut,highcut)
            if region_out_file!=None:rg.save(region_out_file)
        else:
            sys.stdout.write('calculate total signal in each sample in genomic regions defined by' + region_file + "\n")
            rg=Wig(region_file)
        for name in names:sampling_total[name]=self.data[name].multiply(rg).sum()
        sys.stdout.write(str(rg.sum()) + ' (' + str(functions.div(rg.sum()*100.0,rg.gsize())) + 
            ' %) of ' + str(rg.gsize()) + ' base pairs calculated:\n')
        for name in names:
            sys.stdout.write(name + str(sampling_total[name]) + ' (' + 
                str(functions.div(sampling_total[name]*100.0,self.data[name].sum())) + '% of total)\n')

        return sampling_total

    def get(self,k):
        '''
        Description:
            retrieve Wig class instance by name
        Parameter:
            k: the name of the Wig class instance
        Value:
            Wig class instance
        '''
        return self.data[k]
    def keys(self):
        '''
        Description:
            Retrieve the names of all Wig class instances
        Parameter:
            None
        Value:
            a list of names
        '''
        return list(self.data.keys())
    def load(self,path,suppress=False):
        '''
        Description:
            Load multiple Wig class instances from wiggle format files located in one directory
        Parameter:
            path: a path to the directory that contain the wiggle format files
        Value:
            None
        '''
        paths=path
        for path in paths.split(','):
            #wigs={}
            if os.path.isdir(path):
                for infile in glob.glob( os.path.join(path, '*.wig') ): 
                    fname=os.path.split(infile)[-1]
                    if fname[-4:]=='.wig':fname=fname[:-4]
                    self.set(fname,Wig(infile,step=self.step,suppress=suppress)) ########## ---add--- by kaifu on Aug 15,2012 ##########
                    #wigs[infile]=Wig(infile,step=self.step) ########## ---delete--- by kaifu on Aug 15,2012 ##########
            elif os.path.isfile(path):
                fname=os.path.split(path)[-1]
                if fname[-4:]=='.wig':fname=fname[:-4]
                self.set(fname,Wig(path,step=self.step,suppress=suppress)) ########## ---add--- by kaifu on Aug 15,2012 ##########
                #wigs[path]=Wig(path,step=self.step) ########## ---delete--- by kaifu on Aug 15,2012 ##########
            #self.data=wigs
    def nor(self,nor='F',exclude_low_percent=0,exclude_high_percent=0,scalepairs=None,sampling_total=None,nonzero=False):
        '''
        Description:
            normalize among Wig class instances.
        Parameter:
            nor: the normalization method, can be 'Q': quantile normalization, 'F':fold change scaling, 'S': sampling to same coverage, or 'N':no normalization
        Value:
            None
        '''
        if nor!='N' and len(list(self.data.keys()))<2 and scalepairs==None:
            sys.stdout.write("less than 2 datasets, no normalization to be done\n")
            return 1
        #if nor=='Q':
        #sys.stdout.write "quantile normalization"
        #return self.quantileNormalize()
        elif nor=='F':
            #sys.stdout.write "fold normalization"
            return self.foldNormalize(scalepairs=scalepairs,sampling_total=sampling_total,nonzero=nonzero)
        #elif nor=='S':
        #    return self.samplingNormalize()
        elif nor=='N':return 0
        #sys.stdout.write("\nNo normalization method appointed to be done here.\n")
        else: sys.stdout.write("Normalization method "+str(nor)+" not applicable now\n")
        return 0

    def pop(self,k):
        '''
        Description:
            remove Wig class instance by name
        Parameter:
            k: the name of the Wig class instance that is to be removed
        Value:
            None
        '''
        return self.data.pop(k)
    def quantileNormalize(self):
        '''
        Description:
            Normalize between Wig class instances by Quantile
        Parameter:
            None
        Value:
            None
        '''
        ss=time()
        self.ensureSameChrsByRemove()
        wigs=self.data
        r('require("preprocessCore")')
        normq=r('normalize.quantiles')
        
        chrs={}#now it is a dictionary, but will be change to a list later
        for wig in wigs:
            for chr in wigs[wig].data:
                if chr in chrs:chrs[chr]+=1
                else:chrs[chr]=1
        wnum=len(list(wigs.keys()))
        '''
        pops=[]
        for chr in chrs:
            if chrs[chr]<wnum:pops.append(chr)
        for chr in pops:chrs.pop(chr)
        '''
        chrs=list(chrs.keys())#now chrs is a list
        names=list(wigs.keys())
        num=len(names)
        sizes={}
        size=0
        for chr in chrs:
            for name in names:
                if chr not in wigs[name].data:wigs[name].data[chr]=numpy.array([0.0])
                if chr not in sizes:sizes[chr]=wigs[name].data[chr].size
                elif sizes[chr]<wigs[name].data[chr].size:sizes[chr]=wigs[name].data[chr].size
            size+=sizes[chr]
        lst=numpy.array([0.0])
        lst.resize(size*num,refcheck=0)
        
        for i in range(0,num):
            name = names[i]
            tsize=size*i
            for j in range(0,len(chrs)):
                chr = chrs[j]
                wigs[name].data[chr].resize(sizes[chr],refcheck=0)
                ttsize=tsize+sizes[chr]
                lst[tsize:ttsize]+=wigs[name].data[chr][:sizes[chr]]
                tsize=ttsize
        mtr=r.matrix(FloatVector(lst),nrow = size, ncol = num)
        nmtr=normq(mtr)
        for i in range(0,num):
            name=names[i]
            tsize=size*i
            for j in range(0,len(chrs)):
                chr = chrs[j]
                ttsize=tsize+sizes[chr]
                wigs[name].data[chr][:sizes[chr]]=lst[tsize:ttsize]
                tsize=ttsize
        sys.stdout.write('time cost',str(time()-ss)+"\n")
        return 1

    def ajClonal(self,cut=1e-10,extend=1):  ###### add by Kaifu on Nov 14, 2012
        '''
        Description:
            Adjust clonal reads count, fold change between samples will not be altered in this process.
        Parameter:
            cut: the cutoff used to define clonal reads.
                When it is interger,  a read count larger than cut will be defined as clonal;
                when it is float, a read count that is larger than mean count by a Poisson test P value < cut will be defined as clonal.
            fsz: the extension length of each read that is used to calculate the wiggle file. Extension length means the length from 5' end to 3' end,
                e.g. a read may be 36bp when it is generated by the sequencing machine, but it might have been extended to be 80bp or cutted to be 1 bp,
                so the extension length will then be 80bp or 1bp.
        Value:
            None
        Note:
            all wiggle file in a Wigs object must have the same step size.
        '''
        
        #sys.stdout.write '\nremoving clonal singal ...'
        ks=list(self.keys())
        m=deepcopy(self.get(ks[0]))
        if len(ks)>1:
            for k in ks[1:]:m.add(self.get(k))
        m.foldChange(functions.div(1.0,len(ks)))
        avg=m.mean()
        from math import log10,log
        if cut=='0':return
        else:
            if float(cut)>=1:cut=float(cut)
            else:
                co=cut.split('-')
                if len(co)==2:cut=float(co[1])-log10(float(co[0][:-1]))
                else:cut=0-log10(float(co[0]))
                lgpcut=cut
                cut=int(avg+0.5)
                while(0-(functions.div(float(str(functions.ppois(functions.div(cut*1.0,extend),functions.div(avg*1.0,extend),lower_tail=False, log_bool = True)).split()[-1]),log(10)))<lgpcut):cut+=1
        sys.stdout.write('aveage density is ' + str(avg) + ', use clonal signal cutoff ' + str(cut) +  "\n")
        ks=list(self.keys())
        for chr in m.getChrs():
            tchrv=deepcopy(m.data[chr])
            tchrv-=cut#all positive values are count of clonal reads
            tchrv=functions.div(((tchrv**2)**0.5+tchrv),2)#remove all neative values
            tchrv=m.data[chr]-tchrv+numpy.log(tchrv+1)# the addition of '1' is to avoid log(0),"+numpy.log(tchrv+1)" is used to keep the rank order values
            sys.stdout.write(chr+":")
            for k in ks:
                twg=self.get(k)
                if chr not in twg.data:continue
                temp=twg.data[chr].sum()
                sys.stdout.write('\t'+str(k),'reduced from',temp,'to ')
                if chr not in twg.data:twg.data[chr]=numpy.array([0.0])
                if twg.data[chr].size!=tchrv.size:twg.data[chr].resize(tchrv.size,refcheck=0)
                twg.data[chr]=functions.div(tchrv*twg.data[chr],(m.data[chr]+1e-100)) # the addition of '1e-100' is to avoid devide by 0
                sys.stdout.write(twg.data[chr].sum()+ ', percent removed: '+str(100-functions.div(twg.data[chr].sum()*100.0,temp))+"\n")
        
        
    def samplingNormalize(self):
        '''
        Description:
            Normalize between Wig class instances by sampling to same coverage
        Parameter:
            None
        Value:
            None
        '''
        from random import randint
        ss=time()
        wigs=self.data
        wsum={}
        tarray=numpy.array([0.0])
        for wig in wigs:
            wsum[wig]=wigs[wig].sum()
        asum=functions.div(sum(wsum.values()),len(list(wsum.keys())))
        for k in wigs:
            wig=wigs[k]
            oldsum=wig.sum()
            num=oldsum-asum
            if num<0:num=0-num
            else:num=asum
            for chr in wig.data:
                tarray.resize(int(wig.data[chr].sum()),refcheck=0)
                tarray,csz,i,tsz=tarray*0,wig.data[chr].size,0,0
                while i<csz:
                    newtsz=int(tsz+wig.data[chr][i]+0.5)
                    if newtsz>=tarray.size:tarray.resize(newtsz+1000,refcheck=0)
                    while tsz<newtsz:
                        tarray[tsz]=i
                        tsz+=1
                    i+=1
                i,tnum,tsz=0,int(functions.div((functions.div(num*wig.chrSum(chr),oldsum)),wig.step)),tsz
                if oldsum>asum:wig.data[chr]*=0
                while i<tnum:
                    i+=1
                    wig.data[chr][tarray[randint(0,tsz-1)]]+=1
        sys.stdout.write('time cost'+ str(time()-ss) + "\n")
        return 1
    def set(self,k,wig):
        '''
        Description:
            add a Wig class instance
        Parameter:
            k: the name of the new Wig class instance
            wig: the new Wig class instance
        '''
        self.data[k]=wig
    def save(self,path,step=None,format="fixed"):
        '''
        Description:
            save all Wig class instances to wiggle format files in a directory.
            
        Parameter:
            path: the path to the directory that will contain the wiggle format files
            step: the stp size of the wiggle format files
            format: the format of the wiggle files, can be 'fixed' or 'var'
        Value:
            None
        '''
        if step==None:step=self.step
        wigs=self.data
        for k in wigs:
            if path!="":
                if not os.path.isdir(path):os.mkdir(path)
                df=os.path.split(str(k))
                tpath=os.path.join(path,df[-1])
            else:tpath=k
            if tpath[-4:]!='.wig':tpath=tpath+'.wig'
            wigs[k].save(file=tpath,format=format,step=step)
if __name__ == "__main__":
    sys.stdout.write('\n')
    import sys

