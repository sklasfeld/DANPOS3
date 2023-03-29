#!/usr/bin/env python
import numpy,os
from copy import deepcopy
from wig import Wig
from rpy2.robjects import r,FloatVector
from time import time
from random import randint
from math import log
import gzip
import sys
import numpy as np
import functions
import pysam

class reads:
    def __init__(self,file="",rlen=0,step=10,paired=False,cut=1e-10,format='bed'):
        '''
        Description:
            Initialize an instance of the class "reads".
        Parameters:
            rlen: read length
            step: each chrosome will present as an vector, each value in the vector represent a short region of the chrosome, the step value is the size of the short region
            paired: is the reads paired-end
            cut: the cutoff for removing clonal reads, could be P value larger than 0 and small than 1, or count as a positive integer.
            format: the file format of the reads data, can be 'bed','bam',or 'sam'
        Value:
            None
        '''
        self.data={} # an dictionary with chrosome name as key and an numpy array as value, each value in the array represent 1 step in the genome
        self.rlen=rlen #read length
        self.step=step
        if file!="" and format=='bed' and paired:self.loadBedPaired(file=file,step=self.step,cut=cut)
        elif file!="" and format=='bed':self.loadBed(file=file,step=self.step,cut=cut)
        elif file!="" and format=='bam' and paired:self.loadBamPaired(file=file,step=self.step,cut=cut)
        elif file!="" and format=='bam':self.loadBam(file=file,step=self.step,cut=cut)
        elif file!="" and format=='sam' and paired:self.loadSamPaired(file=file,step=self.step,cut=cut)
        elif file!="" and format=='sam':self.loadSam(file=file,step=self.step,cut=cut)
        elif file!="" and format=='bowtie' and paired:self.loadBowtiePaired(file=file,step=self.step,cut=cut)
        elif file!="" and format=='bowtie':self.loadBowtie(file=file,step=self.step,cut=cut)
    def autocorrelation(self,ofile,minsize=10,maxsize=300):
        '''
        Description:
            Calculate correlation between the two strands allowing shift distances.
            
        Parameters:
            ofile: the output file used to write the result to
            minsize: minimal shift distance
            maxsize: maximal shift distance
        
        Value:
            an dictionary with shift distance as key and correlation coefficient  correspondint to the distance
        '''
        step=self.step
        minsize= functions.div(minsize, step)
        maxsize=functions.div(maxsize,step)+1
        print('calculating ...')
        dic={}
        for i in range(int(minsize),int(maxsize)):
            dic[i]=0
            szsum=0
            for chm in self.data:
                sz=self.data[chm]['+'].size-maxsize
                c=self.data[chm]['+'][:sz]*self.data[chm]['-'][i:(sz+i)]
                dic[i]+=sz*functions.div((c.mean()-self.data[chm]['+'][:sz].mean()*self.data[chm]['-'][i:(sz+i)].mean()),(self.data[chm]['+'][:sz].std()*self.data[chm]['-'][i:(sz+i)].std()))
                szsum+=sz
            dic[i]=functions.div(dic[i],(szsum*1.0))
        ks=list(dic.keys())
        ks.sort()
        fo=open(ofile,'w')
        for k in ks:fo.write(str(k)+'\t'+str(dic[k])+'\n')
        return dic

    def clearEmptyEnd(self):
        '''
        Description:
            remove zero values at the end of each chrosome
        
        Parameters:
            None
        
        Value:
            None
        '''
        for chm in self.data:
            final_size=0
            for str in self.data[chm]:
                size=numpy.size(self.data[chm][str])
                if size<=0:continue
                while self.data[chm][str][size-1]==0 and size>0:size-=1#self.data[chm].pop()
                if final_size<size:final_size=size
            for str in self.data[chm]:self.data[chm][str].resize(final_size,refcheck=0)
    def fragSizeDis(self,minsize=10,maxsize=300):
        '''
        Description:
            Calculate most probable size of DNA fragments from which a set of sequencing reads are sequenced.
            
        Parameters:
            minsize: minimal shift distance
            maxsize: maximal shift distance
        
        Value:
            An interger representing the most probable fragment size
            
        '''
        cut=1e-10
        step=self.step
        minsize=functions.div(minsize,step)
        maxsize=functions.div(maxsize,step)+1
        
        avg=self.mean()*self.step ##### '*self.step' is added by Kaifu on Aug 1st, 2012 ##### 
        lgpcut=0-functions.div(log(cut),log(10))
        cut=int(avg+0.5)
        # removed this line below since rpy is not efficient
        while(0-(functions.div(float(str(functions.ppois(cut,avg.item(),lower_tail=False, log_bool = True)).split()[-1]),log(10)))<lgpcut):cut+=1
        if cut<1:cut=1
        
        print('calculating fragment size ...')
        dic={}
        for i in range(int(minsize),int(maxsize)):dic[i]=0
        for chm in self.data:
            sz=self.data[chm]['+'].size-maxsize
            if sz<=0:continue
            tchr=deepcopy(self.data[chm])
            
            for stra in list(tchr.keys()): #remove clonal reads, only necessary in danpos 2.2.0 and later versions
                tchr[stra]-=cut#all positive values are count of clonal reads
                tchr[stra]=functions.div(((tchr[stra]**2)**0.5+tchr[stra]),2)#remove all neative values
                tchr[stra]=self.data[chm][stra]-tchr[stra]
            
            for i in range(int(minsize),int(maxsize)):
                c=tchr['+'][:int(sz)]*tchr['-'][i:(int(sz)+i)]
                dic[i]+=c.sum()
        p=[]
        m=max(dic.values())
        if m<1:m=1
        print('sizes distribution:')
        for i in range(int(minsize),int(maxsize)):
            oline=""
            for j in range(0,int(functions.div(100*dic[i],m)),):oline+='-'
            print(oline,str(i*step)+'bp',str(dic[i]))
            if dic[i]>=m*0.95:p.append(i)
        warning=False
        for i in range(int(minsize),int(minsize)+3):
            if i in p:warning=True
        for i in range(int(maxsize)-3,int(maxsize)):
            if i in p:warning=True
        if warning: print('warning: the probilities of calculated size and up/bottom sizes are too close, we suggest to change up/bottom limit and try again!')
        upv,dpv=0.0,0.0        
        for j in p:
            upv+=j*dic[j]
            dpv+=dic[j]
        finalsize=int(functions.div(upv,dpv)+0.5)
        if (finalsize-minsize)<3 :print('warning: the calculated fragment size seems too close the the bottom limit, we suggest to change bottom limit and try again!')
        if (maxsize-finalsize)<3 :print('warning: the calculated fragment size seems too close the the up limit, we to change up limit and try again!')
        print('potential size:', end=' ')
        for t in p:print(t*step, end=' ')
        print('')
        print('calculated fragment size:',finalsize*step)
        return finalsize*step

    def loadBowtie(self,file="",step=10,cut=1e-10):
        '''
        Description:
            load single-end reads data from a file in '.bed' format
        
        parameter:
            file: a path to the file containing the sequencing reads
            step: each chrosome will present as an vector, each value in the vector represent a short region of the chrosome, the step value is the size of the short region
            cut: the cutoff for removing clonal reads, could be P value larger than 0 and small than 1, or count as a positive integer.
        Value:
            None
        '''
        print('\nparsing from bowtie file',file,'...')

        oldchr=""
        sizes={}
        num=0
        infile=open(file)
        if file[-2:]=='gz':infile=gzip.open(file)
        for line in infile:
            try:
                col=line[:-1].split('\t')
                chm,start,stra=col[2],int(col[3]),col[1]
                end=start+len(col[4])
                #chr,start,end,name,score,stra=col[0:6]
                #start,end=int(start),int(end)
            except:
                print(line.split())#'wrong line:',line[:-1]
                continue
            num+=1
            '''
            if loadcount>0:
                if num>loadcount:
                    self.clearEmptyEnd()
                    for chm in self.data:
                        lth=max(self.data[chm]['+'].size,self.data[chm]['-'].size)
                        self.data[chm]['+'].resize(lth,refcheck=0)
                        self.data[chm]['-'].resize(lth,refcheck=0)
                    if cut>0:self.rvClonal(cut=cut)
                    print 'parsing finished,',num-1,'reads parsed'
                    return
            '''

            if num%1000000==0: print(num,'reads parsed')
            if stra=='+':mid=functions.div(start,step)
            elif stra=='-':mid=functions.div(end,step)
            if chm not in self.data:
                self.data[chm]={'+':numpy.array([0.0]),'-':numpy.array([0.0])}
                sizes[chm]=0
            if mid>=sizes[chm]:
                sizes[chm]=mid+1000
                self.data[chm]['+'].resize(sizes[chm],refcheck=0)
                self.data[chm]['-'].resize(sizes[chm],refcheck=0)
            if mid>=0:self.data[chm][stra][mid]+=1.0
        self.clearEmptyEnd()
        for chm in self.data:
            lth=max(self.data[chm]['+'].size,self.data[chm]['-'].size)
            self.data[chm]['+'].resize(lth,refcheck=0)
            self.data[chm]['-'].resize(lth,refcheck=0)
        if cut>0:self.rvClonal(cut=cut)
        print('parsing finished,',num,'reads parsed')
    
    def loadBowtiePaired(self,file="",step=10,cut=1e-10):
        '''
        Description:
            load paired-end reads data from a file in '.bed' format, the "name" fields of each reads pair must be the same except for the last one character in each field, which should be either '1' or '2', e.g. 'reads/1' or 'reads/2', each pair of reads must be arranged in two neighboring lines.
        
        parameter:
            file: a path to the file containing the sequencing reads
            step: each chrosome will present as an vector, each value in the vector represent a short region of the chrosome, the step value is the size of the short region
            cut: the cutoff for removing clonal reads, could be P value larger than 0 and small than 1, or count as a positive integer.
        Value:
            None
        '''
        print('\nparsing from bowtie file',file,'...')

        oldchr=""
        sizes={}
        fragsizes={}
        num=0
        serr,nerr=0,0
        infile=open(file)
        if file[-2:]=='gz':infile=gzip.open(file)
        end1,end2=[],[]
        for line in infile:
            try:
                col=line[:-1].split('\t')
                chm,start,name,score,stra=col[2],int(col[3]),col[0],0,col[1]
                tnames=name.split()
                if len(tnames)==2:name=tnames[0]+name[-1] ########## add by kaifu on sep 5, 2012 ########## some time the reads name have two words seperated by a space (not '\t')
                end=start+len(col[4])
                if len(end1)<1:
                    end1=[chm,start,end,name,score,stra]#col
                    num+=1
                    if num%1000000==0: print(num,'reads parsed')
                    continue
                else:
                    end2=[chm,start,end,name,score,stra]#col
                    num+=1
                    if num%1000000==0: print(num,'reads parsed')
            except:
                print('wrong line:',line[:-1])
                continue
            if end1[3][:-1]==end2[3][:-1]:
                if end1[5]=='+' and end2[5]=='-':
                    chm = end1[0]
                    mid = int(functions.div((end1[1]+end2[2]),(2*step)))
                    fragsize = end2[2]-end1[1]
                elif end1[5]=='-' and end2[5]=='+':
                    chm = end1[0]
                    mid = int(functions.div((end1[2]+end2[1]),(2*step)))
                    fragsize = end2[1]-end1[2]
                else:
                    print('pair error --- reads from same strand:\n',end1,'\n',end2,'\n')
                    serr+=1
                    end1,end2=[],[]
                    continue
                if fragsize not in fragsizes:fragsizes[fragsize]=1
                else:fragsizes[fragsize]+=1
            else:
                print('pair error --- single end reads:\n',end1,'\n')
                nerr+=1
                end1=end2
                end2=[]
                continue
             
            if chm not in self.data:
                self.data[chm]={'+':numpy.array([0.0]),'-':numpy.array([0.0])}
                sizes[chm]=0
            if mid>=sizes[chm]:
                sizes[chm]=mid+1000
                self.data[chm]['+'].resize(sizes[chm],refcheck=0)
                self.data[chm]['-'].resize(sizes[chm],refcheck=0)
            if mid>0:
                self.data[chm]['+'][mid]+=1.0
                self.data[chm]['-'][mid]+=1.0
            end1,end2=[],[]
        self.clearEmptyEnd()
        for chm in self.data:
            lth=max(self.data[chm]['+'].size,self.data[chm]['-'].size)
            self.data[chm]['+'].resize(lth,refcheck=0)
            self.data[chm]['-'].resize(lth,refcheck=0)
        print('parsing finished,',num,'reads parsed')
        maxv,lths=max(fragsizes.values()),list(fragsizes.keys())
        lths.sort()
        tlth,count=0,0
        maxlth=[]
        for lth in lths:
            if fragsizes[lth]==maxv:maxlth.append(lth)
            tlth+=lth*fragsizes[lth]
            count+=fragsizes[lth]
            dcount=np.ceil(functions.div(100*fragsizes[lth],maxv))
            if dcount>4: print('-'*int(dcount),lth,fragsizes[lth])
        print('average fragment size:',functions.div(tlth*1.0,count))
        print('most enriched fragment size:',maxlth)
        print(serr,'pairs failed due to locations on same strands')
        print(nerr,'reads have no mate reads')
        if cut>0:self.rvClonal(cut=cut)


    def loadBed(self,file="",step=10,cut=1e-10):#,loadcount=0):
        '''
        Description:
            load single-end reads data from a file in '.bed' format
        
        parameter:
            file: a path to the file containing the sequencing reads
            step: each chrosome will present as an vector, each value in the vector represent a short region of the chrosome, the step value is the size of the short region
            cut: the cutoff for removing clonal reads, could be P value larger than 0 and small than 1, or count as a positive integer.
        Value:
            None
        '''
        print('\nparsing from bed file',file,'...')

        oldchr=""
        sizes={}
        num=0
        infile=open(file)
        if file[-2:]=='gz':infile=gzip.open(file)
        for line in infile:
            try:
                col=line.split()
                #if len(col)==5:chm,start,end,name,stra=col[0:5]
                #else:
                chm,start,end,name,score,stra=col[0:6]
                start,end=int(start),int(end)
            except:
                print('wrong format:',line.split())#'wrong line:',line[:-1]
                continue
            #if col[0]!='chr1':continue ################################################ just for test ################################################
            num+=1
            '''
            if loadcount>0:
                if num>loadcount:
                    self.clearEmptyEnd()
                    for chm in self.data:
                        lth=max(self.data[chm]['+'].size,self.data[chm]['-'].size)
                        self.data[chm]['+'].resize(lth,refcheck=0)
                        self.data[chm]['-'].resize(lth,refcheck=0)
                    if cut>0:self.rvClonal(cut=cut)
                    print 'parsing finished,',num-1,'reads parsed'
                    return
            '''
            if num%1000000==0: print(num,'reads parsed')
            if stra=='+':mid=int(functions.div(start,step))
            elif stra=='-':mid=int(functions.div(end,step))
            if chm not in self.data:
                self.data[chm]={'+':numpy.array([0.0]),'-':numpy.array([0.0])}
                sizes[chm]=0
            if mid>=sizes[chm]:
                sizes[chm]=mid+1000
                self.data[chm]['+'].resize(sizes[chm],refcheck=0)
                self.data[chm]['-'].resize(sizes[chm],refcheck=0)
            if mid>=0:self.data[chm][stra][mid]+=1.0
        self.clearEmptyEnd()
        for chm in self.data:
            lth=max(self.data[chm]['+'].size,self.data[chm]['-'].size)
            self.data[chm]['+'].resize(lth,refcheck=0)
            self.data[chm]['-'].resize(lth,refcheck=0)
        if cut>0:self.rvClonal(cut=cut)
        print('parsing finished,',num,'reads parsed')
    def loadBedPaired(self,file="",step=10,cut=1e-10):
        '''
        Description:
            load paired-end reads data from a file in '.bed' format, the "name" 
            fields of each reads pair must be the same except for the last one 
            character in each field, which should be either '1' or '2', e.g. 
            'reads/1' or 'reads/2', each pair of reads must be arranged in two 
            neighboring lines.
        
        parameter:
            file: a path to the file containing the sequencing reads
            step: each chrosome will present as an vector, each value in the vector represent a short region of the chrosome, the step value is the size of the short region
            cut: the cutoff for removing clonal reads, could be P value larger than 0 and small than 1, or count as a positive integer.
        Value:
            None
        '''
        print('\nparsing from bed file',file,'...')

        oldchr=""
        sizes={}
        fragsizes={}
        num=0
        serr,nerr=0,0
        infile=open(file)
        if file[-2:]=='gz':infile=gzip.open(file)
        end1,end2=[],[]
        for line in infile:
            try:
                col=line[:-1].split('\t')
                tnames=col[3].split()
                if len(tnames)==2:col[3]=tnames[0]+col[3][-1] ########## add by kaifu on sep 5, 2012 ########## some time the reads name have two words seperated by a space (not '\t')
                col[1],col[2]=int(col[1]),int(col[2])
                if len(end1)<1:
                    end1=col
                    num+=1
                    if num%1000000==0: print(num,'reads parsed')
                    continue
                else:
                    end2=col
                    num+=1
                    if num%1000000==0: print(num,'reads parsed')
            except:
                print('wrong line:',line[:-1])
                continue
            if end1[3][:-1]==end2[3][:-1]:
                if end1[5]=='+' and end2[5]=='-':
                    chm,mid,fragsize=end1[0],int(functions.div((end1[1]+end2[2]),(2*step))),end2[2]-end1[1]
                elif end1[5]=='-' and end2[5]=='+':
                    chm,mid,fragsize=end1[0],int(functions.div((end1[2]+end2[1]),(2*step))),end2[1]-end1[2]
                else:
                    #print 'pair error --- reads from same strand:\n',end1,'\n',end2,'\n'
                    serr+=1
                    end1,end2=[],[]
                    continue
                if fragsize not in fragsizes:fragsizes[fragsize]=1
                else:fragsizes[fragsize]+=1
            else:
                #print 'pair error --- single end reads:\n',end1,'\n'
                nerr+=1
                end1=end2
                end2=[]
                continue
             
            if chm not in self.data:
                self.data[chm]={'+':numpy.array([0.0]),'-':numpy.array([0.0])}
                sizes[chm]=0
            if mid>=sizes[chm]:
                sizes[chm]=mid+1000
                self.data[chm]['+'].resize(sizes[chm],refcheck=0)
                self.data[chm]['-'].resize(sizes[chm],refcheck=0)
            if mid>0:
                self.data[chm]['+'][mid]+=1.0
                self.data[chm]['-'][mid]+=1.0
            end1,end2=[],[]
        self.clearEmptyEnd()
        for chm in self.data:
            lth=max(self.data[chm]['+'].size,self.data[chm]['-'].size)
            self.data[chm]['+'].resize(lth,refcheck=0)
            self.data[chm]['-'].resize(lth,refcheck=0)
        print('parsing finished,',num,'reads parsed')
        maxv,lths=max(fragsizes.values()),list(fragsizes.keys())
        lths.sort()
        tlth,count=0,0
        maxlth=[]
        for lth in lths:
            if fragsizes[lth]==maxv:maxlth.append(lth)
            tlth+=lth*fragsizes[lth]
            count+=fragsizes[lth]
            dcount=np.ceil(functions.div(100*fragsizes[lth],maxv))
            print(dcount)
            if dcount>4: print('-'*int(dcount),lth,fragsizes[lth])
        print('average fragment size:',functions.div(tlth*1.0,count))
        print('most enriched fragment size:',maxlth)
        print(serr,'pairs failed due to locations on same strands')
        print(nerr,'reads have no mate reads')
        if cut>0:self.rvClonal(cut=cut)

    def loadSam(self,file="",step=10,cut=1e-10):
        '''
        Description:
            load single-end reads data from a file in '.sam' format
        
        parameter:
            file: a path to the file containing the sequencing reads
            step: each chrosome will present as an vector, each value in the vector represent a short region of the chrosome, the step value is the size of the short region
            cut: the cutoff for removing clonal reads, could be P value larger than 0 and small than 1, or count as a positive integer.
        Value:
            None
        '''
        print('\nparsing from sam file',file,'...')

        oldchr=""
        sizes={}
        num=0
        sam_reader = pysam.AlignmentFile(file, "r")
        for sam in sam_reader.fetch():
            try:
                #col=line.split('\t')
                if sam.is_unmapped:continue
                name = sam.query_name
                chm = sam.reference_name
                stra = '+'
                start=int(sam.reference_start + 1)
                rlen=sam.reference_length
                if sam.is_reverse:stra='-'
                end,score=start+rlen,'1'
            except:
                print(sam)#'wrong line:',line[:-1]
                continue
            num+=1
            if num%1000000==0: print(num,'reads parsed')
            if stra=='+':mid=int(functions.div(start,step))
            elif stra=='-':mid=int(functions.div(end,step))
            if chm not in self.data:
                self.data[chm]={'+':numpy.array([0.0]),'-':numpy.array([0.0])}
                sizes[chm]=0
            if mid>=sizes[chm]:
                sizes[chm]=mid+1000
                self.data[chm]['+'].resize(sizes[chm],refcheck=0)
                self.data[chm]['-'].resize(sizes[chm],refcheck=0)
            if mid>=0:self.data[chm][stra][mid]+=1.0
        if num==0: sys.exit("ERROR: Empty Sam File")
        self.clearEmptyEnd()
        for chm in self.data:
            lth=max(self.data[chm]['+'].size,self.data[chm]['-'].size)
            self.data[chm]['+'].resize(lth,refcheck=0)
            self.data[chm]['-'].resize(lth,refcheck=0)
        if cut>0:self.rvClonal(cut=cut)
        print('parsing finished,',num,'reads parsed')
    def loadSamPaired(self,file="",step=10,cut=1e-10):
        '''
        Description:
            load paired-end reads data from a file in '.sam' format.
        
        parameter:
            file: a path to the file containing the sequencing reads
            step: each chrosome will present as an vector, each value in the vector represent a short region of the chrosome, the step value is the size of the short region
            cut: the cutoff for removing clonal reads, could be P value larger than 0 and small than 1, or count as a positive integer.
        Value:
            None
        '''
        oldchr=""
        sizes={}
        fragsizes={}
        num=0
        print('\nparsing from sam file',file,'...')
        sam_reader = pysam.AlignmentFile(file, "r")
        for sam in sam_reader.fetch():
            try:
                if not sam.is_proper_pair:continue
                ref_start=int(sam.reference_start + 1)
                next_ref_start=int(sam.next_reference_start + 1)
                rlen=sam.reference_length
                name = sam.query_name
                chm = sam.reference_name
                stra = '+'
                score = '1'
                mid=int(np.floor(functions.div((ref_start+next_ref_start+rlen),(2*step))))
                if sam.is_reverse:stra='-'
                if next_ref_start>ref_start:fragsize=next_ref_start-ref_start+rlen
                else:fragsize=ref_start-next_ref_start+rlen
                if fragsize not in fragsizes:fragsizes[fragsize]=1
                else:fragsizes[fragsize]+=1
            except:
                print(sam)#'wrong line:',line[:-1]
                continue
            num+=1
            if num%1000000==0: print(num,'reads parsed')
            if chm not in self.data:
                self.data[chm]={'+':numpy.array([0.0]),'-':numpy.array([0.0])}
                sizes[chm]=0
            if mid>=sizes[chm]:
                sizes[chm]=mid+1000
                self.data[chm]['+'].resize(sizes[chm],refcheck=0)
                self.data[chm]['-'].resize(sizes[chm],refcheck=0)
            if mid>=0:self.data[chm][stra][mid]+=1.0
        if num==0: sys.exit("ERROR: Empty Sam File")
        self.clearEmptyEnd()
        for chm in self.data:
            lth=max(self.data[chm]['+'].size,self.data[chm]['-'].size)
            self.data[chm]['+'].resize(lth,refcheck=0)
            self.data[chm]['-'].resize(lth,refcheck=0)
        if cut>0:self.rvClonal(cut=cut)
        print('parsing finished,',num,'reads parsed')
        maxv,lths=max(fragsizes.values()),list(fragsizes.keys())
        lths.sort()
        tlth,count=0,0
        maxlth=[]
        for lth in lths:
            if fragsizes[lth]==maxv:maxlth.append(lth)
            tlth+=lth*fragsizes[lth]
            count+=fragsizes[lth]
            dcount=np.ceil(functions.div(100*fragsizes[lth],maxv))
            if dcount>4: print('-'*int(dcount),lth,fragsizes[lth])
        print('average fragment size:',functions.div(tlth*1.0,count))
        print('most enriched fragment size:',maxlth)

    def loadBam(self,file="",step=10,cut=1e-10):
        '''
        Description:
            load single-end reads data from a file in '.bam' format
        
        parameter:
            file: a path to the file containing the sequencing reads
            step: each chrosome will present as an vector, each value in the vector represent a short region of the chrosome, the step value is the size of the short region
            cut: the cutoff for removing clonal reads, could be P value larger than 0 and small than 1, or count as a positive integer.
        Value:
            None
        '''
        print('\nparsing from Bam file',file,'...')

        oldchr=""
        sizes={}
        num=0
        sam_reader = pysam.AlignmentFile(file, "rb")
        for sam in sam_reader.fetch():
            try:
                # skip unmapped reads
                if sam.is_unmapped:continue
                name = sam.query_name
                chm = sam.reference_name
                stra = '+'
                start=int(sam.reference_start + 1)
                rlen=sam.reference_length
                # set `stra` to "-" for reverse reads
                if sam.is_reverse:stra='-'
                end,score=start+rlen,'1'
            except:
                print(sam)#'wrong line:',line[:-1]
                continue
            num+=1
            if num%1000000==0: print(num,'reads parsed')
            if stra=='+':mid=int(functions.div(start,step))
            elif stra=='-':mid=int(functions.div(end,step))
            if chm not in self.data:
                self.data[chm]={'+':numpy.array([0.0]),'-':numpy.array([0.0])}
                sizes[chm]=0
            if mid>=sizes[chm]:
                sizes[chm]=mid+1000
                self.data[chm]['+'].resize(sizes[chm],refcheck=0)
                self.data[chm]['-'].resize(sizes[chm],refcheck=0)
            if mid>=0:self.data[chm][stra][mid]+=1.0
        if num==0: sys.exit("ERROR: Empty Bam File")
        self.clearEmptyEnd()
        for chm in self.data:
            lth=max(self.data[chm]['+'].size,self.data[chm]['-'].size)
            self.data[chm]['+'].resize(lth,refcheck=0)
            self.data[chm]['-'].resize(lth,refcheck=0)
        if cut>0:self.rvClonal(cut=cut)
        print('parsing finished,',num,'reads parsed')
    def loadBamPaired(self,file="",step=10,cut=1e-10):
        '''
        Description:
            load paired-end reads data from a file in '.bam' format
        
        parameter:
            file: a path to the file containing the sequencing reads
            step: each chrosome will present as an vector, each value in the vector represent a short region of the chrosome, the step value is the size of the short region
            cut: the cutoff for removing clonal reads, could be P value larger than 0 and small than 1, or count as a positive integer.
        Value:
            None
        '''
        print('\nparsing from bam file',file,'...')

        oldchr=""
        sizes={}
        fragsizes={}
        num=0
        sam_reader = pysam.AlignmentFile(file, "rb")
        for sam in sam_reader.fetch():
            try:
                if not sam.is_proper_pair:continue
                ref_start=int(sam.reference_start + 1)
                next_ref_start=int(sam.next_reference_start + 1)
                rlen=sam.reference_length
                name = sam.query_name
                chm = sam.reference_name
                stra = '+'
                score = '1'
                mid=int(functions.div((ref_start+next_ref_start+rlen),(2*step)))
                if sam.is_reverse:stra='-'
                if next_ref_start>ref_start:fragsize=next_ref_start-ref_start+rlen
                else:fragsize=ref_start-next_ref_start+rlen
                if fragsize not in fragsizes:fragsizes[fragsize]=1
                else:fragsizes[fragsize]+=1
            except:
                print(sam)#'wrong line:',line[:-1]
                sys.exit("SOME REASON DID NOT WORK...")
                continue
            num+=1
            if num%1000000==0: print(num,'reads parsed')
            if chm not in self.data:
                self.data[chm]={'+':numpy.array([0.0]),'-':numpy.array([0.0])}
                sizes[chm]=0
            if mid>=sizes[chm]:
                sizes[chm]=mid+1000
                self.data[chm]['+'].resize(sizes[chm],refcheck=0)
                self.data[chm]['-'].resize(sizes[chm],refcheck=0)
            if mid>=0:self.data[chm][stra][mid]+=1.0
        if num==0: sys.exit("ERROR: Empty Bam File")
        self.clearEmptyEnd()
        for chm in self.data:
            lth=max(self.data[chm]['+'].size,self.data[chm]['-'].size)
            self.data[chm]['+'].resize(lth,refcheck=0)
            self.data[chm]['-'].resize(lth,refcheck=0)
        if float(cut)>0:self.rvClonal(cut=cut)
        print('parsing finished,',num,'reads parsed')
        maxv,lths=max(fragsizes.values()),list(fragsizes.keys())
        lths.sort()
        tlth,count=0,0
        maxlth=[]
        for lth in lths:
            if fragsizes[lth]==maxv:maxlth.append(lth)
            tlth+=lth*fragsizes[lth]
            count+=fragsizes[lth]
            dcount=int(np.ceil(functions.div(100*fragsizes[lth],maxv)))
            if dcount>4: print('-'*int(dcount),lth,fragsizes[lth])
        print('average fragment size:',functions.div(tlth*1.0,count))
        print('most enriched fragment size:',maxlth)

    def mean(self):
        '''
        Description:
            calculate the average reads count per nucleotide
        
        Parameter:
            None
        
        Value:
            None
        '''
        
        return functions.div(self.sum()*1.0,self.size())

    def rvClonal(self,cut=1e-10):
        '''
        Description:
            Remove clonal reads, not that this process may change the difference between samples.
            If dno't want to change the difference, transfer all samples to wigs and use the rvClonal function of the Wigs class.
            
        Parameter:
            cut: the cutoff for removing clonal reads, could be P value larger than 0 and small than 1, or count as a positive integer.
        
        Value:
            None
        '''
        ss=time()
        cut=float(cut)
        if cut<=0:return
        print('removing clonal reads ...')
        avg=functions.div(self.mean(),self.step) ##### '*self.step' is added by Kaifu on Aug 1st, 2012 ##### and is further changed to be /self.step on May29, 2014
        if  cut>0 and cut<1:
            lgpcut=0-functions.div(log(cut),log(10))
            cut=float(avg+0.5)
            # removed this line below since rpy is not efficient
            pois_dist = functions.ppois(cut,float(avg), lower_tail=False, log_bool = True)
            while(0-(functions.div(pois_dist,log(10)))<lgpcut):
                cut+=1
                pois_dist = functions.ppois(cut,float(avg), lower_tail=False, log_bool = True)
            if cut<1:cut=1
        print('whole genome average reads density is',avg,'use cutoff:',cut)
        cut*=self.step ##### added by Kaifu on May29, 2014
        cnum=0
        rreads=0
        treads=0
        tchrv=numpy.array([0.0])
        before=self.sum()
        print('before removing:',before,'reads')
        for chm in self.data:
            for stra in list(self.data[chm].keys()):
                tchrv=deepcopy(self.data[chm][stra])
                tchrv-=cut#all positive values are count of clonal reads
                tchrv=functions.div(((tchrv**2)**0.5+tchrv),2)#remove all neative values
                self.data[chm][stra]-=tchrv
        after=self.sum()
        print('after removing:',after,'reads')
        print(functions.div((before-after)*100,before),'percent removed.')
        print('time cost:',time()-ss)

    def size(self):
        '''
        Description:
            Calculate the total genome size.
        
        Parameter:
            None
        
        value:
            Interger value representing genome size.
        '''
        lth=0
        for chm in self.data:
            for stra in self.data[chm]:
                lth+=self.data[chm][stra].size
        return lth*self.step ##### '*self.step' is added by Kaifu on Aug 1st,2012 #####
    def sizeAdjust(self,gfile):
        '''
        Description:
            Adjust the size of each chrosome.
        
        Parameter:
            gfile: path to the file containing the size of each chrosome, each line in the file would be in the format "chrosome_name size", in which size is an integer value, and chrosome_name should contain no empty space
        
        value:
            None.
        '''
        sizes={}
        for line in open(gfile):
            col=line.split()
            sizes[col[0]]=functions.div(int(col[1]),self.step)
        for chm in self.data:
            if chm not in sizes:self.data.pop(chm)
            for str in self.data[chm]:
                self.data[chm][str].resize(sizes[chm],refcheck=0)
    
    def sum(self):
        '''
        Description:
            calculate the total reads in the whole genome
        
        Parameter:
            None
        
        Value:
            an interger value representing the reads count.
        '''
        v=0
        for chm in self.data:
            for stra in self.data[chm]:
                v+=self.data[chm][stra].sum()
        return v

    def toWig(self,fs=None,extend=0,mifrsz=10,mafrsz=300):
        '''
        Description:
            Calculate nucleosome occupancy from the reads data
        
        Parameter:
            fs: average size of fragments that are subject to sequencing and generate the reads, only for signgle-end reads. When this value is not given, a fs value will be infered by the program. For paired-end reads loaded buy the function loadBedPaired(), set fs to 0.
            extend: a interger value, each read will be extend to this length.
            mifrsz: the minimal estimated average fragment size, only for single-end reads 
            mafrsz: the maximal estimated average fragment size, only for single-end reads
        
        Value: a Wig class instance
        '''
        step=self.step
        if fs==None:fs=self.fragSizeDis(minsize=mifrsz,maxsize=mafrsz)
        if extend<=0:extend=fs
        print('extend to',extend)
        old_extend=extend
        fragsize,extend=functions.div(fs,(2*step)),functions.div(extend,(2*step))
        wg=Wig(step=step)
        print('generating wig ...')
        for chm in self.data:
            tmax=max(1000,fragsize*4,extend*4)
            if self.data[chm]['+'].size<tmax:self.data[chm]['+'].resize(tmax,refcheck=0)
            if self.data[chm]['-'].size<tmax:self.data[chm]['-'].resize(tmax,refcheck=0)
            wg.addChr(chm)
            lth=int(self.data[chm]['+'].size)
            wg.resizeChr(chm,lth*step)
            self.data[chm]['+'][fragsize:lth]=self.data[chm]['+'][0:(lth-fragsize)]
            for i in range(fragsize):self.data[chm]['+'][i]=0
            self.data[chm]['+'][0:(lth-fragsize)]+=self.data[chm]['-'][fragsize:lth]
            for p in range(-extend,extend+1):wg.data[chm][extend:(lth-extend)]+=self.data[chm]['+'][(extend+p):(lth-extend+p)]
        wg.foldChange(functions.div(old_extend*1.0,wg.step)) ##### added by Kaifu on May29, 2014
        return wg
    def foldSampling(self,fold=1.0):
        '''
        Description:
            Randomly sample the reads count to a specified fold
            
        Parameter:
            fold: randomly sample the total reads count to the original reads count multiply by this fold value
        
        Value:
            None
        '''
        if fold>1:fold=fold-1
        wig=self
        tarray=numpy.array([0.0])
        for chm in wig.data:
            for str in wig.data[chm]:
                print(chm)
                tarray.resize(int(wig.data[chm][str].sum()),refcheck=0)
                tarray,csz,i,tsz=tarray*0,wig.data[chm][str].size,0,0
                while i<csz:
                    newtsz=int(tsz+wig.data[chm][str][i])
                    if newtsz>=tarray.size:tarray.resize(newtsz+1000,refcheck=0)
                    while tsz<newtsz:
                        tarray[tsz]=i
                        tsz+=1
                    i+=1
                i,tnum,tsz=0,int(wig.data[chm][str].sum()*fold),tsz
                if fold<=1:wig.data[chm][str]*=0
                while i<tnum:
                    i+=1
                    wig.data[chm][str][tarray[randint(0,tsz-1)]]+=1


if __name__ == "__main__":
    print('')
    import sys
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # This allow DANPOS to print each message on screen immediately.
