#!/usr/bin/env python
import numpy,re
from copy import deepcopy
from rpy2.robjects import r,FloatVector
from math import log10,sqrt,log
from summits import Summits
from time import time
import sys
import functions

class Wig:
    def __init__(self,file="",gfile='',step=0,suppress=False):
        '''
        Parameter:
            file: a path to a file in Wiggle format
            step: each chrosome will present as an vector, each value in the vector represent a short region of the chrosome, the step value is the size of the short region.
        '''
        self.data = {} # a dictionary in the format self.data[chr_name]=numpy_array
        self.step=step
        if gfile!='':
            for line in open(gfile):
                col=line.split()
                self.data[col[0]]=numpy.array([0.0])
                if step>0:self.data[col[0]].resize(functions.div(int(col[1]),step)+1,refcheck=0)
                else:self.data[col[0]].resize(int(col[1])+1,refcheck=0)
        if file !="": self.load(file,gfile=gfile,suppress=suppress)
    def absSum(self):
        '''
        Description:
            the sum of absolute value at each data point.
        
        Parameter:
            None
            
        Value:
            float value
        '''
        wig=self.data
        sum=0.0
        for chrom in wig:
            for v in wig[chrom]:sum+=abs(v)
        return sum
    def add(self,wig2):
        '''
        Description:
            add value at each data point between two Wig class instances.
        
        Parameter:
            wig2: a Wig class instance whose value will not change but will be add to the Wig class instance calling this method.
        
        Value:
            None
        '''
        if self.step!=wig2.step:
            wig2=deepcopy(wig2)
            wig2.changeStep(self.step)
        chrs={}
        for chrom in self.data:chrs[chrom]=1
        for chrom in wig2.data:chrs[chrom]=1
        chrs=list(chrs.keys())
        for chrom in chrs:
            if chrom not in self.data:self.addChr(chrom)
            if chrom not in wig2.data:wig2.addChr(chrom)
            lth1=self.data[chrom].size
            lth2=wig2.data[chrom].size
            if lth1>lth2:wig2.data[chrom].resize(lth1,refcheck=0)
            elif lth1<lth2:self.data[chrom].resize(lth2,refcheck=0)
            self.data[chrom]+=wig2.data[chrom]
        return True
    def addChr(self,chrom):
        '''
        Description:
            add a new chrosome to the existance Wig class instance
        
        Parameter:
            chrom: the name of the chrosome to be added
            
        Value: none
        '''
        self.data[chrom]=numpy.array([0.0])
    def ajust_size(self,gfile):
        step=self.step
        crs={}
        for line in open(gfile):
            col=line.split()
            size=functions.div(int(col[1]),step)+1
            crs[col[0]]=size
            if size<self.data[col[0]].size:self.data[col[0]]=self.data[col[0]][:size] #resize(int(col[1])/step+1,refcheck=0)
            if size>self.data[col[0]].size:self.data[col[0]].resize(size,refcheck=0)
        tcrs=list(self.data.keys())
        for cr in tcrs:
            if cr not in crs:self.data.pop(cr)
        
    def getChr(self,chrom):
        '''
        Description:
            retrieve a chrosome by name.
            
        Parameter:
            chrom: the name of the chrosome to be retrived.
        
        Value:
            a list as an instance of the numpy.array.
        '''
        return self.data[chrom]
    def correlation(self,wg2,ofile,mind=0,maxd=50):
        '''
        Description:
            Calculate correlation between the two Wig class instances allowing shift distances.
            
        Parameters:
            ofile: the output file used to save the result.
            mind: minimal shift distance
            maxd: maximal shift distance
        
        Value:
            None
        '''
        d={}
        chrs=self.getChrs()
        chrs2=wg2.getChrs()
        mind = functions.div(mind,self.step)
        maxd = functions.div(maxd,self.step)
        for chrom in chrs:
            d[chrom]={}
            if chrom in chrs2:
                cs=max(self.chrSize(chrom),wg2.chrSize(chrom))
                if self.chrSize(chrom)!=cs:self.resizeChr(chrom,cs)
                if wg2.chrSize(chrom)!=cs:wg2.resizeChr(chrom,cs)
                for td in range(int(mind),int(maxd+1)):
                    v=r.cor(FloatVector(self.getChr(chrom)[:(functions.div(cs,self.step)-td)]),FloatVector(wg2.getChr(chrom)[td:functions.div(cs,self.step)]))
                    d[chrom][td]=float(str(v).split()[-1])
        fo=open(ofile,'w')
        fo.write('Shift_disance\tcorrelation_coefficient\n')
        for td in range(round(mind),round(maxd+1)):
            v=0.0
            s=0
            gs=self.gsize()
            for chrom in d:
                v+=self.chrSize(chrom)*d[chrom][td]
            fo.write(str(td*self.step)+'\t'+str(functions.div(v,gs))+'\n')

    def bgsub(self,wig2,lmd=1000):
        '''
        Description:
            subtract the value of each data point in wig2 from self after smoothing of wig2
        
        Parameter:
            wig2: the Wig class instance used to do the subtraction
            lmd: the bin size used to smooth wig2
        
        Value:
            None
        
        Note:
            a copy of wig2 is smoothed and used to do the subtraction, wig2 will not change.
        '''
        wig2=deepcopy(wig2)
        if self.step!=wig2.step:wig2.changeStep(self.step)
        chrs={}
        for chrom in self.data:chrs[chrom]=1
        for chrom in wig2.data:chrs[chrom]=1
        chrs=list(chrs.keys())
        for chrom in chrs:
            lth1=len(self.data[chrom])
            lth2=len(wig2.data[chrom])
            if lth1>lth2:wig2.data[chrom].resize(lth1,refcheck=0)
            else:self.data[chrom].resize(lth2,refcheck=0)
        if lmd>0:wig2.smooth(lmd=lmd)
        wig2.foldChange(functions.div(self.sum(),wig2.sum()))
        sys.stdout.write('before subtracting: %f\n' % self.sum())
        self.subtract(wig2)
        sys.stdout.write('after subtracting: %f\n' % self.sum())
        self.rvNeg()
        sys.stdout.write('after removing negtive values: %f\n' % self.sum())
        return True
    def callRegions(self,ofile=None,width=0,distance=165,pheight=0,height=0,calculate_P_value=1,mode='w',title_line=1,pos_only=False,fold=0,suppress=False,fdr=False,fdrSampleSize=1000000):
        '''
        Description:
            This fuction is designed to call broad peaks, such as histone modification peaks.
        
        Parameter:
            ofile: a path to the file used to save the peaks.
            width: minimal width of peaks
            distance: minimal distance between peaks, neighboring peaks with distance shorter than this value will be merged.
            pheight: a P value cutoff used to call peaks.
            height: the occupancy cutoff used to call peaks. valide only when pheight is set to 1.
            calculate_P_value: calculate P value for each peak if set to 1, else set to 0.
            mode: the mode to write result to ofile, could be either 'w' to create new file or 'a' to append to a existing file
            title_line: set to 1 if need a title line in the result file ofile
        
        Value:
            pks[chrosome_name][start_position]=end_position
            
        '''
        #ppois=r('''function(q,avg){return(ppois(q,avg,lower.tail=FALSE,log=TRUE))}''')
        m=self.mean()
        if height==0 and pheight!=0:
            #if fold!=0:height=m*fold
            #else:
            height=int(m+0.5)
            while( (0-functions.div(float(functions.ppois(height,m.item()).split()[-1]),log(10))) < pheight):height+=1
        #if not suppress:
        sys.stdout.write('whole genome aveage value is ' +
            str(m)+', use calling cutoff '+str(height) + "\n")
        dic=self.data
        step=self.step
        twidth=functions.div(width,step)
        if twidth<1:twidth=1
        pks={}
        if not suppress: sys.stdout.write('calling ...\n')
        for chrom in dic:
            pks[chrom]={}
            lth=dic[chrom].size#.chrSize(chrom)/twstep-twidth
            if lth<=0:continue
            if not suppress: sys.stdout.write(chrom + ' ')
            start,pos=-1,0
            while pos <lth:
                if start<0:
                    if dic[chrom][pos]>=height:start=pos
                elif dic[chrom][pos]<height:
                    if pos-start>=twidth:pks[chrom][start*step]=(pos)*step
                    start=-1
                pos+=1
            if not suppress: sys.stdout.write("%i\n" % len(pks[chrom]))
        if distance>0: 
            if not suppress: sys.stdout.write('mering\n')
            crs=list(pks.keys())
            crs.sort()
            for chrom in crs :
                ps=list(pks[chrom].keys())
                if not suppress: sys.stdout.write(chrom + ' from ' + str(len(ps)) +' ')
                ps.sort()
                lth=len(ps)-1
                i=0
                while i<lth:
                    #j=i+1
                    #while ps[j]-pks[chrom][ps[i]]<=distance:
                    if ps[i+1]-pks[chrom][ps[i]]<=distance:# or pks[chrom][ps[i+1]]<=pks[chrom][ps[i]]:
                        if pks[chrom][ps[i]]<pks[chrom][ps[i+1]]:pks[chrom][ps[i]]=pks[chrom][ps[i+1]]
                        pks[chrom].pop(ps[i+1])
                        ps[i+1]=ps[i]
                    i+=1
                if not suppress: sys.stdout.write('to '+str(len(pks[chrom])) + "\n")
        if ofile!=None:
            outf=open(ofile,mode)
            if title_line:
                if mode=='w':
                    if calculate_P_value:outf.write("chr\tstart\tend\tsummit_pos\tsummit_value\tstrand\ttotal_signal\twidth_above_cutoff\tsummit_minus_logP\n")
                    else:outf.write("chr\tstart\tend\tsummit_pos\tsummit_value\tstrand\ttotal_signal\twidth_above_cutoff\n")
            crs=list(pks.keys())
            crs.sort()
            for chrom in crs :
                #out[chrom]={}
                ps=list(pks[chrom].keys())
                ps.sort()
                for p in ps:
                    width_above_cutoff=0
                    s,e=functions.div(p,step),functions.div(pks[chrom][p],step)
                    v=max(dic[chrom][s:e])
                    auc=0#dic[chrom][s:e].sum()*step
                    if ((e-s)>=twidth):# and (v>=sheight):
                        smts=[]
                        for i in range(s,e):
                            if dic[chrom][i]==v:smts.append(str(i*step))
                            if dic[chrom][i]>=height:width_above_cutoff+=step
                            if pos_only:
                                if dic[chrom][i]>0:auc+=dic[chrom][i]
                            else:auc+=dic[chrom][i]
                        auc=auc*step
                        smt=','.join(smts)
                        lth=functions.div(len(smts),2)
                        if calculate_P_value:
                            pvl=functions.div(float(functions.ppois(v.item(),m.item()).split()[-1]),log(10))
                            outf.write(chrom+"\t"+str(p)+"\t"+str(pks[chrom][p])+"\t"+smt+'\t'+str(v)+'\t+\t'+str(auc)+'\t'+str(width_above_cutoff)+'\t'+str(0-pvl)+"\n")
                        else:outf.write(chrom+"\t"+str(p)+"\t"+str(pks[chrom][p])+"\t"+smt+'\t'+str(v)+'\t+\t'+str(auc)+'\t'+str(width_above_cutoff)+'\n')
                    else:pks[chrom].pop(p)
        return pks
        
    def fillRegions(self,regions={},file=None,pheight=1,height=0,width=0,calculate_P_value=1,pos_only=False,suppress=False):
        '''
        Add description here
        '''
        #ppois=r('''function(q,avg){return(ppois(q,avg,lower.tail=FALSE,log=TRUE))}''')
        m=self.mean()
        if height==0 and pheight!=0:
            #if fold!=0:height=m*fold
            #else:
            height=int(m+0.5)
            while( (0-functions.div(float(functions.ppois(height,m.item())).split()[-1],log(10))) < pheight):height+=1
        #if not suppress:
        sys.stdout.write('whole genome aveage value is '+str(m)+', use cutoff '+ str(height) + "\n")
        #lines=open(file).readlines()
        outf=open(file,'w')
        if calculate_P_value==1:outf.write('chr\tstart\tend\tcenter\twidth_above_cutoff\ttotal_signal\theight\theight_logP\n')
        else:outf.write('chr\tstart\tend\tcenter\twidth_above_cutoff\ttotal_signal\theight\n')
        step=self.step
        dic=self.data
        total_width_above_cutoff=0
        #for line in lines[1:]:
        #    col=line.split()
        for chrom in regions:
            starts=list(regions[chrom].keys())
            starts.sort()
            for start in starts:
                if regions[chrom][start]-start<width:continue
                s,e=functions.div(start,step),functions.div(regions[chrom][start],step)+1
                if dic[chrom].size<=e:dic[chrom].resize(e+1,refcheck=0)
                width_above_cutoff=0
                v=max(dic[chrom][s:e])
                if pos_only:v=max(v,0)
                auc=0
                #substart=-1
                #subpeaks=''
                i=s
                while i<e:
                    '''
                    if substart<0:
                        if dic[chrom][i]>=height:substart=i
                    elif dic[chrom][i]<height:
                        subpeaks+=str(substart*step)+':'+str(i*step)+','
                        substart=-1
                    '''
                    if dic[chrom][i]>=height:width_above_cutoff+=step
                    if pos_only:
                        if dic[chrom][i]>0:auc+=dic[chrom][i]
                    else:auc+=dic[chrom][i]
                    i+=1
                #if substart>=0:subpeaks+=str(substart*step)+':'+str((e-1)*step)+','
                auc=auc*step
                total_width_above_cutoff+=width_above_cutoff
                if calculate_P_value:
                    pvl=functions.div(float(functions.ppois(v.item(),m.item()).split()[-1]),log(10))
                    outf.write(chrom + '\t' + str(start) + '\t' + str(regions[chrom][start]) +
                        '\t' + str(functions.div((regions[chrom][start] + start),2)) + '\t' + str(width_above_cutoff) +
                         '\t' + str(auc) + '\t' + str(v) + '\t' + str(0-pvl) + "\n")
                else:
                    outf.write(chrom + '\t' + str(start) + '\t' + str(regions[chrom][start]) + '\t' + 
                        str(functions.div((regions[chrom][start]+start),2)) + '\t' + str(width_above_cutoff) + '\t' + 
                        str(auc) + '\t' + str(v) + '\n')
        sys.stdout.write('total_width_above_cutoff: '+str(total_width_above_cutoff) + "\n")
        outf.close()
        
    def callPositions(self,ofile,width=40,distance=165,edge=1,pcut=1e-5,height=0,fill_gap=False,fill_value=1,calculate_P_value=1,mode='w',title_line=1,poscal=0,regions=None,rd=None):
        '''
        Description:
            This fuction is designed to call nucleosome positions
            
        Parameter:
            ofile: a path to the file used to save the positions.
            width: minimal width of positions
            distance: minimal distance between positions, neighboring positions with distance shorter than this value will be merged.
            edge: set to 1 if need to search for position edges, else set to 0
            pcut: a P value cutoff used to call positions.
            height: the occupancy cutoff used to call positions. valide only when pheight is set to 1.
            fill_gap: fill the gap between two neighboring nuclesomes with a new nucleosome if the gap size is reasonable.
            fill_value: the default value to be set to a filled nucleosome
            calculate_P_value: calculate P value for each position if set to 1, else set to 0.
            mode: the mode to write result to ofile, could be either 'w' to create new file or 'a' to append to a existing file
            title_line: set to 1 if need a title line in the result file ofile
            poscal: set to 1 if need to calculate nucleosome positioning score and P value,else set to 0
        
        Value:
            None.
            
        '''
        outf=open(ofile,mode)
        if mode=='w':
            if poscal>0:
                if calculate_P_value:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\tsmt_log10pval\tfuzziness_score\tfuzziness_log10pval\n")
                else:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\tfuzziness_score\n")
            else:
                if calculate_P_value:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\tsmt_log10pval\n")
                else:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\n")
        twig=self
        #ppois=r('''function(q,avg){return(ppois(q,avg,lower.tail=FALSE,log.p=TRUE)/log(10))}''')
        m=self.mean()
        if height==0 and pcut!=0:
            height=int(m+0.5)
            while( (0-float(functions.ppois(height,m.item())/log(10).split()[-1])) < pcut):height+=1
        #if not suppress:
        sys.stdout.write('whole genome aveage value is '+str(m)+', use calling cutoff' + str(height) + "\n")
            
        sys.stdout.write('calling summits ...\n')
        smts=twig.callSummits(width=width,pcut=1,height=height,regions=regions)
        sys.stdout.write('merging summits ...\n')
        smts.merge(wg=twig,distance=distance)
        if fill_gap:
            sys.stdout.write('filling gaps ...\n')
            smts.fillgap(wg=twig,height=height,distance=distance)
        
        
        if poscal>0:
            sys.stdout.write('calculating positioning score\n')
            smts.positioning(twig,rd=rd)
        if edge:sys.stdout.write('searching position edges ...\n')
        else:sys.stdout.write('saving positions ...\n')
        width= functions.div(width,2*twig.step)
        dic=smts.data
        rhalfdis=functions.div(distance,2)
        halfdis=functions.div(distance,(2*twig.step))
        for chrom in dic:
            lth=functions.div(twig.chrSize(chrom),twig.step)-1
            sys.stdout.write(chrom+"\n")
            if lth==0:continue
            poses=dic[chrom]['p']
            valus=dic[chrom]['v']
            if poscal>0:
                positioning=dic[chrom]['s']
                ppos=dic[chrom]['ppos']
            tlen=poses.size
            if calculate_P_value:
                pvs=functions.ppois(FloatVector(dic[chrom]['v']),m.item())/log(10)
            i=0
            while i<tlen:
                pos=poses[i]
                start,end=0,0
                if edge==0:start,end=pos-74+(74%self.step)+1,pos+74-(74%self.step)+1
                else:
                    ppp=functions.div(pos,twig.step)
                    p=ppp-1
                    while(start==0):
                        if p<=width:start=1
                        if ppp-p>=functions.div(halfdis,2):
                            if twig.data[chrom][p]==min(twig.data[chrom][(p-width):(p+width+1)]):start=p*twig.step+1 
                            elif twig.data[chrom][p]>valus[i]:start=p*twig.step+1
                            elif twig.data[chrom][p]<height:start=p*twig.step+1
                        elif i>0 and p*twig.step==poses[i-1]:start=p*twig.step+1 
                        p-=1
                    p=ppp+1
                    while(end==0):
                        if (p+width)>=lth:end=lth*twig.step+1
                        elif p-ppp>=functions.div(halfdis,2):
                            if twig.data[chrom][p]==min(twig.data[chrom][(p-width):(p+width+1)]):end=p*twig.step+1
                            elif twig.data[chrom][p]>valus[i]:end=p*twig.step+1
                            elif twig.data[chrom][p]<height:end=p*twig.step+1
                        elif i<(tlen-1) and p*twig.step==poses[i+1]:end=p*twig.step+1 
                        p+=1
                if end>start:
                    if poscal>0:
                        if calculate_P_value:
                            outf.write(chrom+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\t"+str(pvs[i])+"\t"+str(positioning[i])+"\t"+str(ppos[i])+"\n")
                        else:outf.write(chrom+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\t"+str(positioning[i])+"\n")
                    else:
                        if calculate_P_value:
                            outf.write(chrom+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\t"+str(pvs[i])+"\n")
                        else:outf.write(chrom+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\n")
                i+=1
        outf.close()
        return smts
        

    def fillPositions(self,dic,file,width=40,distance=165,edge=1,pcut=1e-5,height=5,calculate_P_value=1,mode='w',title_line=1,poscal=0,rd=None):
        '''
        Description:
            This fuction is designed to call nucleosome positions
            
        Parameter:
            ofile: a path to the file used to save the positions.
            width: minimal width of positions
            distance: minimal distance between positions, neighboring positions with distance shorter than this value will be merged.
            edge: set to 1 if need to search for position edges, else set to 0
            pcut: a P value cutoff used to call positions.
            height: the occupancy cutoff used to call positions. valide only when pheight is set to 1.
            fill_gap: fill the gap between two neighboring nuclesomes with a new nucleosome if the gap size is reasonable.
            fill_value: the default value to be set to a filled nucleosome
            calculate_P_value: calculate P value for each position if set to 1, else set to 0.
            mode: the mode to write result to ofile, could be either 'w' to create new file or 'a' to append to a existing file
            title_line: set to 1 if need a title line in the result file ofile
            poscal: set to 1 if need to calculate nucleosome positioning score and P value,else set to 0
        
        Value:
            None.
            
        '''
        twig=self
        smts=Summits()
        smts.data=deepcopy(dic)
        '''
        smts=Summits()
        for line in open(file).readlines()[1:]:
            col=line.split()
            if not smts.data.has_key(col[0]):
                smts.data[col[0]]={}
                smts.data[col[0]]['p']=[]
            smts.data[col[0]]['p'].append(int(col[3])-1)
        for cr in smts.data:
            smts.data[cr]['p'].sort()
            smts.data[cr]['p']=numpy.array(smts.data[cr]['p'])
        '''
        
        smts.fetchValueFromWig(twig)
        if poscal>0:
            sys.stdout.write('calculating positioning score\n')
            smts.positioning(twig,rd=rd)
            
        outf=open(file,mode)
        
        if mode=='w':
            if poscal>0:
                if calculate_P_value:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\tsmt_log10pval\tfuzziness_score\tfuzziness_log10pval\n")
                else:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\tfuzziness_score\n")
            else:
                if calculate_P_value:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\tsmt_log10pval\n")
                else:outf.write("chr\tstart\tend\tsmt_pos\tsmt_value\n")
        #ppois=r('''function(q,avg){return(ppois(q,avg,lower.tail=FALSE,log.p=TRUE)/log(10))}''')
        m=self.mean()
        if height==0 and pcut!=0:
            height=int(m+0.5)
            while( (0-float(functions.ppois(height,m.item())/log(10).split()[-1])) < pcut):height+=1
        #print 'calling summits ...'
        #smts=twig.callSummits(width=width,pcut=1,height=height,regions=regions)
        #print 'merging summits ...'
        #smts.merge(wg=twig,distance=distance)
        #if fill_gap:
        #    print 'filling gaps ...'
        #    smts.fillgap(wg=twig,height=height,distance=distance)
        if edge:sys.stdout.write('searching position edges ...\n')
        else:sys.stdout.write('saving positions ...\n')
        width=functions.div(width,2*twig.step)
        dic=smts.data
        rhalfdis=functions.div(distance,2)
        halfdis=functions.div(distance,(2*twig.step))
        for chrom in dic:
            lth=functions.div(twig.chrSize(chrom),twig.step)-1
            sys.stdout.write(chrom+"\n")
            if lth==0:continue
            poses=dic[chrom]['p']
            valus=dic[chrom]['v']
            if poscal>0:
                positioning=dic[chrom]['s']
                ppos=dic[chrom]['ppos']
            tlen=poses.size
            if calculate_P_value:
                pvs=functions.ppois(FloatVector(dic[chrom]['v']),m.item())/log(10)
            i=0
            while i<tlen:
                pos=poses[i]
                '''
                if valus[i]<height:
                    i+=1
                    continue
                '''
                start,end=0,0
                if edge==0:start,end=pos-74+(74%self.step)+1,pos+74-(74%self.step)+1
                else:
                    ppp=functions.div(pos,twig.step)
                    p=ppp-1
                    while(start==0):
                        if p<=width:start=1
                        if ppp-p>=functions.div(halfdis,2):
                            if twig.data[chrom][p]==min(twig.data[chrom][(p-width):(p+width+1)]):start=p*twig.step+1 
                            elif twig.data[chrom][p]>valus[i]:start=p*twig.step+1
                            elif twig.data[chrom][p]<height:start=p*twig.step+1
                        elif i>0 and p*twig.step==poses[i-1]:start=p*twig.step+1 
                        p-=1
                    p=ppp+1
                    while(end==0):
                        if (p+width)>=lth:end=lth*twig.step+1
                        elif p-ppp>=functions.div(halfdis,2):
                            if twig.data[chrom][p]==min(twig.data[chrom][(p-width):(p+width+1)]):end=p*twig.step+1
                            elif twig.data[chrom][p]>valus[i]:end=p*twig.step+1
                            elif twig.data[chrom][p]<height:end=p*twig.step+1
                        elif i<(tlen-1) and p*twig.step==poses[i+1]:end=p*twig.step+1 
                        p+=1
                if end>start:
                    if poscal>0:
                        if calculate_P_value:
                            outf.write(chrom+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\t"+str(pvs[i])+"\t"+str(positioning[i])+"\t"+str(ppos[i])+"\n")
                        else:outf.write(chrom+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\t"+str(positioning[i])+"\n")
                    else:
                        if calculate_P_value:
                            outf.write(chrom+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\t"+str(pvs[i])+"\n")
                        else:outf.write(chrom+"\t"+str(start)+"\t"+str(end)+"\t"+str(pos+1)+"\t"+str(valus[i])+"\n")
                i+=1
        outf.close()
        return smts
        

    def callSummits(self,width=40,pcut=1,height=5,regions=None):
        '''
        Description:
            call occupancy summits using a sliding window.
        
        Parameter:
            width: the width of the sliding window used to call summits.
            pcut:  a P value cutoff used to call positions.
            height: the occupancy cutoff used to call positions. valide only when pheight is set to 1.
            regions: A set of regions in which the summits are to be defined, regions[chromatin_name][start_position]=end_position
        '''
        if pcut!=1:
            qp=r('''function(p,avg){return(qpois(log(p),avg,lower.tail=FALSE,log=TRUE))}''')
            m=self.mean()
            height=int(str(qp(pcut,m)).split()[-1])
            sys.stdout.write('set summit calling cutoff to'+str(height)+"\n")
        smts=Summits()
        width=functions.div(width,2)
        
        if regions==None:
            regions={}
            for cr in self.data:
                regions[cr]={}
                regions[cr][width]=self.data[cr].size*self.step-1
        
        for cr in regions:#self.data:
            sys.stdout.write(str(cr) + "\n")
            poses=numpy.array([0])
            valus=numpy.array([0.0])
            dic={}
            lst=self.data[cr]
            step=self.step
            width=functions.div(width,step)
            backwidth=width+1
            #sz=lst.size
            region_poses=list(regions[cr].keys())
            region_poses.sort()
            num=0
            for p in region_poses:
                i,dlth=functions.div(p,step),functions.div(regions[cr][p],step)
                while i<dlth:#lst.size:
                    v=lst[int(i-width):int(i+backwidth)].max()
                    if v<height:
                        i+=backwidth#continue
                    elif v==lst[int(i)]:
                        #dic[i*step]=lst[i]
                        if num>=poses.size-1:
                            poses.resize(num+1000,refcheck=0)
                            valus.resize(num+1000,refcheck=0)
                        ti,tstart,tend,v1,v2=i,i,i,v,v
                        
                        while v1==v:
                            tstart-=1
                            if tstart>=0:
                                v1=lst[int(tstart)]
                            else:
                                v1,tstart=v-1,0
                        while v2==v:
                            tend+=1
                            if tend<dlth:
                                v2=lst[int(tend)]
                            else:
                                v2,tend=v-1,dlth-1
                        
                        ti=functions.div((tstart+tend),2)
                        i=tend
                        poses[num]=ti*step
                        valus[num]=v
                        num+=1
                        i+=backwidth
                        if i<tend:i=tend
                    else:i+=1
            dic['p']=poses[:num]
            dic['v']=valus[:num]
            smts.data[cr]=dic#self[cr].summit(width=width,height=height)#dic
        return smts
    def changeStep(self,step):
        '''
        Description:
            change the step size.
            Note: The new value for each step is determined by sampling an old value within the step!
        
        Parameter:
            step: the new step that is to be setted as.
        
        Value:
            None
        '''
        if step==self.step:return True
        sys.stdout.write('change wiggle step from ' + str(self.step) + 
            ' to ' + str(step) + "\n")
        if step>self.step:
            if step%self.step!=0:
                sys.stdout.write('Wrong: the fold change between new step ' +
                    str(step) + ' and old step ' + str(self.step) + ' is not integer!' +
                     "\n")
                return False
            fd=functions.div(step,self.step)
            for chrom in self.data:
                lth=self.data[chrom].size
                nlth=functions.div(lth,fd)
                lst=numpy.array([0.0])
                lst.resize(nlth,refcheck=0)
                for npos in range(0,nlth):
                    for pos in range(npos*fd,(npos+1)*fd):lst[npos]+=self.data[chrom][pos]
                self.data[chrom]=lst
                self.data[chrom]=functions.div(self.data[chrom],fd)
        if step<self.step:
            if self.step%step!=0:
                sys.stdout.write('Wrong: the fold change between new step ' + 
                    str(step) + ' and old step ' + str(self.step) + 
                    ' is not integer! ' + str(int(self.step)%step) + "\n")
                return False
            fd=functions.div(self.step,step)
            for chrom in self.data:
                lth=self.data[chrom].size
                nlth=lth*fd
                lst=numpy.array([0.0])
                lst.resize(nlth,refcheck=0)
                for pos in range(0,lth):
                    for npos in range(pos*fd,(pos+1)*fd):lst[npos]=self.data[chrom][pos]
                self.data[chrom]=lst
                self.data[chrom]=self.data[chrom]
        self.step=step
        return True
    
    def chisqTest(self,wig2,tchr=''):
        '''
        Description:
            do Chi-square test to calculate differential signial for each data point between two Wig class instances.
        
        Parameter:
            wig2: a Wig class instance to be compared to
            tchr: specify a chrosome that is to be compared, leave it to '' if want to do for all chrosomes.
            
        Value:
            A Wig class instance that cantain data for the differential signial
        
        '''
        r.options(warn=-2)
        wig1=deepcopy(self)
        wig2=deepcopy(wig2)
        if wig1.step!=wig2.step:wig2.changeStep(wig1.step)
        out=Wig(step=wig1.step)
        sum1=wig1.sum()
        sum2=wig2.sum()
        tlen=0
        for chrom in wig1.data:tlen+=len(wig1.data[chrom])
        donenum=0
        ctime=time()
        for chrom in wig1.data:
            if tchr!='' and chr!=tchr:continue
            if chrom in wig2.data:
                wig1.data[chrom]+=1
                wig2.data[chrom]+=1
                out.data[chrom]=numpy.array(0.0)
                lth1=len(wig1.data[chrom])
                lth2=len(wig2.data[chrom])
                lth=max(lth1,lth2)
                if lth>lth2:wig2.data[chrom].resize(lth,refcheck=0)
                elif lth>lth1:wig1.data[chrom].resize(lth,refcheck=0)
                out.data[chrom].resize(lth,refcheck=0)
                for p in range(0,lth):
                    #if wig1.data[chrom][p]<1:wig1.data[chrom][p]=1
                    #if wig2.data[chrom][p]<1:wig2.data[chrom][p]=1
                    vec=r.c(wig1.data[chrom][p],wig2.data[chrom][p],sum1-wig1.data[chrom][p],sum2-wig2.data[chrom][p])
                    if (wig1.data[chrom][p]==wig2.data[chrom][p]):out.data[chrom][p]=0
                    else:
                        test=r['chisq.test'](r.matrix(vec,nrow=2))
                        pvl=float(str(test[2]).split()[-1])
                        if pvl<1e-323:pvl=1e-323
                        if wig1.data[chrom][p]>wig2.data[chrom][p]:out.data[chrom][p]=-log10(pvl)
                        else:out.data[chrom][p]=log10(pvl)
                    donenum+=1
                    if time()-ctime>10:
                        ctime=time()
                        sys.stdout.write(str(functions.div(donenum*100.0,tlen)) + 
                            " percent of " + str(tlen) +" done\n")
        wig1.clearEmptyEnd()
        wig2.clearEmptyEnd()
        return out

    def chrSize(self,chrom):
        '''
        Description:
            retrive chrosome size by name
            
        Parameter:
            chr: the name of chrosome whose size is to be retrived
        Value:
            Interger value (step size has been multiplied)
            
        '''
        if chrom in self.data:return self.data[chrom].size*self.step
        else:return 0
    def chrSum(self,chrom):
        '''
        Description:
            Retrieve the sum of occupancy by chrosome
        Parameter:
            chrom: the name of chrosome whose occupancy sum is to be retrieved
        Value:
            Float value
        '''
        return self.data[chrom].sum()*self.step
    def clearEmptyEnd(self):
        '''
        Description:
            Clear the 0 values at the end of each chrosome
        Parameter:
            None
        Value:
            None
        '''
        for chrom in self.data:
            size=len(self.data[chrom])
            while self.data[chrom][size-1]==0 and size>0:size-=1#self.data[chrom].pop()
            self.data[chrom].resize(size,refcheck=0)
    def divideAndLog2(self,wig2):
        '''
        Description:
            divid by wig2 at each data point and then transform the resultant value by log2
        
        Parameter:
            wig2: the Wig class instance that will be used to devide
            
        Value:
            A Wig class instance
        '''
        wig1=deepcopy(self)
        wig2=deepcopy(wig2)
        if wig1.step!=wig2.step:
            wig2.changeStep(wig1.step)
        chrs={}
        for chrom in wig1.data:chrs[chrom]=1
        for chrom in wig2.data:chrs[chrom]=1
        chrs=list(chrs.keys())
        for chrom in chrs:
            lth1=len(wig1.data[chrom])
            lth2=len(wig2.data[chrom])
            if lth1>lth2:wig2.data[chrom].resize(lth1,refcheck=0)
            else:wig1.data[chrom].resize(lth2,refcheck=0)
            lth=max(lth1,lth2)
            wig1.data[chrom]+=1
            wig2.data[chrom]+=1
            wig1.data[chrom] = functions.div(wig1.data[chrom], wig2.data[chrom])
            wig1.data[chrom][0:lth]=r.log2(FloatVector(wig1.data[chrom]))[0:lth]
        return wig1
    def divide(self,wig2):
        '''
        Description:
            divid by wig2 at each data point and then transform the resultant value by log2
        
        Parameter:
            wig2: the Wig class instance that will be used to devide
            
        Value:
            A Wig class instance
        '''
        wig1=deepcopy(self)
        wig2=deepcopy(wig2)
        if wig1.step!=wig2.step:
            wig2.changeStep(wig1.step)
        chrs={}
        for chrom in wig1.data:chrs[chrom]=1
        for chrom in wig2.data:chrs[chrom]=1
        chrs=list(chrs.keys())
        for chrom in chrs:
            lth1=len(wig1.data[chrom])
            lth2=len(wig2.data[chrom])
            if lth1>lth2:wig2.data[chrom].resize(lth1,refcheck=0)
            else:wig1.data[chrom].resize(lth2,refcheck=0)
            lth=max(lth1,lth2)
            #wig1.data[chrom]+=1
            wig2.data[chrom]+=1
            wig1.data[chrom]= functions.div(wig1.data[chrom], wig2.data[chrom])
            #wig1.data[chrom][0:lth]=r.log2(FloatVector(wig1.data[chrom]))[0:lth]
        return wig1
    def dfTest(self,cwig,test='C'):
        '''
        Description:
            Do differential test wit cwig
        Parameter:
            cwig: the Wig class instance which is to be tesed to.
            test: the statistical method that will be used to do the differential test
        Value:
            a Wig class instance containing the differential signal data
        '''
        if test=='C':
            sys.stdout.write('Chi-square test\n')
            pwig=self.chisqTest(cwig)
            return pwig
        elif test=='F':
            pwig=self.divideAndLog2(cwig)
            return pwig
        elif test=='P':
            sys.stdout.write('Poisson test\n')
            pwig=self.ppois(cwig)
            return pwig
        elif test=='S':
            sys.stdout.write('subtraction\n')
            pwig=deepcopy(self)
            pwig.subtract(cwig)
            return pwig
        elif test=='N':
            sys.stdout.write("No test method appointed, will not do any differential test.\n")
            return Wig(step=self.step)
        else:
            sys.stdout.write("Normalization method "+str(test)+" not applicable now\n")
            return Wig(step=self.step)

    def fisherTest(self,wig2):
        '''
        Description:
            Do Fisher's exact test to wig2.
            
        Parameter:
            wig2: the Wig class instance which is to be tesed to.
            
        Value:
            a Wig class instance containing the differential signal data
        '''
        r.options(warn=-1)
        if self.step!=wig2.step:
            wig2=deepcopy(wig2)
            wig2.changeStep(self.step)
        out=Wig(step=wig1.step)
        sum1=int(self.sum())
        sum2=int(wig2.sum())
        tlen=0
        for chrom in self.data:tlen+=len(self.data[chrom])
        donenum=0
        for chrom in self.data:
            if chr in wig2.data:
                out.data[chrom]=numpy.array(0.0)
                lth1=len(self.data[chrom])
                lth2=len(wig2.data[chrom])
                lth=max(lth1,lth2)
                if lth>lth2:wig2.data[chrom].resize(lth1,refcheck=0)
                elif lth>lth1:self.data[chrom].resize(lth2,refcheck=0)
                out.data[chrom].resize(lth,refcheck=0)
                for p in range(0,lth):
                    vec=r.c(int(self.data[chrom][p]),int(wig2.data[chrom][p]),sum1-int(self.data[chrom][p]),sum2-int(wig2.data[chrom][p]))
                    if (self.data[chrom][p]==wig2.data[chrom][p]):out.data[chrom][p]=1
                    else:
                        test=r['fisher.test'](r.matrix(vec,nrow=2))
                        pvl=float(str(test).split()[17])
                        if pvl<1e-323:pvl=1e-323
                        if self.data[chrom][p]>wig2.data[chrom][p]:out.data[chrom][p]=-log10(pvl)
                        else:out.data[chrom][p]=log10(pvl)
                    donenum+=1
                    if donenum%1000==0:sys.stdout.write(str(functions.div(donenum*100.0,tlen)) +
                        " percent of " + str(tlen) + "done\n")
        self.clearEmptyEnd()
        wig2.clearEmptyEnd()
        return out
    def foldChange(self,fold):
        '''
        Description:
            Do fold change at each data point.
            
        Parameter:
            fold: the value that will be multiplied by each data point
            
        Value:
            None
        '''
        for chrom in self.data:self.data[chrom]*=fold

    def getChrs(self):
        '''
        Description:
            Retrive all chrosome names
        Parameter:
            None
        Value:
            a list of chrosome names
        '''
        return list(self.data.keys())
        
    def gsize(self):
        '''
        Description:
            Calculate genome size
        Parameter:
            None
        Value:
            Interger value
        '''
        lth=0
        for chrom in self.data:
            lth+=self.chrSize(chrom)#self.step*self.data[chrom].size
        return lth
    def histogram(self,bnum=100000,nonzero_end=False):
        ma,mi=self.maxmin(nonzero=nonzero_end)
        counts,bins=numpy.histogram([],bins=bnum,range=(mi,ma))
        for cr in self.data:counts+=numpy.histogram(self.data[cr],bins=bnum,range=(mi,ma))[0]
        return [counts,bins]

    def load(self,file,gfile,suppress=False):
        '''
        Description:
            Load Wig class instance from Wiggle format file
        Parameter:
            file: a path to the file containing the data
        Value:
            None
        '''
        if file[-3:]=='wig':fi=open(file)
        if file[-6:]=='wig.gz':
            import gzip
            fi=gzip.open(file)
        for line in fi:
            if line[0]=='t':continue
            elif line[0:9]=='fixedStep':
                self.loadFixed(file,gfile=gfile,suppress=suppress)
                return
            elif line[0:3]=='var':
                self.loadVar(file,gfile=gfile,suppress=suppress)
                return
            else:sys.stdout.write("Load failure: format not recoganized! " + line + "\n")
    
    def loadFixed(self,file,gfile,suppress=False):#add by Kaifu April 4, 2012
        '''
        Description:
            load data from Fixed wiggle format file
        Parameter:
            file:a path to the file containging the data
            suppress: suppress waring message? True or False
        Value:
            None
        '''
        sss=time()
        
        ########## ---start---add by kaifu on Aug 14, 2012 ##########
        if self.step<1:
            tempf=open(file)
            if not suppress: sys.stdout.write('detecting step size ...\n')
            while self.step<1:
                line=tempf.readline()
                if line=='\n':continue
                if line[0]=='f':
                    col=line.split()
                    for term in col[1:]:
                        kv=term.split("=")
                        if kv[0]=='step':self.step=int(kv[1])
        ########## ---end---add by kaifu on Aug 14, 2012 ##########
        
        if not suppress: sys.stdout.write('parsing from ' + file + "\n")
        chrom=''
        chrs=['chr1']############ test #############
        pn=1
        pos=-1
        for line in open(file):
            if line=='\n':continue
            if line[0]=='t':continue
            elif line[0]=='f':
                col=line.split()
                for term in col[1:]:
                    kv=term.split("=")
                    if kv[0]=='chrom':
                        newchr=kv[1]
                    elif kv[0]=='start':newstart=int(kv[1])
                    elif kv[0]=='step':instep=int(kv[1])
                
                if not suppress:sys.stdout.write(str(newchr) + 
                    ' start from position '+ str(newstart) + "\n")#/instep
                if chrom != '':
                    if chrom not in self.data:self.data[chrom]=numpy.array([0.0])
                    nplst=numpy.array(lst)#[0:]
                    if instep!=self.step:
                        lth=len(lst)*functions.div(instep*1.0,self.step)
                        nplst=numpy.array([0.0])
                        nplst.resize(int(lth),refcheck=0)
                        for pos in range(len(lst)):
                            #nplst[pos*instep/self.step]+=lst[pos]########## deleted by kaifu on Sep 6, 2012 ##########
                            if instep<self.step:nplst[functions.div(pos*instep,self.step)]+=lst[pos] ########## add by kaifu on Sep 6, 2012 ##########
                            else:
                                for tpos in range(functions.div(pos*instep,self.step),functions.div((pos+1)*instep,self.step)):nplst[functions.div(tpos,self.step)]+=functions.div(lst[pos]*1.0*self.step,instep) ########## add by kaifu on Sep 6, 2012 ##########
                    self.data[chrom].resize(functions.div(start,self.step)+len(nplst),refcheck=0)
                    self.data[chrom][(functions.div(start,self.step)):]=nplst
                lst=[]
                chrom=newchr
                start=newstart-1
            else:
                try:
                    #if line=='nan\n':line='0\n'
                    lst.append(float(line.split()[0]))
                except:
                    if line[:3]=='nan':lst.append(0.0)
                    else:sys.stdout.write('wrong line: ' + str(line[:-1]) + "\n")
                    continue
                #print line.split()[0],lst[-1]
                #lst.append(value)
        #if chrom != '':
        if chrom not in self.data:self.data[chrom]=numpy.array([0.0])
        nplst=numpy.array(lst)#[0:]
        if instep!=self.step:
            lth=len(lst)*(functions.div(instep*1.0,self.step))
            nplst=numpy.array([0.0])
            nplst.resize(int(lth),refcheck=0)
            for pos in range(len(lst)):
                #nplst[pos*instep/self.step]=lst[pos] ########## deleted by kaifu on Sep 6, 2012 ##########
                if instep<self.step:nplst[functions.div(pos*instep,self.step)]+=lst[pos] ########## add by kaifu on Sep 6, 2012 ##########
                else:
                    for tpos in range(functions.div(pos*instep,self.step),functions.div((pos+1)*instep,self.step)):nplst[functions.div(tpos,self.step)]+=functions.div(lst[pos]*1.0*self.step,instep) ########## add by kaifu on Sep 6, 2012 ##########
        self.data[chrom].resize(functions.div(start,self.step)+len(nplst),refcheck=0)
        self.data[chrom][(functions.div(start,self.step)):]=nplst

    def loadVar(self,file,gfile,suppress=False):
        '''
        Description:
            load data from Fixed wiggle format file
        Parameter:
            file:a path to the file containging the data
        Value:
            None
        '''
        
        ########## ---start---add by kaifu on Aug 14, 2012 ##########
        if self.step<1:
            sys.stdout.write('set step size to 10\n')
            self.step=10
        ########## ---end---add by kaifu on Aug 14, 2012 ##########
        sys.stdout.write('parsing from ' + file + "\n")
        step=self.step

        starttime=time()
        chrom,size='',0
        if file[-3:]=='wig':fi=open(file)
        if file[-6:]=='wig.gz':
            import gzip
            fi=gzip.open(file)
        for line in fi:
            #print line
            if line[0]=='t':continue
            elif line[0]=='v':
                #if chrom != '':
                #    if not self.data.has_key(chrom):self.data[chrom]=numpy.array([0.0])#lst
                #    print line[:-1]
                if not suppress: sys.stdout.write(line[:-1] + "\n")
                col=line.split()
                right=0
                for term in col[1:]:
                    kv=term.split("=")
                    if kv[0]=='chrom':
                        right+=1
                        chrom=kv[1]
                        #print chrom,#line[:-1]
                        if chrom not in self.data:
                            self.data[chrom]=numpy.array([0.0])#lst
                        lst=self.data[chrom]#chrom=kv[1]
                        size=lst.size
                    elif kv[0]=='span':
                        inspan=int(kv[1])
                        right+=1
                if right<2:
                    sys.stdout.write('wrong format: ' + line + ' chrom and span must be provided\n')
            else:
                col=line.split()
                tstart,value=int(col[0]),float(col[1])
                tend=functions.div((tstart+inspan),step)
                tstart=functions.div(tstart,step)
                if tend>=size:
                    #while(tend>=size):
                    size=tend+1000
                    lst.resize(size,refcheck=0)
                for pos in range(tstart,tend):lst[pos]=value
        self.clearEmptyEnd()
        #print 'time cost:',time()-starttime
    def maxmin(self,nonzero=False):
        ma,mi=0,0
        for cr in self.data:
            tma,tmi=self.data[cr].max(),self.data[cr].min()
            if nonzero:
                tpos=self.data[cr].nonzero()
                tma,tmi=self.data[cr][tpos].max(),self.data[cr][tpos].min()
                if mi==0:mi=tmi
                if ma==0:ma=tma
            if tma>ma:ma=tma
            if tmi<mi:mi=tmi
        return [ma,mi]
    
    def mean(self):
        '''
        Description:
            Calculate the mean occupancy value
        Parameter:
            None
        Value:
            float value
        '''
        #return self.sum()*1.0/self.gsize()
        size,value=0,0
        for chrom in self.data:
            value+=self.data[chrom].sum()
            size+=self.data[chrom].size
        return functions.div(value,size)
    def multiply(self,wig2):
        '''
        Description:
            multiply by wig2 at each data point
        Parameter:
            wig2: the Wiggle class instance to multiply by
        Value:
            an Wiggle class instance
        '''
        if self.step!=wig2.step:
            wig2=deepcopy(wig2)
            wig2.changeStep(self.step)
        wg=Wig(step=self.step)
        chrs={}
        for chrom in self.data:chrs[chrom]=1
        for chrom in wig2.data:chrs[chrom]=1
        chrs=list(chrs.keys())
        for chrom in chrs:
            if chrom not in self.data:self.addChr(chrom)
            if chrom not in wig2.data:wig2.addChr(chrom)
            lth1=self.data[chrom].size
            lth2=wig2.data[chrom].size
            if lth1>lth2:wig2.data[chrom].resize(lth1,refcheck=0)
            elif lth1<lth2:self.data[chrom].resize(lth2,refcheck=0)
            wg.data[chrom]=self.data[chrom]*wig2.data[chrom]
        return wg
    def non0to1(self):
        '''
        change all non-zero value to value 1
        #can be further improved by using the numpy.nonzero() function
        '''
        for cr in self.data:
            a=self.data[cr]
            a[a.nonzero()]=1

    def non0size(self):
        '''
        the number of data point with a non-zero value.
        '''
        size=0
        for cr in self.data:
            size+=self.data[cr].nonzero()[0].size
        return size*self.step

    def percentile(self, p=[0,25,50,75,100],bnum=100000,nonzero_end=False):
        '''
        bnum: number of histogram bins to be calculated between maximal value and minimal value
        '''
        hs=self.histogram(bnum=bnum,nonzero_end=nonzero_end)
        total=hs[0].sum()
        p=deepcopy(p)
        p.sort()
        plth=len(p)
        i,j=0,0
        count=0
        while j<bnum:
            count+=hs[0][j]
            while functions.div(count*100.0,total)>=p[i] and i<plth:
                p[i]=functions.div(hs[1][j]*(100-p[i]),100.0)+functions.div(hs[1][j+1]*(p[i]),100.0)
                i+=1
                if i>=plth:
                    p.append(110)
                    j=bnum
            j+=1
        p=p[:plth]
        #print p
        return p
    
    def pop(self,k):
        '''
        Description:
            remove a chrosome 
        Parameter:
            k: the name of the chrosome that is to be removed
        Value:
            None
        '''
        return self.data.pop(k)
    def power(self,p):
        '''
        Description:
            Self multiply by p times
        Parameter:
            p: the times to do self multiply
        Value:
            None
        '''
        chrs={}
        for chrom in self.data:chrs[chrom]=1
        chrs=list(chrs.keys())
        for chrom in chrs:
            self.data[chrom]=self.data[chrom]**p
        return True

    def ppois(self,wig2,tchr=''):
        '''
        Description:
            do Poisson test to calculate differential signial for each data point between two Wig class instances.
        
        Parameter:
            wig2: a Wig class instance to be compared to
            tchr: specify a chrosome that is to be compared, leave it to '' if want to do for all chrosomes.
            
        Value:
            A Wig class instance that cantain data for the differential signial
        
        '''
        #pp=r('''function(q,avg){return(ppois(q,avg,lower.tail=FALSE,log=TRUE)/log(10))}''')###### add by kaifu on Oct 18,2011
        wig1=self#deepcopy(self)
        if wig1.step!=wig2.step:wig2.changeStep(wig1.step)
        out=Wig(step=wig1.step)
        donenum=0
        for chrom in wig1.data:
            sys.stdout.write(chrom + "\n")
            if tchr!='' and chrom!=tchr:continue
            if chrom in wig2.data:
                out.data[chrom]=numpy.array(0.0)
                lth1=wig1.data[chrom].size
                lth2=wig2.data[chrom].size
                lth=max(lth1,lth2)
                if lth>lth2:wig2.data[chrom].resize(lth,refcheck=0)
                elif lth>lth1:wig1.data[chrom].resize(lth,refcheck=0)
                out.data[chrom].resize(lth,refcheck=0)
                tstr1=numpy.array(0.0)
                tstr2=numpy.array(0.0)
                tstr1.resize(lth,refcheck=0)
                tstr2.resize(lth,refcheck=0)
                p=0
                while p<lth:
                    if wig1.data[chrom][p]>wig2.data[chrom][p]:
                        if wig1.data[chrom][p]<1:tstr1[p]=1
                        else:tstr1[p]=wig1.data[chrom][p]
                        if wig2.data[chrom][p]<1:tstr2[p]=1
                        else:tstr2[p]=wig2.data[chrom][p]
                    else:
                        if wig1.data[chrom][p]<1:tstr2[p]=1
                        else:tstr2[p]=wig1.data[chrom][p]
                        if wig2.data[chrom][p]<1:tstr1[p]=1
                        else:tstr1[p]=wig2.data[chrom][p]
                    p+=1
                frglen=1000000   # Fragmentation was implemented by Kaifu Chen on Feb 22, 2013. The reason is that R can not allocate vector space larger than 1Gb.
                frags={}
                for end in range(frglen,lth,frglen):frags[end-frglen]=end
                frags[lth-lth%frglen]=lth
                starts=list(frags.keys())
                starts.sort()
                for start in starts :
                    end=frags[start]
                    sys.stdout.write('region: ' + str(start) + '-' + str(end) + "\n")
                    result=functions.ppois(FloatVector(tstr1[start:end]),FloatVector(tstr2[start:end]))/log(10)
                    p=start
                    while p<end:
                        if tstr1[p]==tstr2[p]:out.data[chrom][p]=0
                        elif wig1.data[chrom][p]>wig2.data[chrom][p]:out.data[chrom][p]=0-result[p-start]
                        elif wig1.data[chrom][p]<wig2.data[chrom][p]:out.data[chrom][p]=result[p-start]
                        p+=1
        return out

    def regionWithinValueRange(self, lowValue=None,highValue=None):
        '''
        Each value between lowValue and highValue will be set to 1, other value will be set to 0
        #may be further improved by using the numpy.where() or numpy.nonzero() function
        '''
        twig=deepcopy(self)
        if lowValue==None and highValue==None:
            twig.foldChange(0)
            return twig
        ma,mi=twig.maxmin()
        if lowValue==None:lowValue=mi-1
        if highValue==None:highValue=ma+1
        for cr in twig.data:
            twig.data[cr][numpy.where([((self.data[cr]<lowValue) | (self.data[cr]>highValue))])[1]]=0        
            twig.data[cr][numpy.where([((self.data[cr]>=lowValue) & (self.data[cr]<=highValue))])[1]]=1       
        return twig
    
    
    
    def resizeChr(self,chrom,size):
        '''
        Description:
            change chrosome size
        Parameter:
            chrom: the name of the chrosome whose size is going to be changed
            size: the new size to be setted for the chrosome
        Value:
            None
        '''
        self.data[chrom].resize(functions.div(size,self.step),refcheck=0)
    def rvNeg(self):
        '''
        Description:
            set negative value at each data point to 0
        Parameter:
            None
        Value:
            None
        '''
        for chrom in self.data:
            end=functions.div(self.chrSize(chrom),self.step)
            t=self.data[chrom]
            self.data[chrom]=functions.div(((t**2)**0.5+t),2)


    def save(self,file,format="fixed",step=None,suppress=False):
        '''
        Description:
            Save data to wiggle format file
        Parameter:
            file: a path to the output file
            format: the format of the output wiggle file, could be 'fixed' or 'var'
            step: the step size of the ouput wiggle file
        Value:
            None
        '''
        if step==None:step=self.step
        if not suppress:sys.stdout.write('saving  to '+file + "\n")
        outf=open(file,"w")
        if format==None:#add by Kaifu on Feb 10, 2014
            from random import randint
            chrs=list(self.data.keys())
            lth=len(chrs)-1
            count,total=0,1000
            for i in range(total):
                chrom=chrs[randint(0,lth)]
                val=self.data[chrom][randint(0,sef.data[chrom].size-1)]
                if val!=0:count+=1
            if functions.div(count*1.0,total)>=0.5:format='fixed'
            else:format='var'
            sys.stdout.write('Will save in .'+format+' format\n')
        if format=='wiq':
            tlen=0
            for cr in self.data:tlen+=self.data[cr].size
            from numpy.random import seed,randint
        acculen=0
        for chrom in self.data:
            #if self.data[chrom].sum()==0:continue
            if not suppress:sys.stdout.write(chrom + "\n")
            if format=='fixed':
                outf.write("fixedStep chrom="+chrom+" start=1  step="+str(step)+" span="+str(step)+"\n")
                lth=len(self.data[chrom])
                ot=[]
                for v in self.data[chrom]:ot.append(str(v))#str(self.data[chrom][i]))
                outf.write('\n'.join(ot)+"\n")
            elif format=='var':
                lth=len(self.data[chrom])
                outf.write("variableStep chrom="+chrom+" span="+str(step)+"\n")
                #for i in range(0,lth):
                i=0
                while i<lth:
                    if self.data[chrom][i]!=0:outf.write(str(i*step)+'\t'+str(self.data[chrom][i])+'\n')
                    i+=1
            elif format=='wiq':
                lth=len(self.data[chrom])
                #outf.write("variableStep chrom="+chrom+" span="+str(step)+"\n")
                #for i in range(0,lth):
                i=0
                while i<lth:
                    #if self.data[chrom][i]!=0:outf.write(str(i*step)+'\t'+str(self.data[chrom][i])+'\n')
                    acculen+=1
                    seed(acculen)
                    outf.write(str(self.data[chrom][i])+'\t'+str(randint(0,tlen))+'\t'+chrom+'\t'+str(i*step)+'\n')
                    i+=1
            else:
                sys.stdout.write(',format not recogonized, will be saved in fixed format, ')
                outf.write("fixedStep chrom="+chrom+" start=1  step="+str(step)+" span="+str(step)+"\n")
                lth=len(self.data[chrom])
                for i in range(0,lth):outf.write(str(self.data[chrom][i])+"\n")
        outf.close()
        sys.stdout.write('completed\n')
    def saveChr(self,file,chrom=None,lth=None,format="fixed",step=10):
        '''
        Description:
            Save data to wiggle format file by chrosome name
        Parameter:
            file: a path to the output file
            format: the format of the output wiggle file, could be 'fixed' or 'var'
            step: the step size of the ouput wiggle file
            chrom: the name of the chrosome that is going to be saved
            lth: the length of the chrosome to be saved, start from 1 to lth
        Value:
            None
        '''
        tchr=chrom
        sys.stdout.write('saving wig to' + file + '\n')
        outf=open(file,"w")
        for chrom in self.data:
            if chrom!=tchr:continue
            if format=='fixed':
                outf.write("fixedStep chrom="+chrom+" start=1  step="+str(step)+" span="+str(step)+"\n")
                if lth==None:lth=len(self.data[chrom])*self.step
                for i in range(0,lth,step):outf.write(str(self.data[chrom][functions.div(i,self.step)])+"\n")
        outf.close()
        sys.stdout.write('completed\n')

    def sd(self):
        '''
        Description:
            Calculate standard deviation of occupancy
        Parameter:
            None
        Value:
            None
        '''
        sz=self.size()
        avg=functions.div(self.sum(),sz)
        sqm=0
        for chrom in self.data:
            for v in self.data[chrom]:sqm+=(v-avg)*(v-avg)
        sqm=functions.div(sqm,sz)
        return sqrt(sqm)
    

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
            sizes[col[0]]=int(col[1])
        for chrom in self.data:
            if chrom not in sizes:self.data.pop(chrom)
            else:
                sys.stdout.write(str(self.data[chrom].size) + ' ')
                self.data[chrom].resize(functions.div(sizes[chrom],self.step),refcheck=0)
                sys.stdout.write(str(self.data[chrom].size)  + "\n")
    def smooth(self,lmd=100):
        '''
        Description:
            Smooth occupancy by a sliding window
        Parameter:
            the size of the smooth window
        Value:
            None
        '''
        sys.stdout.write('smooth width: '+str(lmd)+"\n")
        ss=time()
        lmd=int(functions.div(lmd,self.step))
        if lmd<=0:return True
        hlmd=int(functions.div(lmd,2))
        tlmd=lmd-hlmd
        wg2=deepcopy(self.data)
        self.foldChange(0.0)
        wg1=self.data
        for chrom in wg1:
            lth=wg1[chrom].size
            if tlmd!=0-hlmd:
                for p in range(0-hlmd,tlmd):wg1[chrom][hlmd:(lth-tlmd)]+=wg2[chrom][(hlmd+p):(lth-tlmd+p)]
            else:
                for p in range(0-hlmd+1,tlmd):wg1[chrom][hlmd-1:(lth-tlmd)]+=wg2[chrom][(hlmd-1+p):(lth-tlmd+p)]
                for p in [0-hlmd,tlmd]:wg1[chrom][hlmd:(lth-tlmd)]+=wg2[chrom][(hlmd+p):(lth-tlmd+p)]*0.5
            wg1[chrom][hlmd:(lth-tlmd)]= functions.div(wg1[chrom][hlmd:(lth-tlmd)],(lmd)*1.0)
        self.data=wg1
        return True
    def sqrt(self):
        '''
        Description:
            translate each data point to its square root value
        Parameter:
            None
        Value:
            None
        '''
        chrs={}
        for chrom in self.data:chrs[chrom]=1
        chrs=list(chrs.keys())
        for chrom in chrs:
            self.data[chrom]=numpy.sqrt(self.data[chrom])
        return True
    def subtract(self,wig2):
        '''
        Description:
            Subtract wig2 at each data point
        Parameter:
            wig2: the Wig class instance to be subtracted
        Value:
            None
        '''
        if self.step!=wig2.step:
            wig2=deepcopy(wig2)
            wig2.changeStep(self.step)
        chrs={}
        for chrom in self.data:chrs[chrom]=1
        for chrom in wig2.data:chrs[chrom]=1
        chrs=list(chrs.keys())
        for chrom in chrs:
            lth1=len(self.data[chrom])
            lth2=len(wig2.data[chrom])
            if lth1>lth2:wig2.data[chrom].resize(int(lth1),refcheck=0)
            else:self.data[chrom].resize(int(lth2),refcheck=0)
            self.data[chrom]-=wig2.data[chrom]
        wig2.clearEmptyEnd()
        return True
    def sum(self):
        '''
        Description:
            return the sum of occupancy across the whole genome
        Parameter:
            None
        Value:
            None
        '''
        return self.mean()*self.gsize()

    def sumWithinValueRange(self, lowValue,highValue):
        '''
        add values between lowValue and highValue, ignore other values
        '''
        twig=self#deepcopy(self)
        total=0
        for cr in twig.data:
            size=twig.data[cr].size()
            i=0
            while i<size:
                if twig.data[cr][i]>=lowValue and twig.data[cr][i]<=highValue:total+=twig.data[cr][i]
                i+=1        
        return total
    
    def sam_coverage(self,sam_file,step=None):
        if step==None:step=self.step
        else:self.step=step
        infile=open(sam_file)#os.popen('samtools view -XS '+sam_file)
        line=infile.readline()
        hlines=[]
        while line[0]=='@':
            hlines.append(line)
            col=line.split()
            if col[0]=='@SQ':
                chrom,clen=col[1][3:],int(col[2][3:])
                self.data[chrom]=numpy.array([0.0])
                self.data[chrom].resize(functions.div(clen,step),refcheck=0)
            line=infile.readline()
            
        infile=open(sam_file)
        for i in range(len(hlines)):infile.readline()
        for line in infile:
            col=line[:-1].split('\t')
            chrom=col[2]
            t1=re.findall('\d+',col[5])#.split('\d+'))
            t2=re.findall('\D+',col[5])
            start=int(col[3])
            for i in range(len(t2)):
                if t2[i]=='M' :
                    end=start+int(t1[i])
                    self.data[chrom][functions.div(start,step):functions.div(end,step)]+=1
                elif t2[i]=='D':start=start+int(t1[i])
                elif t2[i]=='N':start=start+int(t1[i])
                
                '''
                elif t2[i]=='S':continue
                elif t2[i]=='H':continue
                elif t2[i]=='P':continue
                elif t2[i]=='I':continue
                '''
    def std(self):
        v=0
        m=self.mean()
        gs=0
        for chrom in self.data:
            ts=self.data[chrom].size
            i=0
            while i<ts:
                s=self.data[chrom][i]-m
                v+=s*s
                i+=1
            gs+=ts
        v=sqrt(functions.div(v,gs))
        return v
if __name__ == "__main__":
    import sys,re,os
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # This allow DANPOS to print each message on screen immediately.
