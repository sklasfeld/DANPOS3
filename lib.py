#!/usr/bin/env python
import os,sys,re
from rpy2.robjects import r,FloatVector,StrVector
from glob import glob
from wig import Wig
from wigs import Wigs
from reads import reads
import numpy
from copy import deepcopy
from math import log10
from functions import div

#from functions import merge_peaks_by_head_tail_distance
def batchOccInRegions(wgs,outname=None,groupname='',outmode='w',chrColID=1,nameColID=0,startColIDpos=3,startColIDneg=4,endColIDpos=4,endColIDneg=3,straColID=2,sep='\t',second_sep=None,step=0,\
                   lines=None,heatMap=True,flankup=3000,flankdn=3000,vcal='mean',excludeP=1,region_size=3000):
    '''
    parameters:
        wgs: the Wigs class object
    '''
    #calculate average wiggle density in regions and their flanking regions,e.g., gene 
    #for gene body: occInRegions(wg=wg,chrColID=1,nameColID=0,startColIDpos=3,startColIDneg=4,endColIDpos=4,endColIDneg=3,straColID=2,step=1000,sep='\t',second_sep=None,\
    #             lines=lines,heatmapname=None,avgcurvename=None,flankup=1000000,flankdn=1000000,vcal='mean',excludeP=0.01,bin_count=100)
    #vcal: the method to calculate plot value, could be median or mean
    keys=list(wgs.keys())
    keys.sort()
    if len(keys)<1:
        print('at least one wiggle data need to be specified!\n')
        return
    if step<1:
        steps={}
        for k in keys:steps[wgs.get(k).step]=1
        steps=list(steps.keys())
        if len(steps)>1:
            steps.sort()
            print('step sizes in wiggles are not the same, will be set to a common step size',steps[0])
            step=steps[0]
        else:step=steps[0]
    dic={}
    for k in keys:
        print('\ncalculating for ',k,'...')
        wg=wgs.get(k)
        tHeatMapName=None
        #if heatMap:
        #if groupname!=None:tHeatMapName=os.path.join(os.path.split(groupname)[0],os.path.split(groupname)[-1]+'_'+os.path.split(k)[-1]+'.heatmap')
        if outname!=None:tHeatMapName=outname+'_heatmap'
        else:tHeatMapName='heatmap'
        if not os.path.isdir(tHeatMapName):os.mkdir(tHeatMapName)
        if groupname!=None:tHeatMapName=os.path.join(tHeatMapName,os.path.split(groupname)[-1]+'.'+os.path.split(k)[-1]+'.heatmap')
        else:tHeatMapName=k+'.heatmap'
        dic[k]=occInRegions(wg=wg,chrColID=chrColID,nameColID=nameColID,startColIDpos=startColIDpos,startColIDneg=startColIDneg,endColIDpos=endColIDpos,endColIDneg=endColIDneg,straColID=straColID,step=step,sep=sep,second_sep=second_sep,lines=lines,heatmapname=tHeatMapName,avgcurvename=outname,flankup=flankup,flankdn=flankdn,vcal=vcal,excludeP=excludeP,region_size=region_size)
    #print outmode,'outmode'
    if outname!=None:
        if outmode!='w':
            fi=open(outname+'.xls')
            lines=fi.readlines()
            lines[0]=lines[0][:-1]
            fi.close()
        else:
            lines=['pos']
        #print outname
        #print lines[0]
        for k in keys:lines[0]+='\t'+k+'.'+groupname
        #print lines[0]
        fo=open(outname+'.xls','w')
        fo.write(lines[0]+'\n')
        poses=list(dic[keys[0]].keys())
        poses.sort()
        if outmode!='w' and len(poses)!=len(lines)-1:return {}
        for i in range(len(poses)):
            p=poses[i]
            if outmode=='w':
                if p<region_size:oline=str(p)
                else:oline='+'+str(p-region_size)
            else:oline=lines[i+1][:-1]
            for k in keys:oline+='\t'+str(dic[k][p])
            fo.write(oline+'\n')
        fo.close()
    return dic

def batchOccAroundPoints(wgs,outname=None,groupname='',outmode='w',chrColID=1,nameColID=0,posColIDpos=3,posColIDneg=4,straColID=2,sep='\t',second_sep=None,step=0,lines=None,heatMap=True,flankup=3000,flankdn=3000,vcal='mean',excludeP=1):
    '''
    parameters:
        wgs: the Wigs class object
    '''
    #calculate average wiggle density along the flanking regions of Transcription Start sites, Transcription Terminal Sites, or exon/intron junctions
    #for TSS:chrColID=1,nameColID=0,posColIDpos=3,posColIDneg=4,straColID=2,sep='\t',second_sep=None,
    #for TTS:chrColID=1,nameColID=0,posColIDpos=4,posColIDneg=3,straColID=2,sep='\t',second_sep=None,
    #for CDS_SS:chrColID=1,nameColID=0,posColIDpos=5,posColIDneg=6,straColID=2,sep='\t',second_sep=None,
    #for CDS_TS:chrColID=1,nameColID=0,posColIDpos=6,posColIDneg=5,straColID=2,sep='\t',second_sep=None,
    #for exonSS:chrColID=1,nameColID=0,posColIDpos=8,posColIDneg=9,straColID=2,sep='\t',second_sep=',',
    #for exonTs:chrColID=1,nameColID=0,posColIDpos=9,posColIDneg=8,straColID=2,sep='\t',second_sep=',',
    keys=list(wgs.keys())
    keys.sort()
    if len(keys)<1:
        print('at least one wiggle data need to be specified!\n')
        return
    if step<1:
        steps={}
        for k in keys:steps[wgs.get(k).step]=1
        steps=list(steps.keys())
        if len(steps)>1:
            steps.sort()
            print('step sizes in wiggles are not the same, will be set to a common step size',steps[0])
            step=steps[0]
        else:step=steps[0]
    dic={}
    for k in keys:
        print('\ncalculating for ',k,'...')
        wg=wgs.get(k)
        tHeatMapName=None
        #if heatMap:
        #if groupname!=None:tHeatMapName=os.path.join(os.path.split(groupname)[0],os.path.split(groupname)[-1]+'_'+os.path.split(k)[-1]+'.heatmap')
        if outname!=None:tHeatMapName=outname+'_heatmap'
        else:tHeatMapName='heatmap'
        if not os.path.isdir(tHeatMapName):os.mkdir(tHeatMapName)
        if groupname!=None:tHeatMapName=os.path.join(tHeatMapName,os.path.split(groupname)[-1]+'.'+os.path.split(k)[-1]+'.heatmap')
        else:tHeatMapName=k+'.heatmap'
        dic[k]=occAroundPoints(wg=wg,chrColID=chrColID,nameColID=nameColID,posColIDpos=posColIDpos,posColIDneg=posColIDneg,straColID=straColID,step=step,sep=sep,second_sep=second_sep,lines=lines,heatmapname=tHeatMapName,avgcurvename=None,flankup=flankup,flankdn=flankdn,vcal=vcal,excludeP=excludeP)
        #dic[k]=occAroundPoints(wg=wg,chrColID=chrColID,nameColID=nameColID,posColIDpos=posColIDpos,posColIDneg=posColIDneg,straColID=straColID,step=step,sep=sep,second_sep=second_sep,lines=lines,heatmapname=tHeatMapName,avgcurvename=outname,flankup=flankup,flankdn=flankdn,vcal=vcal,excludeP=excludeP)
    if outname!=None:
        if outmode!='w':
            fi=open(outname+'.xls')
            lines=fi.readlines()
            lines[0]=lines[0][:-1]
            fi.close()
        else:
            lines=['pos']
        for k in keys:lines[0]+='\t'+k+'.'+groupname
        fo=open(outname+'.xls','w')
        fo.write(lines[0]+'\n')
        poses=list(dic[keys[0]].keys())
        poses.sort()
        if outmode!='w' and len(poses)!=len(lines)-1:return {}
        for i in range(len(poses)):
            p=poses[i]
            if outmode=='w':oline=str(p)
            else:oline=lines[i+1][:-1]
            for k in keys:oline+='\t'+str(dic[k][p])
            fo.write(oline+'\n')
        fo.close()
    return dic


def batchOccPSD(wgs,outname=None):
    keys=list(wgs.keys())
    if len(keys)<1:
        print('at least one wiggle data need to be specified!\n')
        return
    steps={}
    for k in keys:
        print('calculating for',k,'...')
        steps[wgs.get(k).step]=1
    steps=list(steps.keys())
    if len(steps)>1:
        steps.sort()
        print('step sizes in wiggles are not the same, will be set to a common step size',steps[0])
        for k in keys:wgs.get(k).changeStep(step=steps[0])
    dic={}
    for k in keys:dic[k]=occPSD(wgs.get(k),outname=None)
    if outname!=None:
        fo=open(outname+'.xls','w')
        fo.write('Periodicity\t'+'\t'.join(keys)+'\n')
        poses=list(dic[keys[0]].keys())
        poses.sort()
        for i in poses:
            oline=str(i)
            for k in keys:oline+='\t'+str(dic[k][i])
            fo.write(oline+'\n')
    return dic
def batchPositionDistanceDistribution(data,outname=None,min=100,max=250,step=1):#={'name':[]},outname=None,crColID=0,posColID=3,min=100,max=250,step=1):
    keys=list(data.keys())
    dic={}
    for k in keys:dic[k]=positionDistanceDistribution(dic=data[k],min=min,max=max,step=step)#lines=dic[k],outname=None,crColID=crColID,posColID=posColID,min=min,max=max,step=step)
    if outname!=None:
        fo=open(outname+'.xls','w')
        fo.write('Distance\t'+'\t'.join(keys)+'\n')
        poses=list(dic[keys[0]].keys())
        poses.sort()
        for i in poses:
            oline=str(i)
            for k in keys:oline+='\t'+str(dic[k][i])
            fo.write(oline+'\n')
    return dic
def batchPositionValDistribution(data,outname=None,min=0,max=1500,step=3):
    keys=list(data.keys())
    dic={}
    for k in keys:dic[k]=positionValDistribution(dic=data[k],outname=None,min=min,max=max,step=step)
    if outname!=None:
        fo=open(outname+'.xls','w')
        fo.write('Value\t'+'\t'.join(keys)+'\n')
        poses={}
        for name in dic:#[keys[0]].keys()
            for pos in dic[name]:poses[pos]=0
        poses=list(poses.keys())
        poses.sort()
        for i in poses:
            oline=str(i)
            for k in keys:
                if i not in dic[k]:dic[k][i]=0
                oline+='\t'+str(dic[k][i])
            fo.write(oline+'\n')
    return dic
def batchPositionAroundPoints(smts,outname=None,flankup=2500,flankdn=2500,step=10,chrColID=1,nameColID=0,posColIDpos=3,posColIDneg=4,straColID=2,sep='\t',\
                          second_sep=None,lines=[]):
    dic={}
    keys=list(smts.keys())
    keys.sort()
    for k in smts:dic[k]=positionAroundPoint(smts[k],outname=outname,flankup=flankup,flankdn=flankdn,step=step,chrColID=chrColID,nameColID=nameColID,\
                                         posColIDpos=posColIDpos,posColIDneg=posColIDneg,straColID=straColID,sep=sep,second_sep=second_sep,\
                                         lines=[])
    if outname!=None:
        fo=open(outname+'.xls','w')
        fo.write('pos\t'+'\t'.join(keys)+'\n')
        poses=list(dic[keys[0]].keys())
        poses.sort()
        for i in poses:
            oline=str(i)
            for k in keys:oline+='\t'+str(dic[k][i])
            fo.write(oline+'\n')
    return dic

def occAroundPoints(wg,chrColID,nameColID,posColIDpos,posColIDneg,straColID,step=0,sep='\t',second_sep=None,\
                    lines=[],heatmapname=None,avgcurvename=None,flankup=3000,flankdn=3000,vcal='mean',excludeP=0):
    #calculate average wiggle density along the flanking regions of Transcription Start sites, Transcription Terminal Sites, or exon/intron junctions
    #for TSS:chrColID=1,nameColID=0,posColIDpos=3,posColIDneg=4,straColID=2,sep='\t',second_sep=None,
    #for TTS:chrColID=1,nameColID=0,posColIDpos=4,posColIDneg=3,straColID=2,sep='\t',second_sep=None,
    #for CDS_SS:chrColID=1,nameColID=0,posColIDpos=5,posColIDneg=6,straColID=2,sep='\t',second_sep=None,
    #for CDS_TS:chrColID=1,nameColID=0,posColIDpos=6,posColIDneg=5,straColID=2,sep='\t',second_sep=None,
    #for exonSS:chrColID=1,nameColID=0,posColIDpos=8,posColIDneg=9,straColID=2,sep='\t',second_sep=',',
    #for exonTs:chrColID=1,nameColID=0,posColIDpos=9,posColIDneg=8,straColID=2,sep='\t',second_sep=',',
    #vcal: the method to calculate plot value, could be median or mean
    if step<1:step=wg.step
    else:wg.changeStep(step=step)
    flankup=div(flankup,step)
    flankdn=div(flankdn,step)
    if heatmapname==None:heatmapname='heatmapname'
    outf=open(heatmapname+'.xls',"w")
    outf.write('name\tmax\tmin\tsum')
    for i in range(0-flankup,flankdn):outf.write('\t'+str(i))
    outf.write('\n')
    #if avgcurvename==None:avgcurvename='avgcurvename'
    #outf2=open(avgcurvename+'.xls',"w")
    #outf2.write('pos\tvalue\n')
        
    #lst=resize([0.0],flankup+flankdn)
    lst={}
    for i in range(flankup+flankdn):lst[i]=[]
    num=0
    for line in lines:
        if line[:-1]=='\n':line=line[:-1]
        col=line.split(sep)
        chr,name,stra=col[chrColID],col[nameColID],col[straColID]
        if chr not in wg.data:continue
        if stra=='+':
            if second_sep!=None:poses=col[posColIDpos].split(second_sep)
            else:poses=[col[posColIDpos]]
            #print poses,line
        elif stra=='-':
            if second_sep!=None:poses=col[posColIDneg].split(second_sep)
            else:poses=[col[posColIDneg]]
        for pos in poses:
            if pos=='':continue
            tlst=[0.0]*(flankup+flankdn)
            tss=div(int(pos),step)
            if stra=='+':
                for i in range(0-flankup,flankdn):
                    try:tlst[flankup+i]=wg.data[chr][tss+i]
                    except:continue
            else:
                for i in range(1-flankdn,flankup+1):
                    try:tlst[flankup-i]=wg.data[chr][tss+i]
                    except:continue
            regionmax,regionmin,regionsum=max(tlst),min(tlst),sum(tlst)
            ol=name+'\t'+'\t'.join([str(regionmax),str(regionmin),str(regionsum)])
            for i in range(flankup+flankdn):ol+='\t'+str(tlst[i])
            outf.write(ol+'\n')
            num+=1
    outf.close()
    print('caculating average curve ...')
    dic={}
    vec=numpy.array([0.0])
    vec.resize(num,refcheck=0)
    for i in range(4,flankup+flankdn+4):
        vec=vec*0
        fi=open(heatmapname+'.xls')
        fi.readline()
        ln=0
        for line in fi:
            vec[ln]=float(line.split()[i])
            ln+=1
        fi.close()
        vec.sort()
        if vcal=='mean':
            s=sum(vec[int(num*excludeP):int(num-num*excludeP)])*1.0
            v=div(s,(num-num*excludeP*2))#len(vec)
        elif vcal=='median':v=vec[div(num,2)]
        #outf2.write(str((i-flankup-4)*step)+'\t'+str(v)+'\n')
        dic[(i-flankup-4)*step]=v
    print('')
    return dic

def occInRegions(wg,chrColID,nameColID,startColIDpos,startColIDneg,endColIDpos,endColIDneg,straColID,step=0,sep='\t',second_sep=None,\
                    lines=[],heatmapname=None,avgcurvename=None,flankup=3000,flankdn=3000,vcal='mean',excludeP=0,region_size=3000):
    #calculate average wiggle density in regions and their flanking regions,e.g., gene 
    #for gene body: occInRegions(wg=wg,chrColID=1,nameColID=0,startColIDpos=3,startColIDneg=4,endColIDpos=4,endColIDneg=3,straColID=2,step=1000,sep='\t',second_sep=None,\
    #             lines=lines,heatmapname=None,avgcurvename=None,flankup=1000000,flankdn=1000000,vcal='mean',excludeP=0.01,bin_count=100)
    #vcal: the method to calculate plot value, could be median or mean
    #print flankup,flankdn
    owg=deepcopy(wg)
    ostep=owg.step
    if step<1:step=wg.step
    else:wg.changeStep(step=step)
    flankup=div(flankup,step)
    flankdn=div(flankdn,step)
    bin_count=div(region_size,step)
    if heatmapname==None:heatmapname='heatmap'
    outf=open(heatmapname+'.xls',"w")
    outf.write('name\tmax\tmin\tsum')
    #if bin_count<max(flankup,flankdn):bin_count=max(flankup,flankdn)
    for i in range(0-flankup,0):outf.write('\t'+str(i))
    for i in range(0,bin_count):outf.write('\t'+str(i))
    for i in range(0,flankdn):outf.write('\t+'+str(i+bin_count))
    outf.write('\n')
    if avgcurvename==None:avgcurvename='avgcurvename'
    #outf2=open(avgcurvename+'.xls',"w")
    #outf2.write('pos\tvalue\n')
        
    #lst={}
    #for i in range(flankup+flankdn+bin_count):lst[i]=[]
    
    num=0
    tlst=numpy.array([0.0])
    tlst.resize(flankup+flankdn+bin_count,refcheck=0)
    for line in lines:
        if line[:-1]=='\n':line=line[:-1]
        col=line.split(sep)
        chr,name,stra=col[chrColID],col[nameColID],col[straColID]
        if chr not in wg.data:continue
        if stra=='+':
            if second_sep!=None:starts,ends=col[startColIDpos].split(second_sep),col[endColIDpos].split(second_sep)
            else:starts,ends=[col[startColIDpos]],[col[endColIDpos]]
        elif stra=='-':
            if second_sep!=None:starts,ends=col[startColIDneg].split(second_sep),col[endColIDneg].split(second_sep)
            else:starts,ends=[col[startColIDneg]],[col[endColIDneg]]
        id,lth=0,len(starts)
        if lth!=len(ends):continue
        while id<lth:
            tlst=tlst*0
            if starts[id]=='':continue
            tss,tes,otss,otes=div(int(starts[id]),step),div(int(ends[id]),step),div(int(starts[id]),ostep),div(int(ends[id]),ostep)
            if stra=='+':
                for i in range(0-flankup,0):
                    try:tlst[flankup+i]=wg.data[chr][tss+i]
                    except:continue
                bstep=div((otes-otss)*1.0,bin_count)
                for i in range(0,bin_count):
                    try:tlst[flankup+i]=owg.data[chr][otss+int(i*bstep)]
                    except:continue
                for i in range(0,flankdn):
                    try:tlst[flankup+bin_count+i]=wg.data[chr][tes+i]
                    except:continue
            else:
                for i in range(0-flankup,0):
                    try:tlst[flankup+i]=wg.data[chr][tss-i]
                    except:continue
                bstep=div((otss-otes)*1.0,bin_count)
                for i in range(0,bin_count):
                    try:tlst[flankup+i]=owg.data[chr][otss-int(i*bstep)]
                    except:continue
                for i in range(0,flankdn):
                    try:tlst[flankup+bin_count+i]=wg.data[chr][tes-i]
                    except:continue
            regionmax,regionmin,regionsum=max(tlst),min(tlst),sum(tlst)
            #if heatmapname!=None:
            ol=name+'\t'+'\t'.join([str(regionmax),str(regionmin),str(regionsum)])
            for i in range(flankup+flankdn+bin_count):ol+='\t'+str(tlst[i])
            outf.write(ol+'\n')
            #outf.write(ol+'\t'+str(otes-otss)+'\t'+str(bstep)+'\n')
            num+=1
            id+=1
    outf.close()
    
    print('caculating average curve ...')
    dic={}
    vec=numpy.array([0.0])
    vec.resize(num,refcheck=0)
    #print flankup,flankdn,bin_count,flankup+flankdn+bin_count
    for i in range(4,flankup+flankdn+bin_count+4):
        vec=vec*0
        fi=open(heatmapname+'.xls')
        fi.readline()
        ln=0
        for line in fi:
            vec[ln]=float(line.split()[i])
            ln+=1
        fi.close()
        vec.sort()
        if vcal=='mean':
            s=sum(vec[int(num*excludeP):int(num-num*excludeP)])*1.0
            v=div(s,(num-num*excludeP*2))#len(vec)
        elif vcal=='median':v=vec[div(num,2)]
        #outf2.write(str((i-flankup-4)*step)+'\t'+str(v)+'\n')
        dic[(i-flankup-4)*step]=v
    print('')
    return dic

def positionSelectorOld(positionLines=[],selection=None,geneFile=None,outGeneFile=None,chrbinsize=10000000000000):
    '''
    positionLines:  The first line in positionLines must be the title line, each line should have a '\n' at the end, positionLines should be in the default output format of DANPOS, see DANPOS documentation for example
    selection:  promoter:-350:50&control_smt_val:0:1000&0-log10(point_diff_pval):0:1e-10
    '''
    if selection==None:return positionLines
    if len(positionLines)<2:return positionLines
    if outGeneFile!=None:
        ogf=open(outGeneFile,'w')
        ogf.write('genes\t'+positionLines[0])
        oglines={}
    retr,tcol,tld=[positionLines[0]],positionLines[0].split('\t'),{}
    for i in range(len(tcol)):tld[tcol[i]]=i
    sels=selection.split('&')
    flank=0
    for i in range(len(sels)):
        sels[i]=sels[i].split(':')
        if sels[i][1]!='':
            sels[i][1]=float(sels[i][1])
            if sels[i][0]=='promoter':
                if geneFile==None:
                    print('Error: gene file is not provided.')
                    return []
                if abs(sels[i][1])>flank:flank=abs(sels[i][1])
        if sels[i][2]!='':
            sels[i][2]=float(sels[i][2])
            if sels[i][0]=='promoter':
                if geneFile==None:
                    print('Error: gene file is not provided.')
                    return []
                if abs(sels[i][2])>flank:flank=abs(sels[i][2])
    if flank>chrbinsize:chrbinsize=flank
    if geneFile!=None:
        gd={}
        for line in open(geneFile).readlines()[1:]:
            col=line.split()
            if col[2]=='+':gname,cr,stra,tss=col[0],col[1],col[2],int(col[3])
            else:gname,cr,stra,tss=col[0],col[1],col[2],int(col[4])
            if cr not in gd:gd[cr]={}
            bin=int(div(tss,chrbinsize))
            if bin not in gd[cr]:gd[cr][bin]={}
            gd[cr][bin][tss]=[stra,gname]
    for line in positionLines[1:]:
        genes,ok,col='',True,line.split('\t')
        for sel in sels:
            if sel[0]!='promoter':
                v=float(col[tld[sel[0]]])
                if sel[1]!='':
                    if v<sel[1]:ok=False
                if sel[2]!='':
                    if v>sel[2]:ok=False
        for sel in sels:
            if ok and sel[0]=='promoter':
                cr,poses=col[0],[]
                for i in ([1,2,3]):
                    try:poses.append(int(col[i]))
                    except:continue
                minpos,maxpos=min(poses),max(poses)
                minbin,maxbin=int(div((minpos-flank),chrbinsize)),int(div((maxpos+flank),chrbinsize))
                bins=list(range(minbin,maxbin+1))
                for bin in bins:
                    if bin not in gd[cr]:continue
                    for tss in gd[cr][bin]:
                        if gd[cr][bin][tss][0]=='+' and (maxpos>tss+sel[1] and minpos<tss+sel[2]):
                            if outGeneFile!=None:
                                if gd[cr][bin][tss][1] not in oglines:oglines[gd[cr][bin][tss][1]]=[]
                                oglines[gd[cr][bin][tss][1]].append(gd[cr][bin][tss][1]+'\t'+line)
                            genes+=gd[cr][bin][tss][1]+','
                        elif gd[cr][bin][tss][0]=='-' and (maxpos>tss-sel[2] and minpos<tss-sel[1]):
                            if outGeneFile!=None:
                                if gd[cr][bin][tss][1] not in oglines:oglines[gd[cr][bin][tss][1]]=[]
                                oglines[gd[cr][bin][tss][1]].append(gd[cr][bin][tss][1]+'\t'+line)
                            genes+=gd[cr][bin][tss][1]+','
                if genes=='':ok=False
        if ok:retr.append(line)
    if outGeneFile!=None:
        for gene in oglines:
            for line in oglines[gene]:ogf.write(line)
        ogf.close()
    return retr

def GREATdomain(bup=-5000,bdn=1000,eup=-1000000,edn=1000000,geneFile=None,posStartCol=3,posEndCol=3,negStartCol=4,negEndCol=4,chrbinsize=1000000):
    '''
    #this function is not finished
    if negStartCol==None:negStartCol=posEndCol
    if negEndCol==None:negEndCol=posStartCol
    tgd,gd,n={},{},0
    #define basal and extenede domains for each element
    for line in open(geneFile).readlines()[1:]:
        n,col=n+1,line.split()
        gname,cr,stra=col[0],col[1],col[2]
        if not tgd.has_key(cr):tgd[cr]={}
        if stra=='+':
            starts,ends=col[posStartCol].split(','),col[posEndCol].split(',')
            if starts[-1]=='':starts,ends=starts[:-1],ends[:-1]
            for i in range(len(starts)):
                start,end=int(starts[i]),int(ends[i])#biological start and end
                bstart,bend,estart,eend=start+bup,end+bdn,start+bup+eup,end+bdn+edn#start and end of basal and extended domain
                if not tgd[cr].has_key(bstart):tgd[cr][bstart]=[]
                tgd[cr][bstart].append([bend,gname,start,end,estart,eend])
            
        else:
            starts,ends=col[negStartCol].split(','),col[negEndCol].split(',') 
            if starts[-1]=='':starts,ends=starts[:-1],ends[:-1]
            for i in range(len(starts)):
                start,end=int(starts[i]),int(ends[i])#biological start and end
                bstart,bend,estart,eend=end-bdn,start-bup,end-bdn-edn,start-bup-eup#start and end of basal and extended domain, note that the physical starts is the biological ends on the negative strand
                if not tgd[cr].has_key(bstart):tgd[cr][bstart]=[]
                tgd[cr][bstart].append([bend,gname,start,end,estart,eend])
    
    bends={}
    for cr in tgd:
        bends[cr]={}
        for s in tgd[cr]:
            e=0
            for i in range(len(tgd[cr][s])):
                if e<tgd[cr][s][i][0]:e=tgd[cr][s][i][0]
            bends[cr][s]=e

    #define final domain for each element
    for cr in tgd:
        if not gd.has_key(cr):gd[cr]={}
        poses=tgd[cr].keys()
        poses.sort()
        i,lth=1,len(poses)-1
        while i <lth:
            bstart,prebstart,nxtbstart=poses[i],poses[i-1],poses[i+1]#basal start
            #print tgd[cr][bstart],len(tgd[cr][bstart])
            #print bends[cr]
            for j in range(len(tgd[cr][bstart])):
                #print j
                #bend,prebend,nxtbend=tgd[cr][bstart][j][0],bends[prebstart],bends[nxtbstart]#basal end
                bend=tgd[cr][bstart][j][0]#basal end
                prebend,nxtbend=bends[cr][prebstart],bends[cr][nxtbstart]#basal end
                gname,start,end,estart,eend=tgd[cr][bstart][j][1],tgd[cr][bstart][j][2],tgd[cr][bstart][j][3],tgd[cr][bstart][j][4],tgd[cr][bstart][j][5]
                minpos,maxpos=min(bstart,max(estart,prebend)),max(bend,min(eend,nxtbstart))
                if minpos<0:minpos=0
                for bin in range(minpos/chrbinsize,maxpos/chrbinsize+1):
                    if not gd[cr].has_key(bin):gd[cr][bin]={}
                    if not gd[cr][bin].has_key(minpos):gd[cr][bin][minpos]=[]
                    gd[cr][bin][minpos].append([maxpos,gname,start,end])
            i+=1
        print cr,poses[0],tgd[cr][poses[0]]
        
        
        if lth>=1:
            #the first gene on chr
            i=0
            bstart,prebstart,nxtbstart=poses[i],0,poses[i+1]#basal start
            for j in range(len(tgd[cr][bstart])):
                bend,prebend,nxtbend=tgd[cr][bstart][j][0],0,bends[cr][nxtbstart]#basal end
                gname,start,end,estart,eend=tgd[cr][bstart][j][1],tgd[cr][bstart][j][2],tgd[cr][bstart][j][3],tgd[cr][bstart][j][4],tgd[cr][bstart][j][5]
                minpos,maxpos=min(bstart,max(estart,prebend)),max(bend,min(eend,nxtbstart))
                if minpos<0:minpos=0
                for bin in range(minpos/chrbinsize,maxpos/chrbinsize+1):
                    if not gd[cr].has_key(bin):gd[cr][bin]={}
                    if not gd[cr][bin].has_key(minpos):gd[cr][bin][minpos]=[]
                    gd[cr][bin][minpos].append([maxpos,gname,start,end])
                
            #the last gene on chr
            i=lth
            bstart,prebstart,nxtbstart=poses[i],poses[i-1],None#basal start
            for j in range(len(tgd[cr][bstart])):
                bend,prebend,nxtbend=tgd[cr][bstart][j][0],bends[cr][prebstart],None#basal end
                gname,start,end,estart,eend=tgd[cr][bstart][j][1],tgd[cr][bstart][j][2],tgd[cr][bstart][j][3],tgd[cr][bstart][j][4],tgd[cr][bstart][j][5]
                minpos,maxpos=min(bstart,max(estart,prebend)),eend#max(bend,min(eend,nxtbstart))
                if minpos<0:minpos=0
                for bin in range(minpos/chrbinsize,maxpos/chrbinsize+1):
                    if not gd[cr].has_key(bin):gd[cr][bin]={}
                    if not gd[cr][bin].has_key(minpos):gd[cr][bin][minpos]=[]
                    gd[cr][bin][minpos].append([maxpos,gname,start,end])
        elif lth>=0:#only one gene on chr
            i=0
            bstart,prebstart,nxtbstart=poses[i],None,None#basal start
            for j in range(len(tgd[cr][bstart])):
                bend,prebend,nxtbend=tgd[cr][bstart][j][0],None,None#basal end
                gname,start,end,estart,eend=tgd[cr][bstart][j][1],tgd[cr][bstart][j][2],tgd[cr][bstart][j][3],tgd[cr][bstart][j][4],tgd[cr][bstart][j][5]
                minpos,maxpos=estart,eend#min(bstart,max(estart,prebend)),max(bend,min(eend,nxtbstart))
                if minpos<0:minpos=0
                for bin in range(minpos/chrbinsize,maxpos/chrbinsize+1):
                    if not gd[cr].has_key(bin):gd[cr][bin]={}
                    if not gd[cr][bin].has_key(minpos):gd[cr][bin][minpos]=[]
                    gd[cr][bin][minpos].append([maxpos,gname,start,end])
    return gd
    '''
def positionSelectorByGreatTSS(positionLines=[],selection='-5000:1000:1000000',geneFile=None,chrbinsize=None):
    '''
    positionLines:  The first line in positionLines must be the title line, each line should have a '\n' at the end, positionLines should be in the default output format of DANPOS, see DANPOS documentation for example
    '''
    if geneFile==None:
        print('Error: gene file is not provided.')
        return []
    if selection==None:return positionLines
    sels=selection.split(':')
    if len(sels)<3:
        print('Wrong! Less than three fields could be detected in the selector:',selection)
        return []
    elif sels[0]=='':
        print('Wrong! Please set a upstream bound in the GREAT selector:',selection)
        return []
    elif sels[1]=='':
        print('Wrong! Please set a downstream bound in the GREAT selector:',selection)
        return []
    else:sels[0],sels[1]=int(sels[0]),int(sels[1])
    if sels[2]=='':sels[2]=0
    else:sels[2]=int(sels[2])
    if chrbinsize==None:chrbinsize=max(sels)
    if len(positionLines)<2:return positionLines
    if positionLines[0][-1]=='\n':positionLines[0]=positionLines[0][:-1]
    retr,tcol,tld=[positionLines[0]+'\trelatedGenes\n'],positionLines[0].split('\t'),{}
    tgd,gd,n={},{},0
    for line in open(geneFile).readlines()[1:]:
        n,col=n+1,line.split()
        gname,cr,stra=col[0],col[1],col[2]
        if cr not in tgd:tgd[cr]={}
        if stra=='+':
            pos=int(col[3])
            tgd[cr][pos+sels[0]]=[gname,pos,pos+sels[1]]
        else:
            pos=int(col[4])
            tgd[cr][pos-sels[1]]=[gname,pos,pos-sels[0]]
    for cr in tgd:
        if cr not in gd:gd[cr]={}
        if 'TSS' not in gd[cr]:gd[cr]['TSS']={}
        poses=list(tgd[cr].keys())
        poses.sort()
        i,lth=1,len(poses)-1
        while i <lth:
            bstart,prebstart,nxtbstart=poses[i],poses[i-1],poses[i+1]
            bend,prebend,nxtbend=tgd[cr][bstart][2],tgd[cr][prebstart][2],tgd[cr][nxtbstart][2]
            pos,gname=tgd[cr][bstart][1],tgd[cr][bstart][0]
            minpos,maxpos=min(bstart,max(prebend,pos-sels[2])),max(bend,min(nxtbstart,pos+sels[2]))
            for bin in range(div(minpos,chrbinsize),div(maxpos,chrbinsize)+1):
                if bin not in gd[cr]['TSS']:gd[cr]['TSS'][bin]={}
                if minpos not in gd[cr]['TSS'][bin]:gd[cr]['TSS'][bin][minpos]=[]
                gd[cr]['TSS'][bin][minpos].append(['',gname,maxpos,pos])
            i+=1
    for line in positionLines[1:]:
        genes,ok,col='',0,line.split('\t')
        tgenes={}
        cr,poses=col[0],[]
        for i in ([1,2,3]):
            try:poses.append(int(col[i]))
            except:continue
        minpos,maxpos=min(poses),max(poses)
        minbin,maxbin=div(minpos,chrbinsize),div(maxpos,chrbinsize)
        bins=list(range(minbin,maxbin+1))
        for bin in bins:
            if cr not in gd:continue
            if bin not in gd[cr]['TSS']:continue
            for start in gd[cr]['TSS'][bin]:
                for i in range(len(gd[cr]['TSS'][bin][start])):
                    gname,end,pos=gd[cr]['TSS'][bin][start][i][1],gd[cr]['TSS'][bin][start][i][2],gd[cr]['TSS'][bin][start][i][3]
                    if maxpos>start and minpos<end:tgenes[gname+'/'+str(pos)+'/'+str(pos)]=1
        if len(tgenes)!=0:
            ok+=1
            genes+='TSS'+':'+','.join(list(tgenes.keys()))+'|'
        if ok>0:retr.append(line[:-1]+'\t'+genes[:-1]+'\n')
    return retr

def positionSelectorByGeneStructure(positionLines=[],selection=None,geneFile=None,chrbinsize=10000000000000):
    '''
    positionLines:  The first line in positionLines must be the title line, each line should have a '\n' at the end, positionLines should be in the default output format of DANPOS, see DANPOS documentation for example
    selection:  promoter:-350:50&control_smt_val:0:1000&0-log10(point_diff_pval):0:1e-10
    '''
    if geneFile==None:
        print('Error: gene file is not provided.')
        return []
    if selection==None:return positionLines
    if len(positionLines)<2:return positionLines
    '''
    if outGeneFile!=None:
        ogf=open(outGeneFile,'w')
        ogf.write('genes\t'+positionLines[0])
        oglines={}
    '''
    if positionLines[0][-1]=='\n':positionLines[0]=positionLines[0][:-1]
    retr,tcol,tld=[positionLines[0]+'\trelatedGenes\n'],positionLines[0].split('\t'),{}
    sels=selection.split(',')
    #print sels
    if sels[-1]=='and':andor='and'
    elif sels[-1]=='or':andor='or'
    elif len(sels)>1:
        print("Error: the selection must be defined with 'and' or 'or' at the end")
        return []
    else:
        sels.append('and')
        andor='and'
    sels=sels[:-1]
    selsdic,flank,tsels={},0,[]
    for i in range(len(sels)):
        sels[i]=sels[i].split(':')
        if not sels[i][0] in ['TSS','TTS','CSS','CTS','ESS','ETS','exon','intron','gene']:
            print('Error: can not do selection for',sels[i][0])
            return []
        selsdic[sels[i][0]]=1
        if sels[i][1]!='':
            sels[i][1]=int(sels[i][1])
            if abs(sels[i][1])>flank:flank=abs(sels[i][1])
        if sels[i][2]!='':
            sels[i][2]=int(sels[i][2])
            if abs(sels[i][2])>flank:flank=abs(sels[i][2])
        tsels.append(sels[i])
    sels=tsels
    if flank>chrbinsize:chrbinsize=flank
    gd={}
    n=0
    for line in open(geneFile).readlines()[1:]:
        n+=1
        col=line.split()
        cr,stra=col[1],col[2]
        if cr not in gd:gd[cr]={}
        if 'TSS' in selsdic:
            if 'TSS' not in gd[cr]:gd[cr]['TSS']={}
            if stra=='+':gname,cr,stra,pos=col[0],cr,stra,int(col[3])
            else:gname,cr,stra,pos=col[0],cr,stra,int(col[4])
            bin=int(div(pos,chrbinsize))
            if bin not in gd[cr]['TSS']:gd[cr]['TSS'][bin]={}
            if pos not in gd[cr]['TSS'][bin]:gd[cr]['TSS'][bin][pos]=[]
            gd[cr]['TSS'][bin][pos].append([stra,gname,pos])
        
        if 'TTS' in selsdic:
            if 'TTS' not in gd[cr]:gd[cr]['TTS']={}
            if stra=='+':gname,cr,stra,pos=col[0],cr,stra,int(col[4])
            else:gname,cr,stra,pos=col[0],cr,stra,int(col[3])
            bin=int(div(pos,chrbinsize))
            if bin not in gd[cr]['TTS']:gd[cr]['TTS'][bin]={}
            if pos not in gd[cr]['TTS'][bin]:gd[cr]['TTS'][bin][pos]=[]
            gd[cr]['TTS'][bin][pos].append([stra,gname,pos])
        if 'CSS' in selsdic:
            if 'CSS' not in gd[cr]:gd[cr]['CSS']={}
            if stra=='+':gname,cr,stra,pos=col[0],cr,stra,int(col[5])
            else:gname,cr,stra,pos=col[0],cr,stra,int(col[6])
            bin=int(div(pos,chrbinsize))
            if bin not in gd[cr]['CSS']:gd[cr]['CSS'][bin]={}
            if pos not in gd[cr]['CSS'][bin]:gd[cr]['CSS'][bin][pos]=[]
            gd[cr]['CSS'][bin][pos].append([stra,gname,pos])
        if 'CTS' in selsdic:
            if 'CTS' not in gd[cr]:gd[cr]['CTS']={}
            if stra=='+':gname,cr,stra,pos=col[0],cr,stra,int(col[6])
            else:gname,cr,stra,pos=col[0],cr,stra,int(col[5])
            bin=int(div(pos,chrbinsize))
            if bin not in gd[cr]['CTS']:gd[cr]['CTS'][bin]={}
            if pos not in gd[cr]['CTS'][bin]:gd[cr]['CTS'][bin][pos]=[]
            gd[cr]['CTS'][bin][pos].append([stra,gname,pos])
        if 'ESS' in selsdic:
            if 'ESS' not in gd[cr]:gd[cr]['ESS']={}
            if stra=='+':gname,cr,stra,poss=col[0],cr,stra,col[8][:-1].split(',')
            else:gname,cr,stra,poss=col[0],cr,stra,col[9][:-1].split(',')
            for pos in poss:
                pos=int(pos)
                bin=int(div(pos,chrbinsize))
                if bin not in gd[cr]['ESS']:gd[cr]['ESS'][bin]={}
                if pos not in gd[cr]['ESS'][bin]:gd[cr]['ESS'][bin][pos]=[]
                gd[cr]['ESS'][bin][pos].append([stra,gname,pos])
        if 'ETS' in selsdic:
            if 'ETS' not in gd[cr]:gd[cr]['ETS']={}
            if stra=='+':gname,cr,stra,poss=col[0],cr,stra,col[9][:-1].split(',')
            else:gname,cr,stra,poss=col[0],cr,stra,col[8][:-1].split(',')
            for pos in poss:
                pos=int(pos)
                bin=int(div(pos,chrbinsize))
                if bin not in gd[cr]['ETS']:gd[cr]['ETS'][bin]={}
                if pos not in gd[cr]['ETS'][bin]:gd[cr]['ETS'][bin][pos]=[]
                gd[cr]['ETS'][bin][pos].append([stra,gname,pos])
        if 'exon' in selsdic:
            if 'exon' not in gd[cr]:gd[cr]['exon']={}
            #if stra=='+':
            gname,cr,stra,starts,ends=col[0],cr,stra,col[8][:-1].split(','),col[9][:-1].split(',')
            #else:gname,cr,stra,starts=col[0],cr,stra,col[9][:-1].split(',')
            num=len(ends)
            for i in range(num):
                #for start in starts:
                start,end=int(starts[i]),int(ends[i])
                #print start,end,chrbinsize
                for pos in range(start,end,chrbinsize):
                    bin=int(div(pos,chrbinsize))
                    if bin not in gd[cr]['exon']:gd[cr]['exon'][bin]={}
                    if start not in gd[cr]['exon'][bin]:gd[cr]['exon'][bin][start]=[]
                    gd[cr]['exon'][bin][start].append([stra,gname,end])
        if 'intron' in selsdic:
            if 'intron' not in gd[cr]:gd[cr]['intron']={}
            #if stra=='+':
            gname,cr,stra,starts,ends=col[0],cr,stra,col[8][:-1].split(','),col[9][:-1].split(',')
            #else:gname,cr,stra,starts=col[0],cr,stra,col[9][:-1].split(',')
            num=len(ends)
            for i in range(1,num):
                #for start in starts:
                start,end=int(ends[i-1]),int(starts[i])
                for pos in range(start,end,chrbinsize):
                    bin=int(div(pos,chrbinsize))
                    if bin not in gd[cr]['intron']:gd[cr]['intron'][bin]={}
                    if start not in gd[cr]['intron'][bin]:gd[cr]['intron'][bin][start]=[]
                    gd[cr]['intron'][bin][start].append([stra,gname,end])
        if 'gene' in selsdic:
            if 'gene' not in gd[cr]:gd[cr]['gene']={}
            gname,cr,stra,start,end=col[0],cr,stra,int(col[3]),int(col[4])
            for pos in range(start,end,chrbinsize):
                bin=int(div(pos,chrbinsize))
                if bin not in gd[cr]['gene']:gd[cr]['gene'][bin]={}
                if start not in gd[cr]['gene'][bin]:gd[cr]['gene'][bin][start]=[]
                gd[cr]['gene'][bin][start].append([stra,gname,end])
            
      
    print(n,'genes')
    '''
    n=0
    for cr in gd:
        for sel in gd[cr]:
            for bin in gd[cr][sel]:
                for start in gd[cr][sel][bin]:
                    n+=len(gd[cr][sel][bin][start])
    #print n, 'genic sites'
    '''
    
    for line in positionLines[1:]:
        genes,ok,col='',0,line.split('\t')
        for sel in sels:
            tgenes={}#sel[0]+":"
            cr,poses=col[0],[]
            for i in ([1,2,3]):
                try:poses.append(int(col[i]))
                except:continue
            minpos,maxpos=min(poses),max(poses)
            minbin,maxbin=int(div((minpos-flank),chrbinsize)),int(div((maxpos+flank),chrbinsize))
            bins=list(range(minbin,maxbin+1))
            for bin in bins:
                if cr not in gd:continue
                if bin not in gd[cr][sel[0]]:continue
                for start in gd[cr][sel[0]][bin]:
                    for i in range(len(gd[cr][sel[0]][bin][start])):
                        end=gd[cr][sel[0]][bin][start][i][2]
                        if gd[cr][sel[0]][bin][start][i][0]=='+' and (maxpos>start+sel[1] and minpos<end+sel[2]):
                            tgenes[gd[cr][sel[0]][bin][start][i][1]+'/'+str(start)+'/'+str(end)]=1
                        elif gd[cr][sel[0]][bin][start][i][0]=='-' and (maxpos>start-sel[2] and minpos<end-sel[1]):
                            tgenes[gd[cr][sel[0]][bin][start][i][1]+'/'+str(start)+'/'+str(end)]=1
            if len(tgenes)!=0:
                ok+=1
                genes+=sel[0]+':'+','.join(list(tgenes.keys()))+'|'
        if andor=='and' and ok==len(sels):
            #print andor,ok
            retr.append(line[:-1]+'\t'+genes[:-1]+'\n')
        elif andor=='or' and ok>0:
            #print andor,ok
            retr.append(line[:-1]+'\t'+genes[:-1]+'\n')
        
    return retr

def positionSelectorByValue(positionLines=[],selection=None):
    '''
    PositionLines:  The first line in positionLines must be the title line, each line should have a '\n' at the end, positionLines should be in the default output format of DANPOS, see DANPOS documentation for example
    selection:  promoter:-350:50&control_smt_val:0:1000&0-log10(point_diff_pval):0:1e-10
    '''
    if selection==None:return positionLines
    if len(positionLines)<2:return positionLines
    retr,tcol,tld=[positionLines[0]],positionLines[0].split('\t'),{}
    if tcol[-1][-1]=='\n':tcol[-1]=tcol[-1][:-1]
    for i in range(len(tcol)):tld[tcol[i]]=i
    sels=selection.split(',')
    #print sels
    if sels[-1]=='and':andor='and'
    elif sels[-1]=='or':andor='or'
    elif len(sels)>1:
        print(sels)
        print("Error: the selection must be defined with 'and' or 'or' at the end")
        return []
    else:
        sels.append('and')
        andor='and'
    sels=sels[:-1]
    #print sels
    for i in range(len(sels)):
        sels[i]=sels[i].split(':')
        if not sels[i][0] in tcol[1:]:
            print("Error:", sels[i][0], "is not a column name in the position file")
            #print tcol
            return []
        if sels[i][1]!='':
            sels[i][1]=float(sels[i][1])
        if sels[i][2]!='':
            sels[i][2]=float(sels[i][2])
    for line in positionLines[1:]:
        ok,col=0,line.split('\t')
        for sel in sels:
            v=float(col[tld[sel[0]]])
            ok1,ok2=False,False
            if sel[1]!='':
                if v>=sel[1]:ok1=True
            else:ok1=True
            if sel[2]!='':
                if v<=sel[2]:ok2=True
            else:ok2=True
            if ok1 and ok2: ok+=1
            #print sel[0],tld[sel[0]],v,sel[1],sel[2],ok1,ok2,ok
        #print ok,len(sels)
        if andor=='and' and ok==len(sels):retr.append(line) #all selection condition must be ok
        elif andor=='or' and ok>0:retr.append(line) #at least one selection condition must be ok
        #print andor, ok
    return retr
 
def retrieve_positions_by_value(in_file=None,out_file=None,cr_col_name='chr',pos_col_name='diff_smt_loca',val_col_name='point_diff_pval',direction_by=['treat_point_val','control_point_val'],top_value=1e-7,bottom_value=0.0,log10trans=False):
    lines=open(in_file).readlines()
    if len(lines)<2:return {}
    col=lines[0].split()
    ids={}
    for i in range(0,len(col)):ids[col[i]]=i
    crid,posid,vid=ids[cr_col_name],ids[pos_col_name],ids[val_col_name]
    if len(direction_by)==2:dirb=[ids[direction_by[0]],ids[direction_by[1]]]###################
    out={}
    retrieve=0
    if out_file:
        fo=open(out_file,'w')
        fo.write(lines[0])
    for line in lines[1:]:
        col=line.split()
        if col[posid]=='-':continue
        try:
            if log10trans:chr,pos,v=col[crid],int(col[posid]),log10(float(col[vid])+1)
            else:chr,pos,v=col[crid],int(col[posid]),float(col[vid])
        except:
            print(line)
            continue
        if top_value!=None and  v>=top_value: continue
        if bottom_value!=None and  v<=bottom_value:continue
        #################
        if len(direction_by)==2:
            if v==0:v=1e-323
            if float(col[dirb[0]])<float(col[dirb[1]]):
                v,col[vid]=0-v,str(0-v)
                line='\t'.join(col)+'\n'
        ###################
        if out_file:fo.write(line)
        if chr not in out:out[chr]={}
        out[chr][pos]=v
        retrieve+=1
    print('\nretrieved',retrieve,'summits out of',len(lines)-1, 'by',val_col_name,bottom_value,'to',top_value)
    return out

def retrieve_positions_by_rank(in_file=None,out_file=None,cr_col_name='chr',pos_col_name='diff_smt_loca',val_col_name='point_diff_pval',toprank=None,bottomrank=None,decreasing=False,direction_by=['treat_point_val','control_point_val']):
    lines=open(in_file).readlines()
    if len(lines)<2:return {}
    col=lines[0].split()
    ids={}
    for i in range(0,len(col)):ids[col[i]]=i
    crid,posid,vid=ids[cr_col_name],ids[pos_col_name],ids[val_col_name]
    if len(direction_by)==2:dirb=[ids[direction_by[0]],ids[direction_by[1]]]###################
    out={}
    if out_file:
        fo=open(out_file,'w')
        fo.write(lines[0])
    tosort={}
    linesdic={}
    for line in lines[1:]:
        col=line.split()
        try:
            tosort[col[crid]+','+col[posid]]=float(col[vid])
            linesdic[col[crid]+','+col[posid]]=line
        except:
            print(line)
            continue
    from operator import itemgetter
    aftersort=sorted(list(tosort.items()),key=itemgetter(1),reverse=decreasing)
    retrieve=0
    if toprank==None:toprank=0#len(aftersort)
    if bottomrank==None:bottomrank=len(aftersort)
    for i in range(toprank,bottomrank):
        chr,pos=aftersort[i][0].split(',')
        v=aftersort[i][1]
        #################
        if len(direction_by)==2:
            if v==0:v=1e-323
            col=linesdic[aftersort[i][0]].split()
            if float(col[dirb[0]])<float(col[dirb[1]]):
                v,col[vid]=0-v,str(0-v)
                linesdic[aftersort[i][0]]='\t'.join(col)+'\n'
        ###################
        if out_file:fo.write(linesdic[aftersort[i][0]])
        if chr not in out:out[chr]={}
        out[chr][int(pos)]=v
        retrieve+=1
    print('\nretrieved',retrieve,'summits out of',len(lines)-1, 'by',col_name)
    return out

def positions2Points(positions={},out_file=None,up=350,down=50,chrColID=1,nameColID=0,posColIDpos=3,posColIDneg=4,straColID=2,sep='\t',second_sep=None,step=1,\
                 neglog=True, rankby='max',lines=[]):
    #positions[chr][pos]=value
    #suggest to download gene table from UCSC genome browser in format:
    #"name    chrom   strand(+/-)  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        proteinID"
    summits=positions
    tsmit={}
    if neglog:
        from math import log10
        for cr in summits:
            tsmit[cr]={}
            for pos in summits[cr]:
                if summits[cr][pos]>0:tsmit[cr][pos]=0-log10(summits[cr][pos])
                elif summits[cr][pos]<0:tsmit[cr][pos]=log10(0-summits[cr][pos])
                else:
                    print('fail to do log transfer at',cr,pos,'due to value',summits[cr][pos])
                    
    summits=tsmit
    
    smts,msmts={},{}
    for chr in summits:
        msmts[chr],smts[chr]={},{}
        for pos in summits[chr]:
            p=pos-(pos%step)#note that the orginal summits located within step will be merged
            if p not in smts[chr]:smts[chr][p]=summits[chr][pos]
            elif rankby=='max':smts[chr][p]=max(smts[chr][p],summits[chr][pos])
            else:smts[chr][p]=min(smts[chr][p],summits[chr][pos])

    if out_file:fo=open(out_file,'w')
    if out_file:fo.write('name\tchr\tpos\tmax_value\tmin_value\tvalues\n')
    out={}
    for line in lines:
        if line[:-1]=='\n':line=line[:-1]
        col=line.split(sep)
        if col[chrColID] not in smts:continue
        tup,tdn,pos,cr=up-up%step,down-down%step,col[posColIDpos],col[chrColID]
        if col[straColID]!='+':tup,tdn,pos,cr=down-down%step,up-up%step,col[posColIDneg],col[chrColID]
        if second_sep==None:poses=[pos]
        else:poses=pos.split(second_sep)
        for pos in poses:
            pos=int(pos)-int(pos)%step
            ps,vs=[],[]
            for p in range(pos-tup,pos+tdn+step,step):
                if p in smts[cr]:
                    msmts[cr][p]=1
                    ps.append(str(p))
                    vs.append(smts[cr][p])
            if len(ps)>0:
                maxv,minv=max(vs),min(vs)
                tvs=[]
                for i in range(0,len(ps)):tvs.append(str(vs[i]))
                line='\t'.join([col[nameColID],cr,str(pos),str(maxv),str(minv),','.join(tvs)])
                if rankby=='min':out[line]=minv
                else:out[line]=maxv
    mcount,tcount,ocount,=0,0,0
    for chr in msmts:mcount+=len(msmts[chr])
    for chr in smts:tcount+=len(smts[chr])
    for chr in summits:ocount+=len(summits[chr])
    print(ocount,'summits merged to',tcount,'by step size',step)
    print('\n',mcount,'of',tcount,'summits mapped to',len(out),'of',len(lines)-1,'genes\n')
    from operator import itemgetter
    if rankby=='min':aftersort=sorted(list(out.items()),key=itemgetter(1))
    else:aftersort=sorted(list(out.items()),key=itemgetter(1),reverse=True)
    rout=[]
    for item in aftersort:
        if out_file:fo.write(item[0]+'\n')
        rout.append(item[0])
    return rout

def plot(dic={'name':{}},outname='',main='',region_size=0,nrow=2,ncol=2,xmin=None,xmax=None,ymin=None,ymax=None,xlab='Relative distance to TSS',ylab='Average occupancy',colors=['black','gray','red','blue','orange','purple','skyblue','cyan','green','blue4','darkgoldenrod'],names=None):
    if main=='':main=outname
    rcode=''
    if names==None:
        names=list(dic.keys())
        names.sort()
    xmincal,xmaxcal,ymincal,ymaxcal=min(dic[names[0]].keys()),max(dic[names[0]].keys()),min(dic[names[0]].values()),max(dic[names[0]].values())
    if len(colors)<len(names):
        print('Wrong:please specify ',len(names),'colors for the curves')
        return ''
    for i in range(len(names)):
        name=names[i]
        poses,vals=list(dic[name].keys()),list(dic[name].values())
        txmin,txmax,tymin,tymax=min(poses),max(poses),min(vals),max(vals)
        if xmincal>txmin:xmincal=txmin
        if xmaxcal<txmax:xmaxcal=txmax
        if ymincal>tymin:ymincal=tymin
        if ymaxcal<tymax:ymaxcal=tymax
        poses,vals=[],[]
        tposes=list(dic[name].keys())
        tposes.sort()
        for pos in tposes:
            poses.append(str(pos))
            vals.append(str(dic[name][pos]))
        rcode+='lines(c('+','.join(poses)+'),c('+','.join(vals)+'),col="'+colors[i]+'")\n'
    rcode+='legend("topright",legend=c("'+'","'.join(names)+'"),col=c("'+'","'.join(colors[0:len(names)])+'"),lty=1)\n'
    if xmin!=None:xmincal=xmin
    if xmax!=None:xmaxcal=xmax
    if ymin!=None:ymincal=ymin
    if ymax!=None:ymaxcal=ymax
    else:ymaxcal+=(ymaxcal-ymincal)*(len(names)*0.12+0.1)

    if region_size==0:rcode='plot(0,0,type="n",main="'+main+'",xlim=c('+str(xmincal)+','+str(xmaxcal)+'),'+'ylim=c('+str(ymincal)+','+str(ymaxcal)+'),xlab="'+str(xlab)+'",ylab="'+str(ylab)+'",)\n'+rcode
    else:
        rcode='plot(0,0,type="n",main="'+main+'",xaxt="n",xlim=c('+str(xmincal)+','+str(xmaxcal)+'),'+'ylim=c('+str(ymincal)+','+str(ymaxcal)+'),xlab="'+str(xlab)+'",ylab="'+str(ylab)+'",)\n'+rcode
        poses=list(dic[names[0]].keys())
        poses.sort()
        at,lb=['0',str(region_size)],['\"Start\"','\"End\"']
        #print at,lb
        rcode+='axis(side=1,at=c('+','.join(at)+'),labels=c('+','.join(lb)+'))\n'
        lth=poses[-1]-poses[0]+poses[1]-poses[0]
        #print poses[-1],poses[0]
        lth=int(div(lth,6))
        #tlth=str(lth)
        #if len(tlth)>3:
        at,lb=[],[]
        if poses[0]<0:
            for pos in range(lth,0-poses[0]+lth,lth):
                at.append('-'+str(pos))
                lb.append('\"-'+str(pos)+'\"')
        rcode+='axis(side=1,at=c('+','.join(at)+'),labels=c('+','.join(lb)+'))\n'
        at,lb=[],[]
        for pos in range(lth,region_size-lth+1,lth):
                at.append(str(pos))
                lb.append('\"'+str(pos)+'\"')
        rcode+='axis(side=1,at=c('+','.join(at)+'),labels=c('+','.join(lb)+'))\n'
        at,lb=[],[]
        for pos in range(region_size+lth,poses[-1]+lth,lth):
                at.append(str(pos))
                lb.append('"+'+str(pos-region_size)+'\"')
        #print at,lb
        rcode+='axis(side=1,at=c('+','.join(at)+'),labels=c('+','.join(lb)+'))\n'

    if outname!='':
        rcode='par(mfrow=c('+str(nrow)+','+str(ncol)+'))\n'+rcode
        rcode='pdf("'+outname+'.pdf")\n'+rcode
        rcode+='dev.off()\n'
        fo=open(outname+'.R','w')
        fo.write(rcode)
        fo.close()
        r(rcode)
    return rcode
def vioplot(dic={'name':[]},outname='',main='',nrow=2,ncol=2,ymin=None,ymax=None,xlab='Relative distance to TSS',ylab='Average occupancy',colors=['black','gray','red','blue','orange','purple','skyblue','cyan','green','blue4','darkgoldenrod'],names=None):
    if main=='':main=outname
    if names==None:
        names=list(dic.keys())
        names.sort()
    #xmincal,xmaxcal,ymincal,ymaxcal=min(dic[names[0]].keys()),max(dic[names[0]].keys()),min(dic[names[0]].values()),max(dic[names[0]].values())
    if len(colors)<len(names):
        print('please specify ',len(names),'colors for the curves')
        return ''
    rcode="library('vioplot')\nvioplot("
    for name in names:
        temp=[]
        for value in dic[name]:temp.append(str(value))
        rcode+="c("+','.join(temp)+"),"
    rcode+="ylim=c("+str(ymin)+","+str(ymax)+"),names=c("
    for name in names:rcode+="'"+name+"',"
    rcode=rcode[:-1]
    rcode+="))\n"
    rcode+="mtext('"+main+"')\n"
    if outname!='':
        rcode='par(mfrow=c('+str(nrow)+','+str(ncol)+'))\n'+rcode
        rcode='pdf("'+outname+'.pdf")\n'+rcode
        rcode+='dev.off()\n'
        fo=open(outname+'.R','w')
        fo.write(rcode)
        fo.close()
        r(rcode)
    return rcode
    
def occPSD0(wg,outname=None):
    psd=r('''function(q){return(spec.pgram(q,plot = FALSE)$spec)}''')
    lth=div(100000,wg.step)
    d=wg.data
    spe=[0.0]*(div(lth,2))
    wn=0
    #print 'calculating spectrum'
    for cr in d:
        print(cr)
        sz=d[cr].size
        for i in range(0,sz-lth,div(lth,2)):
            wn+=1
            if wn%100==0:print(wn,'window calculated ...')
            v=psd(FloatVector(d[cr][i:(i+lth)]))
            for j in range(div(lth,2)):spe[j]+=v[j]
        print(wn,'window calculated.')
    if outname!=None:fo=open(outname+'.xls','w')
    dic={}
    for j in range(int(div(lth*wg.step,250)),int(div(lth*wg.step,100)+1)):
        dic[div(lth*wg.step*1.0,j)]=div(spe[j],wn)
        if outname!=None:fo.write(str(div(lth*wg.step*1.0,j))+'\t'+str(div(spe[j],wn))+'\n')
    return dic
def occPSD(wg,outname=None):
    cor=r('''function(q1,q2){return(cor(q1,q2))}''')
    #lth=100000/wg.step
    d=wg.data
    #spe=[0.0]*(lth/2)
    #wn=0
    #print 'calculating spectrum'
    tsz,mi,ma,dic=0,100,250,{}
    for i in range(mi,ma):dic[i]=0
    
    for cr in d:
        print(cr)
        sz=d[cr].size-ma-1
        tsz+=sz
        for i in range(mi,ma):
            v=cor(FloatVector(d[cr][0:sz]),FloatVector(d[cr][i:(i+sz)]))
            v=float(str(v).split()[1])
            dic[i]+=v*sz
    if outname!=None:fo=open(outname+'.xls','w')
    for i in range(mi,ma):
        dic[i]=div(dic[i],tsz)
        if outname!=None:fo.write(str(i)+'\t'+str(dic[i])+'\n')
    return dic
def positionDistanceDistribution(dic,outname=None,min=100,max=250,step=1):
    max=max+1
    '''
    dic={}
    ct=0
    for line in lines:
        ct+=1
        col=line.split()
        if not dic.has_key(col[crColID]):dic[col[crColID]]={}
        dic[col[crColID]][int(col[posColID])]=1#float(col[valColID])
    '''
    ct,dis=0,{}
    for cr in dic:
        poses=list(dic[cr].keys())
        ct+=len(poses)
        poses.sort()
        lth=len(poses)
        for i in range(1,lth):
            d=int(div((poses[i]-poses[i-1]),step))
            if d not in dis:dis[d]=1
            else:dis[d]+=1
    if outname!=None:fo=open(outname+'.xls','w')
    if min==None:min=min(dis.keys())
    if max==None:max=max(dis.keys())
    odic={}
    for d in range(int(div(min,step)),int(div(max,step))):
        if d not in dis:dis[d]=0
        dis[d]=div(dis[d]*100.0,ct)  #change to percentage
        odic[d*step]=dis[d]
        if outname!=None:fo.write( str(d*step)+'\t'+str(dis[d])+'\n')
    return odic
def positionDistance(dic,outname=None,min=0,max=350):
    max=max+1
    '''
    dic={}
    ct=0
    for line in lines:
        ct+=1
        col=line.split()
        if not dic.has_key(col[crColID]):dic[col[crColID]]={}
        dic[col[crColID]][int(col[posColID])]=1#float(col[valColID])
    '''
    ct,dis=0,{}
    for cr in dic:
        dis[cr]={}
        poses=list(dic[cr].keys())
        poses.sort()
        lth=len(poses)
        for i in range(1,lth):
            #print poses[i],poses[i]-poses[i-1]
            #if poses[i]-poses[i-1]>min and poses[i]-poses[i-1]<max:dis[cr][poses[i-1]]=poses[i]-poses[i-1]
            dis[cr][poses[i-1]]=poses[i]-poses[i-1]
    '''
            if not dis.has_key(d):dis[d]=1
            else:dis[d]+=1
    if outname!=None:fo=open(outname+'.xls','w')
    if min==None:min=min(dis.keys())
    if max==None:max=max(dis.keys())
    odic={}
    for d in range(int(min/step),int(max/step)):
        if not dis.has_key(d):dis[d]=0
        dis[d]=dis[d]*100.0/ct  #change to percentage
        odic[d*step]=dis[d]
        if outname!=None:fo.write( str(d*step)+'\t'+str(dis[d])+'\n')
    '''
    return dis

def positionValDistribution(dic,outname=None,min=0,max=1500,step=3):
    #for occupancy:valColID=4,min=0,max=1500,step=3
    #for fuzziness:valColID=6,min=15,max=25,step=0.01
    
    vdic={}
    ct=0
    for cr in dic:
        ct+=len(dic[cr])
        for v in list(dic[cr].values()):
            v=int(div(v,step)+0.5)
            if v not in vdic:vdic[v]=1
            else:vdic[v]+=1
    if ct<1:return {}
    if outname!=None:fo=open(outname+'.xls','w')
    odic={}
    #print int(min/step),int(max/step)
    for d in range(int(div(min,step)),int(div(max,step))):
        if d not in vdic:vdic[d]=0
        vdic[d]=div(vdic[d]*100.0,ct)  #change to percentage
        if outname!=None:fo.write( str(d*step)+'\t'+str(vdic[d])+'\n')
        odic[d*step]=vdic[d]
    rodic={}
    for k in odic:
        if odic[k]>0:rodic[k]=odic[k]
    return rodic

def positionAroundPoint(smts,outname=None,flankup=2500,flankdn=2500,step=10,chrColID=1,nameColID=0,posColIDpos=3,posColIDneg=4,straColID=2,sep='\t',second_sep=None,lines=None):
    #calculate the percentage of summits distributed in falnking region of TTS
    flankup-=flankup%step
    flankdn-=flankdn%step
    tsmtdic=smts#.data
    smtdic={}
    for chr in list(tsmtdic.keys()):
        smtdic[chr]={}
        for pos in tsmtdic[chr]:
            tpos=pos-pos%step
            smtdic[chr][tpos]=tsmtdic[chr][pos]
    dis={}
    num=0
    for line in lines:
        col=line.split()
        chr,stra=col[chrColID],col[straColID]
        tss=col[posColIDpos]
        if stra!='+':tss=col[posColIDneg]
        if second_sep==None:tsses=[tss]
        else:tsses=tss.split(second_sep)
        for tss in tsses:
            if tss=='':continue
            tss=int(tss)
            tss-=tss%step
            if chr not in smtdic:continue
            num+=1
            if stra=='+':
                for pos in range(tss-flankup,tss+flankdn,step):
                    if pos in smtdic[chr]:
                        if pos-tss in dis:dis[pos-tss]+=1#smtdic[chr][pos]
                        else: dis[pos-tss]=1#smtdic[chr][pos]
            elif stra=='-':
                for pos in range(tss-flankdn+step,tss+flankup+step,step):
                    if pos in smtdic[chr]:
                        if tss-pos in dis:dis[tss-pos]+=1#smtdic[chr][pos]
                        else:dis[tss-pos]=1#smtdic[chr][pos]
    num*=step
    if outname!=None:outf=open(outname+'.xls',"w")
    for k in range(0-flankup,flankdn,step):
        if k not in dis:dis[k]=0
        dis[k]=div(dis[k]*1.0,num)
        if outname!=None:outf.write(str(k)+"\t"+str(dis[k])+"\n")
    return dis

def runall(inpath=None,outpath=None,flankup=3000,flankdn=3000,minfuz=15,maxfuz=25,minocc=0,maxocc=1500,step=10,lines=None):
    
    '''
    if outpath==None:outpath=inpath
    for file in glob(os.path.join(inpath,'*positions.differential.xls')):
        pks=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='point_diff_pval',toprank=None,bottomrank=None,decreasing=True,direction_by=['treat_point_val','control_point_val'])
        positions2Points(positions=pks,out_file=os.path.join(outpath,os.path.split(file)[-1][:-3]+'2tss.xls'),up=350,down=50,chrColID=1,nameColID=0,posColIDpos=3,posColIDneg=4,straColID=2,sep='\t',second_sep=None,step=1,\
                     neglog=True, rankby='max',lines=lines)
    '''
    
    '''
    if outpath==inpath:outpath=None
    pth=os.path.join(inpath,'diff')
    if outpath==None:outpath=pth
    smts={}
    for file in glob(os.path.join(pth,'*positions.xls')):
        fn=os.path.split(file)[-1]
        fn=re.sub('positions.xls','',fn)
        fn=re.sub('\.\w+\_diff','',fn)
        smts[fn]=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='smt_value',toprank=0,bottomrank=2000,decreasing=True,direction_by=[])
    dic=batchPositionAroundPoints(smts,outname=os.path.join(outpath,'diff_position_ard_ESS'),flankup=flankup,flankdn=flankdn,step=step,chrColID=1,nameColID=0,posColIDpos=8,posColIDneg=9,straColID=2,sep='\t',second_sep=',',lines=lines)
    plot(dic=dic,outname=os.path.join(outpath,'diff_position_ard_ESS'),xlab='distance to ESS',ylab='position count')
    
    smts={}
    for file in glob(os.path.join(pth,'*positions.xls')):
        fn=os.path.split(file)[-1]
        fn=re.sub('positions.xls','',fn)
        fn=re.sub('\.\w+\_diff','',fn)
        smts[fn]=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='smt_value',toprank=0,bottomrank=2000,decreasing=True,direction_by=[])
    dic=batchPositionAroundPoints(smts,outname=os.path.join(outpath,'diff_position_ard_ETS'),flankup=flankup,flankdn=flankdn,step=step,chrColID=1,nameColID=0,posColIDpos=9,posColIDneg=8,straColID=2,sep='\t',second_sep=',',lines=lines)
    plot(dic=dic,outname=os.path.join(outpath,'diff_position_ard_ETS'),xlab='distance to ETS',ylab='position count')
    
    smts={}
    for file in glob(os.path.join(pth,'*positions.xls')):
        fn=os.path.split(file)[-1]
        fn=re.sub('positions.xls','',fn)
        fn=re.sub('\.\w+\_diff','',fn)
        smts[fn]=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='smt_value',toprank=0,bottomrank=2000,decreasing=True,direction_by=[])
    dic=batchPositionAroundPoints(smts,outname=os.path.join(outpath,'diff_position_ard_TSS'),flankup=flankup,flankdn=flankdn,step=step,chrColID=1,nameColID=0,posColIDpos=3,posColIDneg=4,straColID=2,sep='\t',second_sep=',',lines=lines)
    plot(dic=dic,outname=os.path.join(outpath,'diff_position_ard_TSS'),xlab='distance to TSS',ylab='position count')
    
    smts={}
    for file in glob(os.path.join(pth,'*positions.xls')):
        fn=os.path.split(file)[-1]
        fn=re.sub('positions.xls','',fn)
        fn=re.sub('\.\w+\_diff','',fn)
        smts[file]=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='smt_value',toprank=0,bottomrank=2000,decreasing=True,direction_by=[])
    dic=batchPositionAroundPoints(smts,outname=os.path.join(outpath,'diff_position_ard_TTS'),flankup=flankup,flankdn=flankdn,step=step,chrColID=1,nameColID=0,posColIDpos=4,posColIDneg=3,straColID=2,sep='\t',second_sep=',',lines=lines)
    plot(dic=dic,outname=os.path.join(outpath,'diff_position_ard_TTS'),xlab='distance to TTS',ylab='position count')  
    
    smts={}
    for file in glob(os.path.join(pth,'*positions.xls')):
        fn=os.path.split(file)[-1]
        fn=re.sub('positions.xls','',fn)
        fn=re.sub('\.\w+\_diff','',fn)
        smts[fn]=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='smt_value',toprank=0,bottomrank=2000,decreasing=True,direction_by=[])
    dic=batchPositionValDistribution(data=smts,outname=os.path.join(outpath,'diff_position_Value_distribution'),min=minocc,max=maxocc,step=1)
    plot(dic=dic,outname=os.path.join(outpath,'diff_position_Value_distribution'),xlab='-log P value',ylab='Percentage')
    smts={}
    for file in glob(os.path.join(pth,'*positions.xls')):
        fn=os.path.split(file)[-1]
        fn=re.sub('positions.xls','',fn)
        fn=re.sub('\.\w+\_diff','',fn)
        smts[fn]=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='smt_value',toprank=0,bottomrank=2000,decreasing=True,direction_by=[])
    dic=batchPositionDistanceDistribution(data=smts,outname=os.path.join(outpath,'diff_position_distance_distribution'),min=100,max=250,step=1)
    plot(dic=dic,outname=os.path.join(outpath,'diff_position_distance_distribution'),xlab='distance',ylab='Percentage')
    
    wgs=Wigs(pth)
    for k in wgs.keys():
        fn=os.path.split(k)[-1]
        fn=re.sub('\.\w+\_diff','',fn)
        print fn
        twg=deepcopy(wgs.get(k))
        twg.rvNeg()
        wgs.set(fn[:-3]+'gain',twg)
        wgs.get(k).foldChange(-1)
        wgs.get(k).rvNeg()
        wgs.set(fn[:-3]+'loss',wgs.pop(k))
    for wg in wgs.keys():print wg,wgs.get(wg).sum()
    
    tts=batchOccAroundPoints(wgs=wgs,outname=os.path.join(outpath,'diff_value_ard_TTS'),flankup=flankup,flankdn=flankdn,chrColID=1,nameColID=0,posColIDpos=4,posColIDneg=3,straColID=2,sep='\t',second_sep=None,lines=lines)
    plot(dic=tts,nrow=2,ncol=2,outname=os.path.join(outpath,'diff_value_ard_TTS'),xlab='Relative distance to TTS',ylab='Average value')
    tss=batchOccAroundPoints(wgs=wgs,outname=os.path.join(outpath,'diff_value_ard_TSS'),flankup=flankup,flankdn=flankdn,chrColID=1,nameColID=0,posColIDpos=3,posColIDneg=4,straColID=2,sep='\t',second_sep=None,lines=lines)
    plot(dic=tss,nrow=2,ncol=2,outname=os.path.join(outpath,'diff_value_ard_TSS'),xlab='Relative distance to TSS',ylab='Average value')
    ess=batchOccAroundPoints(wgs=wgs,outname=os.path.join(outpath,'diff_value_ard_ESS'),flankup=flankup,flankdn=flankdn,chrColID=1,nameColID=0,posColIDpos=8,posColIDneg=9,straColID=2,sep='\t',second_sep=',',lines=lines)
    plot(dic=ess,nrow=2,ncol=2,outname=os.path.join(outpath,'diff_value_ard_ESS'),xlab='Relative distance to ESS',ylab='Average value')
    ets=batchOccAroundPoints(wgs=wgs,outname=os.path.join(outpath,'diff_value_ard_ETS'),flankup=flankup,flankdn=flankdn,chrColID=1,nameColID=0,posColIDpos=9,posColIDneg=8,straColID=2,sep='\t',second_sep=',',lines=lines)
    plot(dic=ets,nrow=2,ncol=2,outname=os.path.join(outpath,'diff_value_ard_ETS'),xlab='Relative distance to ETS',ylab='Average value')
    '''
    
    #if outpath==pth:outpath=None
    pth=os.path.join(inpath,'pooled')
    if outpath==None:outpath=pth
    '''
    smts={}
    for file in glob(os.path.join(pth,'*positions.xls')):
        fn=os.path.split(file)[-1]
        fn=re.sub('positions.xls','',fn)
        smts[file]=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='fuzziness_score',toprank=None,bottomrank=None,decreasing=False,direction_by=[])
    dic=batchPositionAroundPoints(smts,outname=os.path.join(outpath,'position_ard_ESS'),flankup=flankup,flankdn=flankdn,step=step,chrColID=1,nameColID=0,posColIDpos=8,posColIDneg=9,straColID=2,sep='\t',second_sep=',',lines=lines)
    plot(dic=dic,outname=os.path.join(outpath,'position_ard_ESS'),xlab='distance to ESS',ylab='position count')
    smts={}
    for file in glob(os.path.join(pth,'*positions.xls')):
        fn=os.path.split(file)[-1]
        fn=re.sub('positions.xls','',fn)
        smts[file]=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='fuzziness_score',toprank=None,bottomrank=None,decreasing=False,direction_by=[])
    dic=batchPositionAroundPoints(smts,outname=os.path.join(outpath,'position_ard_ETS'),flankup=flankup,flankdn=flankdn,step=step,chrColID=1,nameColID=0,posColIDpos=9,posColIDneg=8,straColID=2,sep='\t',second_sep=',',lines=lines)
    plot(dic=dic,outname=os.path.join(outpath,'position_ard_ETS'),xlab='distance to ETS',ylab='position count')
    smts={}
    for file in glob(os.path.join(pth,'*positions.xls')):
        fn=os.path.split(file)[-1]
        fn=re.sub('positions.xls','',fn)
        smts[file]=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='fuzziness_score',toprank=None,bottomrank=None,decreasing=False,direction_by=[])
    dic=batchPositionAroundPoints(smts,outname=os.path.join(outpath,'position_ard_TSS'),flankup=flankup,flankdn=flankdn,step=step,chrColID=1,nameColID=0,posColIDpos=3,posColIDneg=4,straColID=2,sep='\t',second_sep=',',lines=lines)
    plot(dic=dic,outname=os.path.join(outpath,'position_ard_TSS'),xlab='distance to TSS',ylab='position count')
    smts={}
    for file in glob(os.path.join(pth,'*positions.xls')):
        fn=os.path.split(file)[-1]
        fn=re.sub('positions.xls','',fn)
        smts[file]=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='fuzziness_score',toprank=None,bottomrank=None,decreasing=False,direction_by=[])
    dic=batchPositionAroundPoints(smts,outname=os.path.join(outpath,'position_ard_TTS'),flankup=flankup,flankdn=flankdn,step=step,chrColID=1,nameColID=0,posColIDpos=4,posColIDneg=3,straColID=2,sep='\t',second_sep=',',lines=lines)
    plot(dic=dic,outname=os.path.join(outpath,'position_ard_TTS'),xlab='distance to TTS',ylab='position count')  
    '''
    
    
    smts={}
    for file in glob(os.path.join(pth,'*positions.xls')):
        fn=os.path.split(file)[-1]
        fn=re.sub('positions.xls','',fn)
        smts[file]=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='fuzziness_score',toprank=None,bottomrank=None,decreasing=False,direction_by=[])
    dic=batchPositionValDistribution(data=smts,outname=os.path.join(outpath,'position_fuzziness_distribution'),min=minfuz,max=maxfuz,step=0.01)
    plot(dic=dic,outname=os.path.join(outpath,'position_fuzziness_distribution'),nrow=1,ncol=2,xlab='Fuzziness',ylab='Percentage')
    
    smts={}
    for file in glob(os.path.join(pth,'*positions.xls')):
        fn=os.path.split(file)[-1]
        fn=re.sub('positions.xls','',fn)
        smts[file]=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='smt_value',toprank=None,bottomrank=None,decreasing=False,direction_by=[])
    dic=batchPositionValDistribution(data=smts,outname=os.path.join(outpath,'position_occupancy_distribution'),min=minocc,max=maxocc,step=3)
    plot(dic=dic,outname=os.path.join(outpath,'position_occupancy_distribution'),xlab='occupancy',ylab='Percentage')
    
    smts={}
    for file in glob(os.path.join(pth,'*positions.xls')):
        fn=os.path.split(file)[-1]
        fn=re.sub('positions.xls','',fn)
        smts[file]=retrieve_positions_by_rank(in_file=file,out_file=None,col_name='smt_value',toprank=None,bottomrank=None,decreasing=False,direction_by=[])
    dic=batchPositionDistanceDistribution(data=smts,outname=os.path.join(outpath,'position_distance_distribution'),min=50,max=400,step=step)
    plot(dic=dic,outname=os.path.join(outpath,'position_distance_distribution'),xlab='distance',ylab='Percentage')
    
    '''
    wgs=Wigs(pth,step=step)
    for wg in wgs.keys():print wg,wgs.get(wg).sum()
    tts=batchOccAroundPoints(wgs=wgs,outname=os.path.join(outpath,'occ_ard_random'),flankup=flankup,flankdn=flankdn,chrColID=2,nameColID=1,posColIDpos=5,posColIDneg=4,straColID=3,sep='\t',second_sep=None,lines=lines)
    plot(dic=tts,nrow=2,ncol=2,outname=os.path.join(outpath,'occ_ard_random'),xlab='Relative distance to random sites',ylab='Average occupancy')
        
    for wg in wgs.keys():print wg,wgs.get(wg).sum()
    tts=batchOccAroundPoints(wgs=wgs,outname=os.path.join(outpath,'occ_ard_TTS'),flankup=flankup,flankdn=flankdn,chrColID=2,nameColID=1,posColIDpos=5,posColIDneg=4,straColID=3,sep='\t',second_sep=None,lines=lines)
    plot(dic=tts,nrow=2,ncol=2,outname=os.path.join(outpath,'occ_ard_TTS'),xlab='Relative distance to TTS',ylab='Average occupancy')
    
    tss=batchOccAroundPoints(wgs=wgs,outname=os.path.join(outpath,'occ_ard_TSS'),flankup=flankup,flankdn=flankdn,chrColID=2,nameColID=1,posColIDpos=4,posColIDneg=5,straColID=3,sep='\t',second_sep=None,lines=lines)
    plot(dic=tss,nrow=2,ncol=2,outname=os.path.join(outpath,'occ_ard_TSS'),xlab='Relative distance to TSS',ylab='Average occupancy')
    
    ess=batchOccAroundPoints(wgs=wgs,outname=os.path.join(outpath,'occ_ard_ESS'),flankup=flankup,flankdn=flankdn,chrColID=2,nameColID=1,posColIDpos=9,posColIDneg=10,straColID=3,sep='\t',second_sep=',',lines=lines)
    plot(dic=ess,nrow=2,ncol=2,outname=os.path.join(outpath,'occ_ard_ESS'),xlab='Relative distance to ESS',ylab='Average occupancy')
    
    ets=batchOccAroundPoints(wgs=wgs,outname=os.path.join(outpath,'occ_ard_ETS'),flankup=flankup,flankdn=flankdn,chrColID=2,nameColID=1,posColIDpos=10,posColIDneg=9,straColID=3,sep='\t',second_sep=',',lines=lines)
    plot(dic=ets,nrow=2,ncol=2,outname=os.path.join(outpath,'occ_ard_ETS'),xlab='Relative distance to ETS',ylab='Average occupancy')
    
    psd=batchOccPSD(wgs,outname=os.path.join(outpath,'psd'))
    plot(dic=psd,outname=os.path.join(outpath,'psd'),xlab='Periodicity unit length',ylab='Strength')
    '''   
    print('job done\n')
    
def randomTSS(genomefile=None, genefile=None):
    gd={}
    from random import randint
    for line in open(genomefile):
        col=line.split()
        gd[col[0]]=int(col[1])-1
    '''
    stras=['+','-']
    n=0
    for cr in gd:
        for i in range(gd[col[0]]/10000):
            n+=1
            start,end,stra=str(randint(0,gd[cr])),str(randint(0,gd[cr])),stras[randint(0,1)]
            print '\t'.join([str(n),str(n),cr,stra,start,end,'end'])
    '''
    fi=open(genefile)
    print(fi.readline()[:-1])
    for line in fi:
        col=line[:-1].split('\t')
        col[4]=str(randint(0,gd[col[2]]))
        col[5]=str(randint(0,gd[col[2]]))
        print('\t'.join(col))
    
def positionDicMinMax(dic,lowPercent=0,highPercent=100):
    outmin,outmax,values,dvalues=0,0,[],{}
    for name in dic:
        for chr in dic[name]:
            values+=list(dic[name][chr].values())
    values.sort()
    lth=len(values)
    outmin,outmax=values[int(div(lth*lowPercent,100.0))],values[int(div(lth*highPercent,100.0))-1]

    for name in dic:
        dvalues[name]=[]
        for chr in dic[name]:
            #values+=dic[name][chr].values()
            for value in list(dic[name][chr].values()):
                if value>outmin and value<outmax:dvalues[name].append(value)
    return [outmin,outmax,dvalues]

def translocationReads(file,bindic={},binSize=1000000,outReadsFile='out.sam',outmode='a',pdis=3000,step=1,mapq=30,clipSize=3,inter=True,intra=True,allWigFile=None,transWigFile=None,readsCount={}):
    #if inter and intra:outReadsFile,transWigFile,allWigFile=os.path.join(name,file[:-3]+'trans.sam'),os.path.join(name,file[:-3]+'trans.wig'),os.path.join(name,file[:-3]+'all.wig')
    #elif inter:outReadsFile,transWigFile,allWigFile=os.path.join(name,file[:-3]+'inter.sam'),os.path.join(name,file[:-3]+'inter.wig'),os.path.join(name,file[:-3]+'all.wig')
    #elif intra:outReadsFile,transWigFile,allWigFile=os.path.join(name,file[:-3]+'intra.sam'),os.path.join(name,file[:-3]+'intra.wig'),os.path.join(name,file[:-3]+'all.wig')
    #else:return []
    transwig=Wig(step=step)
    fo=open(outReadsFile,outmode)
    infile=open(file)#os.popen('samtools view -XS '+file)
    line=infile.readline()
    hlines=[]
    while line[0]=='@':
        col=line.split()
        if col[0]=='@SQ':
            chr,clen=col[1][3:],int(col[2][3:])+step
            if chr not in bindic:bindic[chr]={}
            for bin in range(div(clen,binSize)+1):
                if bin not in bindic[chr]:bindic[chr][bin]={}
            transwig.data[chr]=numpy.array([0.0])
            transwig.data[chr].resize(div(clen,step),refcheck=0)
        hlines.append(line)
        line=infile.readline()
    #rdlth=len(line.split('\t')[9])
    infile=open(file)
    for i in range(len(hlines)):infile.readline()
    
    #if allWigFile!=None:
    allwig=deepcopy(transwig)
    nt,ne,na,no,nl,nu=-1,0,0,0,0,0
    line1,line2='',''
    for line in infile:#.readline()
        nt+=1
        if nt%2==0:
            line1=line
            continue
        else:line2=line
        col1,col2=line1[:-1].split('\t'),line2[:-1].split('\t')
        chr1,chr2,mapq1,mapq2,pos1,pos2=col1[2],col2[2],float(col1[4]),float(col2[4]),int(col1[3]),int(col2[3])
        if chr1=='*' or chr2=='*':
            nu+=2 #un-mappable
            continue
        unique=True
        if mapq1<mapq or mapq2<mapq:nl,unique=nl+2,False #low mapping-quality ->non-unique
        tra,ter=False,False
        if chr1==chr2:
            dis=pos1-pos2
            if dis<0:dis=0-dis
            if dis>pdis:
                tra=True
        else:
            ter=True
        trans=False
        if inter and ter:trans=True
        elif intra and tra:trans=True
        bre=False #will be set to True if clipSize<=0 or any of the reads in the pair has soft clip
        if col1[5]=='*' or col2[5]=='*':continue
        for line in [line1,line2]:
            col=line.split('\t')
            rdlth=len(col[9])
            t1=re.findall('\d+',col[5])
            t2=re.findall('\D+',col[5])
            '''
            if not unique:#len(t1)!=len(t2):
                print line
                continue
            '''
            start=int(col[3])
            tbre=False# will set bre to True if tbre become True later.
            if clipSize<=0:tbre=True
            else:
                if t2[0]=='S':
                    if int(t1[0])>=clipSize:tbre=True
                elif t2[-1]=='S':
                    if int(t1[-1])>=clipSize:tbre=True
            if tbre:bre=True#set bre to True if tbre become True.
            mlth=0
            for i in range(len(t2)):
                t1[i]=int(t1[i])
                if t2[i]=='M' :mlth+=t1[i]
            
            for i in range(len(t2)):
                if t2[i]=='M' :#all matched nucleotide will be counted into wiggle format data.
                    end=start+t1[i]
                    allwig.data[col[2]][div(start,step):div(end,step)]+=div(1.0,mlth)
                    if unique and tbre and trans:transwig.data[col[2]][div(start,step):div(end,step)]+=1.0#/rdlth
                    start=end
                elif t2[i]=='D':start=start+t1[i]
                elif t2[i]=='N':start=start+t1[i]
            
        if bre and trans:
                fo.write(line1+line2)
                if ter and bre:ne+=2 #inter-chromosome translocation
                elif tra and bre:na+=2 #intra-chromosome translocation
        else:no+=2
    
    fo.close()
    if allWigFile!=None:allwig.save(allWigFile)
    if transWigFile!=None:transwig.save(transWigFile)
    
    print('all raw reads:',nt+1)
    print('unmappable:',nu,div(nu*100.0,nt),'%')
    print('low map quality (non-unique):',nl,div(nl*100.0,nt),'%')
    if inter:print('inter-chromosome translocated and unique:',ne,div(ne*100.0,nt),'%')
    if intra:print('intra-chromosome translocated and unique:',na,div(na*100.0,nt),'%')
    print('other unique:',no,div(no*100.0,nt),'%')
    print('All unique:',ne+na+no,div((ne+na+no)*100.0,nt),'%')
    if 'all' not in readsCount:readsCount['all']=nt
    else:readsCount['all']+=nt
    if 'unmappable' not in readsCount:readsCount['unmappable']=nu
    else:readsCount['unmappable']+=nu
    if 'non_unique' not in readsCount:readsCount['non_unique']=nl
    else:readsCount['non_unique']+=nl
    if 'unique_inter' not in readsCount:readsCount['unique_inter']=ne
    else:readsCount['unique_inter']+=ne
    if 'unique_intra' not in readsCount:readsCount['unique_intra']=na
    else:readsCount['unique_intra']+=na
    if 'unique_other' not in readsCount:readsCount['unique_other']=no
    else:readsCount['unique_other']+=no
    if 'unique' not in readsCount:readsCount['unique']=ne+na+no
    else:readsCount['unique']+=ne+na+no
    if 'mappable' not in readsCount:readsCount['mappable']=ne+na+no+nl
    else:readsCount['mappable']+=ne+na+no+nl
    if 'trans' not in readsCount:readsCount['trans']=0
    if inter:readsCount['trans']+=ne#readsCount['unique_inter']
    if intra:readsCount['trans']+=na#readsCount['unique_intra']
    #if allWigFile!=None:return [transwig,allwig,]#trans reads count, all reads coverage, trans reads coverage
    return [transwig,allwig]
    
def translocationLinks(peaks,samFile,linkfile,bindic={},binSize=1000000,wsize=500,wstep=0,fold=0,logP=0):
    #fisherTest=r('''function(x,m,k,t){return(phyper(x - 1, m, t-m, k, lower.tail = FALSE,log.p = TRUE)/log(10))}''')
    #fisherTest=r('''function(ov,n1,n2,n){return(phyper(ov - 1, n1, n-n1, n2, lower.tail = FALSE,log.p = TRUE)/log(10))}''')
    fisherTest=r('''function(ov,n1,n2,n){return(phyper(ov - 1, n1, n-n1, n2, lower.tail = FALSE))}''')
    pks,lks=bindic,{}
    infile=open(samFile)#os.popen('samtools view -XS '+file)
    line=infile.readline()
    hlines=[]
    while line[0]=='@':
        hlines.append(line)
        line=infile.readline()
    
    infile=open(samFile)
    for i in range(len(hlines)):infile.readline()
    
    pklines=['-\t-\t-\t-\t-']#'\t'.join(title.split('\t')[:3]+['len'])]
    id=0
    bn={}
    bn[0]=0
    for cr in peaks:
        starts=list(peaks[cr].keys())
        starts.sort()
        for start in starts:
            end = peaks[cr][start]
            if wstep>0:
                for nstart in range(start-div(wsize,2),end-div(wsize,2),wstep):
                    if nstart<0:nstart=0
                    id+=1
                    lks[id],bn[id]={},0
                    pklines.append('\t'.join([cr,str(start),str(end),str(nstart),str(nstart+wsize)]))
                    sBin,eBin=div(nstart,binSize),div((nstart+wsize),binSize)
                    for bin in range(sBin,eBin+1):
                        if bin not in pks[cr]:pks[cr][bin]={}
                        pks[cr][bin][nstart]=[nstart+wsize,id]
            else:
                if end-start<wsize:
                    mid=div((start+end),2)
                    nstart,nend=mid-div(wsize,2),mid+div(wsize,2)
                    if nstart<0:nstart,nend=0,wsize
                else:
                    nstart,nend=start,end
                id+=1
                lks[id],bn[id]={},0
                pklines.append('\t'.join([cr,str(start),str(end),str(nstart),str(nend)]))
                sBin,eBin=div(nstart,binSize),div(nend,binSize)
                for bin in range(sBin,eBin+1):
                    if bin not in pks[cr]:pks[cr][bin]={}
                    pks[cr][bin][nstart]=[nend,id]
    #print 'peaks count:',id
    lks[0]={}
    ids=list(lks.keys())
    ids.sort()
    tn=0
    for line in infile:
        tn+=0.5
        col=line[:-1].split('\t')
        cr1,cr2,pos1,pos2=col[2],col[6],int(col[3]),int(col[7])
        if cr2=='=':cr2=cr1
        bin1,bin2=div(pos1,binSize),div(pos2,binSize)
        id1,id2=[],[]
        for start in pks[cr1][bin1]:
            if pos1>start and pos1<pks[cr1][bin1][start][0]:id1.append(pks[cr1][bin1][start][1])
        for start in pks[cr2][bin2]:
            if pos2>start and pos2<pks[cr2][bin2][start][0]:id2.append(pks[cr2][bin2][start][1])
        if len(id2)==0:id2=[0]
        if len(id1)==0:id1=[0]
        for fid in id1:
            bn[fid]+=0.5
            for tid in id2:
                #if pklines[fid]>pklines[tid]:continue
                #elif pklines[fid]==pklines[tid]:adv=0.5#the value to be added 
                #bn[fid]+=0.5
                if fid!=tid:bn[tid]+=0.5
                if fid>tid:
                    temp=tid
                    tid=fid
                    fid=temp
                if tid not in lks[fid]:lks[fid][tid]=0
                lks[fid][tid]+=0.5
        #if ln-tn-temp>0.5:
        #    temp=ln-tn
        #if len(id1)>1 or len(id2)>1:print tn,ln,cr1,cr2,pos1,pos2,id1,id2,pklines[fid],pklines[tid]
    #tn=int(tn)
    #print 'total edges:',ln,ln*2
    print('total edges:',tn,sum(bn.values()))#tn may be smaller than the sum of values in bn or lks, due to the fact that there are overlapped bins in bn keys 
    
    lkf=open(linkfile,'w')
    lkf.write('chrA\tstartA\tendA\twindowAstart\twindowAend\tedgeA\tchrB\tstartB\tendB\twindowBstart\twindowBend\tedgeB\tedgeAB\texpected\tlog10P\n')
    ids2=deepcopy(ids)
    for fk in ids:
        tks=list(lks[fk].keys())
        tks.sort()
        for tk in ids2:#lks[fk]:
            if tk in lks[fk]:ov=lks[fk][tk]
            else:ov=0
            pv=str(fisherTest(ov,bn[fk],bn[tk],tn)).split()[-1]
            exp=div(bn[tk]*bn[fk]*1.0,tn)
            lkf.write('\t'.join(    [pklines[fk],str(int(bn[fk])),pklines[tk],str(int(bn[tk])),str(int(ov)),str(exp),str(pv)]   )+'\n')
    lkf.close()
    #print tn
    return tn

def mergeLinks(inlinkfile,samFile,outlinkfile,binSize=1000000,wsize=500,wstep=100,fold=100,logP=-100):
    pf=open(inlinkfile)
    pf.readline()
    fdColID,pvColID=9,10
    peaks={}
    for line in pf:
        col=line.split()
        if float(col[fdColID])<fold:continue
        if float(col[pvColID])>logP:continue
        cr,start,end=col[0],int(col[1]),int(col[2])
        if cr not in peaks:peaks[cr]={}#dic[col[0]]
        peaks[cr][start]=end
    peaks=merge_peaks_by_head_tail_distance(peaks,distance=0)
    fo=open(outlinkfile,'w')
    for cr in peaks:
        for start in peaks[cr]:fo.write(cr+'\t'+str(start)+'\t'+str(peaks[cr][start])+'\n')
    fo.close()
    translocationLinks(peakFile=outlinkfile,samFile=samFile,linkfile=outlinkfile,binSize=binSize,wsize=0,wstep=0,fold=0,logP=0)
def peaks2bed(pfile,bfile=None,flank=100):
    if bfile==None:bfile=pfile[:-3]+'bed'
    pf=open(pfile)
    bf=open(bfile,'w')
    pf.readline()
    for line in pf:
        col=line.split()
        cr,ps=col[0],[int(col[1]),int(col[2])]
        for tp in col[3].split(','):ps.append(int(tp))
        ps.sort()
        bf.write('\t'.join([cr,str(max(ps[0]-flank,0)),str(ps[-1]+flank),'1','1','+'])+'\n')
def transPair(linkFile):
    lf=open(linkFile)
    lf.readline()
    dic={}
    for line in lf:
        col=line.split()
        s1,s2='\t'.join(col[:3]),'\t'.join(col[6:9])
        if s1 not in dic:dic[s1]={}
        dic[s1][s2]=[float(col[12]),float(col[15])]
    '''
    for s1 in dic:
        for s2 in dic[s1]:
            if s1!=s2: print s1+'\t'+s2
    '''
    return dic
def pksDic(peakFile,dic={}):
    dic={}
    pf=open(peakFile)
    pf.readline()
    for line in pf:
        col=line.split()
        cr,start,end=col[0],int(col[1]),int(col[2])
        if cr not in dic:dic[cr]={}
        if start not in dic[cr]:dic[cr][start]=end
        else:
            if dic[cr][start]<end:dic[cr][start]=end
    return dic

if __name__ == "__main__":
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # This allow DANPOS to print each message on screen immediately.
    print('')

    