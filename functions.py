#!/usr/bin/env python
from time import time
import numpy,glob,os
from wigs import Wigs
from wig import Wig
from summits import Summits
from reads import reads
from copy import deepcopy
from random import randint
from math import log10,log,sqrt
from rpy2.robjects import r, FloatVector
import re,sys
# to replace having to use R in python
from scipy.stats import poisson,f
import pdb

# python2 divides ints differently
# than python3. Since this code was
# originally written in python2 we
# correct the code using the following
# two commands:

def isint(x):
    if type(x).__module__ == "numpy":
        return(x.dtype == numpy.integer)
    else:
        return(isinstance(x,int))


def div(a,b):
    if isint(a) and isint(b):return(a//b)
    else:return(a/b)


# This converts a numpy number to a standard type
def unnumpy(foo):
    if type(foo).__module__ == "numpy":
        foo = foo.item()
    return(foo)

def ppois(q: int, mean: float) -> float:
    """
    Calculates the logarithmic survival function of the Poisson distribution.
    equivalent of ppois function with lower.tail = False and log.p = True in R

    Args:
        q (int): vector of quantiles
        mean (float): (AKA lambda) vector of (non-negative) means

    Returns:
        float: The logarithmic survival function value.
    """
    logsf_value=poisson.logsf(q,mean)
    return(logsf_value)

def pf(q, df1, df2):
    """
    Calculates the base-10 logarithm of the cumulative distribution function of the F-distribution.
    equivalent of pf function with lower.tail = True and log.p = True in R

    Args:
        q (float): vector of quantiles
        df1 (float): The degrees of freedom of the numerator.
        df2 (float): The degrees of freedom of the denominator.

    Returns:
        float: The logarithm of the CDF value.
    """
    logcdf_value = f.logcdf(q, dfn=df1, dfd=df2)
    log10_cdf_value = logcdf_value / numpy.log10(numpy.e)
    return log10_cdf_value


def danpos(tpath=None, tbg=None, opath='./', nonzero=False, amount=None, nor='Q',
            region_file=None, test="P", save=True, fdr=0, pheight=1e-5, height=5,
            logp=5, call_position=1, ref_position=None, both=True, width=40,
            distance=165, fill_gap=False, fill_value=1, edge=1, fcut=None, rd=None,
            ratio=0.9, call_peak=0, ref_peak=None, peak_distance=3000, peak_width=40,
            call_region=0, ref_region=None, region_distance=3000, region_width=40,
            lmd=300, cut=0, wgfmt='fixed', step=10, exclude_low_percent=1, exclude_high_percent=1,
            fs=None, mifrsz=50, mafrsz=200, extend=74, smooth_width=40, pcfer=0, paired=0):
    """
    Description:
        This is the main function that accepts all input parameters and calls other functions in the danpos package to complete the work.

    Parameters:
        Please see the help messages for DANPOS by command lines "python danpos.py -h".
    """

    # Step 1: Parameter checking and processing
    starttime = time()
    rd = distance
    hdiff = 1 - ratio

    fcut = float(numpy.std(list(range(0 - rd, rd + step, step)))) * (0.95 - div(step, (2.0 * rd)))
    print('rd', rd, ', step', step, ', fcut,', fcut)

    pairs, groups = pathParser(tpath=tpath)
    if len(groups) == 0:
        return False

    bggroups, subpairs = bgPathParser(tbg=tbg, groups=groups)
    if tbg != None and (len(bggroups) == 0 or len(subpairs) == 0):
        return False

    scalepairs = scaleParser(amount=amount, extend=extend, groups=groups)
    if amount != None and len(scalepairs) == 0:
        return False

    while opath[-1] == '/':
        opath = opath[:-1]
    if len(opath) < 1:
        opath = 'result'
    elif opath[-1] == '.':
        opath += '/result'
    if not os.path.isdir(opath):
        os.mkdir(opath)
    addname = '.'

    # Step 1 completed.


    # Step 2: load MNase-seq data sets
    maxgroupsize = 0
    for groupname in groups:
        groups[groupname]=loadinput(groups[groupname], fs=fs, cut=cut, save=False,
                                    wgfmt=wgfmt, step=step, extend=extend, mifrsz=mifrsz,
                                    mafrsz=mafrsz, paired=paired, starttime=starttime)

        if maxgroupsize < len(list(groups[groupname].keys())):
            maxgroupsize = len(list(groups[groupname].keys()))

    # Save raw data if nexessary
    if save:
        if (tbg is not None) or (maxgroupsize > 1) or (smooth_width > 1) or (amount is not None) or ((nor != 'N') and ((len(groups) > 1) or (maxgroupsize > 1))):
            if not os.path.isdir(os.path.join(opath,'raw')):
                os.mkdir(os.path.join(opath, 'raw'))
            print('\nsaving raw data')
            for groupname in groups:
                groups[groupname].save(os.path.join(opath,'raw'), format=wgfmt, step=step)
            print('time elapsed:', time()-starttime,'seconds')
    # Step 2 completed


    # Step 3 load genomic background data sets
    if tbg is not None:
        for bggroupname in bggroups:
            bggroups[bggroupname] = loadinput(bggroups[bggroupname], fs=100, cut=cut,
                            save=False, wgfmt=wgfmt, step=step, extend=100, mifrsz=mifrsz,
                            mafrsz=mafrsz, paired=paired, starttime=starttime)

        # Save raw background data if nexessary
        if save:
            print('\nsaving raw background data')
            for bggroupname in bggroups:
                bggroups[bggroupname].save(os.path.join(opath,'raw'), format=wgfmt, step=step)
            print('time elapsed:', time()-starttime,'seconds')
    # Step 3 completed


    # Step 4: genomic background (genomic input) subtraction
    if tbg is not None:
        print('\nsubtracting input effects... ')
        addname+='bgsub.'
        if not os.path.isdir(os.path.join(opath,'pooled')):
            os.mkdir(os.path.join(opath,'pooled'))

        # Compute pooled background groups
        pooledbggroups={}
        for groupname in subpairs:
            print(groupname)
            bggroupname = subpairs[groupname]
            if bggroupname is not 'None':
                if bggroupname not in pooledbggroups:
                    bgfilenames=list(bggroups[bggroupname].keys())
                    if len(bgfilenames)>1:
                        bggroups[bggroupname].nor(nor='F', exclude_low_percent=0,
                                        exclude_high_percent=0, scalepairs=None,
                                        sampling_total=None, nonzero=False)
                    pooledbggroups[bggroupname] = bggroups[bggroupname].pop(bgfilenames[0])
                    if len(bgfilenames)>1:
                        for bgfilename in bgfilenames[1:]:
                            pooledbggroups[bggroupname].add(bggroups[bggroupname].pop(bgfilename))
                        pooledbggroups[bggroupname].foldChange(div(1.0,len(bgfilenames)))
                    if not os.path.isfile(os.path.join(opath, 'pooled', bggroupname+".wig")):
                        pooledbggroups[bggroupname].save(os.path.join(opath, 'pooled', bggroupname+".wig"),
                                                         format=wgfmt, step=step)

            # Subtract background from each experimental group
            for filename in list(groups[groupname].keys()):
                newfilename = filename[:-3]+'bgsub.wig'
                temp = groups[groupname].pop(filename)
                if bggroupname!='None':
                    temp.bgsub(pooledbggroups[bggroupname], lmd=lmd)
                groups[groupname].set(newfilename,temp)
                temp=''
        print('time elapsed:', time()-starttime,'seconds')

        if save:
            if (maxgroupsize > 1) or (smooth_width > 1) or (amount != None) or ((nor != 'N') and ((len(groups) > 1) or (maxgroupsize > 1))):
                print('\nsaving background subtracted data')
                if (not os.path.isdir(os.path.join(opath,"bgsub"))):
                    os.mkdir(os.path.join(opath,"bgsub"))
                for groupname in groups:
                    groups[groupname].save(os.path.join(opath, "bgsub"), format=wgfmt, step=step)
                print('time elapsed:', time()-starttime, 'seconds')
        pooledbggroups=''
        bggroups=''
    # Step 4 completed.

    # Step 5: normalization
    norname={'Q':'Qnor.','S':'Snor.','F':'Fnor.','N':''}
    if len(groups) < 2:
        groupname = list(groups.keys())[0]
        if len(list(groups[groupname].keys())) < 2 and amount == None:
            save = 0
            nor = 'N'
    if (nor != 'N') and ((len(groups) > 1) or (maxgroupsize > 1) or amount != None):
        print('\nnormalizing wigs...')
        addname+=norname[nor]
        wigs = Wigs(step = step)
        for groupname in groups:
            for filename in list(groups[groupname].keys()):
                newname = filename[:-3] + norname[nor] + 'wig'
                groups[groupname].set(newname,groups[groupname].pop(filename))
                wigs.set(newname,groups[groupname].get(newname))
                if amount != None:
                    scalepairs[newname] = scalepairs[groupname]
        if amount == None:
            scalepairs = None
        sampling_total = wigs.samplingTotal(region_file = region_file,
                                            region_out_file = os.path.join(opath, 'normalize_region.wig'),
                                            exclude_low_percent=exclude_low_percent,
                                            exclude_high_percent=exclude_high_percent)
        if region_file == None and sampling_total != None:
            region_file = os.path.join(opath, 'normalize_region.wig')
        sucnor = wigs.nor(nor = nor, scalepairs = scalepairs,
                          sampling_total = sampling_total, nonzero = nonzero)
        if sucnor != 1:
            return sucnor
        print('time elapsed:', time() - starttime,'seconds')

        if save:
            if (maxgroupsize > 1) or (smooth_width > 1) or (amount != None):
                if not os.path.isdir(os.path.join(opath, 'nor')):
                    os.mkdir(os.path.join(opath, 'nor'))
                print("\nsaving normalized wigs...")
                wigs.save(os.path.join(opath, 'nor'), format=wgfmt, step=step)
                print('time elapsed:', time() - starttime, 'seconds')
        wigs=''
    # Step 5 completed.

    # Step 6: replicate smoothing
    if (smooth_width > 1):#if not save, the smoothing will be done after pooling; else it will be done here
        addname+='smooth.'
        # smoothing of individual samples will be done only if necessary
        if save and (maxgroupsize > 1):
            print('\nsmoothing individual samples...')
            for groupname in groups:
                for filename in list(groups[groupname].keys()):
                    groups[groupname].get(filename).smooth(smooth_width)
                    groups[groupname].set(filename[:-3] + 'smooth.wig', groups[groupname].pop(filename))
                if not os.path.isdir(os.path.join(opath, 'smooth')):
                    os.mkdir(os.path.join(opath,'smooth'))
                groups[groupname].save(os.path.join(opath, 'smooth'), format=wgfmt, step=step)
            print('time elapsed:', time()-starttime,'seconds')
    # Step 6 completed.

    # Step 7: replicate position calling
    # Currently not supported
    '''
    if pcfer:
        if (len(groups) > 1) or (maxgroupsize > 1):
            if not os.path.isdir(opath+'/replicate_positions'):
                os.mkdir(opath+'/replicate_positions')
            for groupname in groups:
                for filename in groups[groupname].keys():
                    groups[groupname].get(filename).callPositions(os.path.join(opath, "replicate_positions/" + filename[:-3] + "positions.xls"),
                                                                 width = width, distance = distance, edge = edge, fill_gap = fill_gap,
                                                                 fill_value = fill_value, pcut = pheight, height = height, calculate_P_value = 1, poscal = 1)
    '''
    # Step 7 completed.

    # Step 8: pooling
    pooledgroups = {}
    all_wig = all_wig_format(path = tpath)
    for groupname in groups:
        print("\npooling group", groupname, "...")
        if not os.path.isdir(os.path.join(opath, 'pooled')):
            os.mkdir(os.path.join(opath, 'pooled'))
        filenames = list(groups[groupname].keys())
        pooledgroups[groupname] = groups[groupname].pop(filenames[0])
        if len(filenames) > 1:
            for filename in filenames[1:]:
                pooledgroups[groupname].add(groups[groupname].pop(filename))
            pooledgroups[groupname].foldChange(div(1.0, len(filenames)))
        if (not save) and (smooth_width>1):
            pooledgroups[groupname].smooth(smooth_width)#if not save, the smoothing will be done here; else it have been done before replicate position calling
        if (len(filenames)>1) or (smooth_width>1) or (cut!='0') or (nor!='N') or (tbg!=None) or (amount!=None) or (not all_wig):
            pooledgroups[groupname].save(os.path.join(opath, 'pooled', groupname+addname+"wig"), format=wgfmt, step=step)
        print('time elapsed:', time()-starttime,'seconds')
    groups=''
    # Step 8 complete


    ###### step 11 start --- differential testing ######
    testname={'C':'.chisq_diff.','F':'.fold_diff.','P':'.pois_diff.','S':'.sub_diff.','L':'.lsub_diff.'}
    dfgroups={}
    for i in range(len(pairs)):
        if len(pairs[i])<2:continue
        if (pairs[i][0]=='None') or (pairs[i][1]=='None'):continue
        dfname=':'.join(pairs[i])
        print('\ndifferential test for '+dfname+'...')
        dfgroups[dfname]=pooledgroups[pairs[i][0]].dfTest(test=test,cwig=pooledgroups[pairs[i][1]])
        #if len(dfgroups)>0:
        if not os.path.isdir(os.path.join(opath,'diff')):os.mkdir(os.path.join(opath,'diff'))
        for dfname in dfgroups:
            dfgroups[dfname].save(
                os.path.join(opath,'diff',re.sub(':','-',dfname+testname[test]+"wig")),
                format=wgfmt,step=step)
        print('time elapsed:', time()-starttime,'seconds')
    #else:print 'no comparison was specified to be done'
    ###### step 11 end --- differential testing --- ######

    ###### step 12 start --- peak calling --- ######

    peaks={}
    if (call_region==1 ) or (call_peak==1):
        if call_peak==1:printhead='\npeak'
        else:printhead='\nregion'
        peakgroups,dfpeakgroups,dfpeakgroups2={},{},{}
        for groupname in pooledgroups:
            print(printhead,'calling for',groupname,'...')
            peakgroups[groupname]=pooledgroups[groupname].callRegions(
                ofile=None,width=0,
                distance=0,
                pheight=pheight,
                height=height,
                calculate_P_value=1,
                mode='w',title_line=1,pos_only=False)
            print('time elapsed:', time()-starttime,'seconds')

        if len(dfgroups)>0 and logp!=0:
            #print('\ncalling differential regions...')
            #print('call for gaining...')
            for dfname in dfgroups:
                print(printhead,'calling for',dfname,'gaining...')
                dfpeakgroups[dfname]=dfgroups[dfname].callRegions(ofile=None,width=0,distance=0,pheight=1,height=logp,calculate_P_value=0,mode='w',title_line=1,pos_only=True)
                print('time elapsed:', time()-starttime,'seconds')
            for dfname in dfgroups:
                print(printhead,'calling for',dfname,'loss...')
                dfgroups[dfname].foldChange(-1)
                dfpeakgroups2[dfname]=dfgroups[dfname].callRegions(ofile=None,width=0,distance=0,pheight=1,height=logp,calculate_P_value=0,mode='w',title_line=1,pos_only=True)
                dfgroups[dfname].foldChange(-1)
                print('time elapsed:', time()-starttime,'seconds')

        if len(pooledgroups)>1:
            if call_peak==1:print('\nmerging peaks from all groups...')
            else:print('\nmerging regions from all groups...')
        peak_group_list=[peakgroups]
        if len(dfgroups)>0 and logp!=0:peak_group_list+=[dfpeakgroups,dfpeakgroups2]
        for temp_peak_group in peak_group_list:#merge positions from all group into regions
            for name in temp_peak_group:
                temp_peaks=temp_peak_group[name]
                for chr in temp_peaks:
                    if chr not in peaks:peaks[chr]={}
                    for pos in temp_peaks[chr]:
                        if pos not in peaks[chr]:peaks[chr][pos]=temp_peaks[chr][pos]
                        elif peaks[chr][pos]<temp_peaks[chr][pos]:peaks[chr][pos]=temp_peaks[chr][pos]
    if call_peak==1:
        if ref_peak!=None:
            print('reading reference peaks from',ref_peak)
            peaks=readPeaks(ref_peak)
        elif len(pooledgroups)>1:
            print('define reference peaks by pooling peaks defined in all samples...')
            peaks=merge_peaks_by_head_tail_distance(peaks=peaks,distance=peak_distance)
            '''
            fo=open(os.path.join(opath,'reference_peaks.xls'),'w')
            fo.write('chr\tstart\tend\n')
            for chr in peaks:
                starts=peaks[chr].keys()
                starts.sort()
                for start in starts:fo.write(chr+'\t'+str(start)+'\t'+str(peaks[chr][start])+'\n')
            '''


        print('\nretriving peak values for each group...')

        for groupname in pooledgroups:
            print(groupname)
            temp_peaks=merge_peaks_by_head_tail_distance(peaks=peakgroups[groupname],distance=peak_distance)
            pooledgroups[groupname].fillRegions(
                regions=temp_peaks,
                file=os.path.join(opath,'pooled',groupname+addname+"peaks.xls"),
                pheight=pheight,
                height=height,
                width=peak_width,
                calculate_P_value=1,
                pos_only=False)
            if len(pooledgroups)>1:pooledgroups[groupname].fillRegions(regions=peaks,file=os.path.join(opath,'pooled',groupname+addname+"refpeaks.xls"),pheight=pheight,height=height,width=peak_width,calculate_P_value=1,pos_only=False)
        if len(dfgroups)>0:
            for dfname in dfgroups:
                print(dfname)
                if logp!=0:temp_peaks=merge_peaks_by_head_tail_distance(peaks=dfpeakgroups[dfname],distance=peak_distance)
                if logp!=0:dfgroups[dfname].fillRegions(regions=temp_peaks,file=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"local_gain.peaks.xls")),pheight=pheight,height=logp,width=peak_width,calculate_P_value=0,pos_only=True)
                dfgroups[dfname].fillRegions(regions=peaks,file=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"local_gain.refpeaks.xls")),pheight=pheight,height=logp,width=peak_width,calculate_P_value=0,pos_only=True)
                dfgroups[dfname].foldChange(-1)
                if logp!=0:temp_peaks=merge_peaks_by_head_tail_distance(peaks=dfpeakgroups2[dfname],distance=peak_distance)
                if logp!=0:dfgroups[dfname].fillRegions(regions=temp_peaks,file=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"local_loss.peaks.xls")),pheight=pheight,height=logp,width=peak_width,calculate_P_value=0,pos_only=True)
                dfgroups[dfname].fillRegions(regions=peaks,file=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"local_loss.refpeaks.xls")),pheight=pheight,height=logp,width=peak_width,calculate_P_value=0,pos_only=True)
                dfgroups[dfname].foldChange(-1)
        if call_region!=1: print('time elapsed:', time()-starttime,'seconds')
    ###### step 12 start --- peak calling --- ######




    ###### step 12 start --- peaks to regions --- ######
    if call_region==1:
        if ref_region!=None:
            print('reading reference regions from',ref_region)
            regions=readPeaks(ref_region)
        elif len(pooledgroups)>1:
            print('define reference regions by pooling regions defined in all sample...')
            regions=merge_peaks_by_head_tail_distance(peaks=peaks,distance=region_distance)
            '''
            fo=open(os.path.join(opath,'reference_regions.xls'),'w')
            fo.write('chr\tstart\tend\n')
            for chr in regions:
                starts=regions[chr].keys()
                starts.sort()
                for start in starts:fo.write(chr+'\t'+str(start)+'\t'+str(regions[chr][start])+'\n')
            '''
        print('\nretriving region values for each group...')
        for groupname in pooledgroups:
            print(groupname)
            temp_peaks=merge_peaks_by_head_tail_distance(peaks=peakgroups[groupname],distance=region_distance)
            pooledgroups[groupname].fillRegions(regions=temp_peaks,file=os.path.join(opath,'pooled',groupname+addname+"regions.xls"),pheight=pheight,height=height,width=region_width,calculate_P_value=1,pos_only=False)
            if ref_region!=None or len(pooledgroups)>1:pooledgroups[groupname].fillRegions(regions=regions,file=os.path.join(opath,'pooled',groupname+addname+"refregions.xls"),pheight=pheight,height=height,width=region_width,calculate_P_value=1,pos_only=False)
        if len(dfgroups)>0:
            for dfname in dfgroups:
                print(dfname)
                if logp!=0:temp_peaks=merge_peaks_by_head_tail_distance(peaks=dfpeakgroups[dfname],distance=region_distance)
                if logp!=0:dfgroups[dfname].fillRegions(regions=temp_peaks,file=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"local_gain.regions.xls")),pheight=pheight,height=logp,width=region_width,calculate_P_value=0,pos_only=True)
                dfgroups[dfname].fillRegions(regions=regions,file=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"local_gain.refregions.xls")),pheight=pheight,height=logp,width=region_width,calculate_P_value=0,pos_only=True)
                dfgroups[dfname].foldChange(-1)
                if logp!=0:temp_peaks=merge_peaks_by_head_tail_distance(peaks=dfpeakgroups2[dfname],distance=region_distance)
                if logp!=0:dfgroups[dfname].fillRegions(regions=regions,file=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"local_loss.regions.xls")),pheight=pheight,height=logp,width=region_width,calculate_P_value=0,pos_only=True)
                dfgroups[dfname].fillRegions(regions=regions,file=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"local_loss.refregions.xls")),pheight=pheight,height=logp,width=region_width,calculate_P_value=0,pos_only=True)
                dfgroups[dfname].foldChange(-1)
        print('time elapsed:', time()-starttime,'seconds')
    ###### step 12 end --- peaks to regions --- ######



    ###### step 12 start --- position calling --- ######
    if call_position==1:
        if call_region==0 and call_peak==0: peaks=None
        smtgroups={}
        positionFiles=[]
        for groupname in pooledgroups:
            print('\nposition calling for',groupname,'...')
            positionFiles.append(os.path.join(opath,'pooled',groupname+addname+"positions.xls"))
            smtgroups[groupname]=pooledgroups[groupname].callPositions(os.path.join(opath,'pooled',groupname+addname+"positions.xls"),width=width,distance=distance,edge=0,fill_gap=fill_gap,fill_value=fill_value,pcut=pheight,height=height,calculate_P_value=0,poscal=1,regions=peaks,rd=rd)
            print('time elapsed:', time()-starttime,'seconds')
        if ref_position!=None:
            print('\nreading reference positions from',ref_position)
            refdic=readPositions(ref_position)
        else:
            if len(pooledgroups)>1:print('\ndefining a set of reference positions by pooling positions defined in all samples...')
            refdic=refPositions(positionFiles=positionFiles,distance=distance)
            if len(pooledgroups)>1:
                fo=open(os.path.join(opath,'reference_positions.xls'),'w')
                fo.write('chr\tpos\n')
                for chr in refdic:
                    #starts=refdic[chr].keys()
                    #starts.sort()
                    i,lth=0,len(refdic[chr]['p'])
                    while i<lth:
                        fo.write(chr+'\t'+str(refdic[chr]['p'][i])+'\n')#+str(refdic[chr]['v'][i])+'\n')
                        i+=1

        for groupname in pooledgroups:
            if len(pooledgroups)<2:continue
            print('\nfine-tuning positions for',groupname,'by comparing to the reference map...')
            tdic=positionAdjust(dic=refdic,infile=os.path.join(opath,'pooled',groupname+addname+"positions.xls"),outfile=None,wg=pooledgroups[groupname],distance=distance,fcut=fcut,hdiff=hdiff)
            pooledgroups[groupname].fillPositions(dic=tdic,file=os.path.join(opath,'pooled',groupname+addname+"positions.ref_adjust.xls"),width=width,distance=distance,edge=edge,pcut=pheight,height=height,calculate_P_value=1,mode='w',title_line=1,poscal=1,rd=rd)

        if len(dfgroups)>0:
            dsmtgroups,dsmtgroups2={},{}
            if logp!=0:
                #print('\ncalling differential positions...')
                for dfname in dfgroups:
                    print('\nposition calling for',dfname,'gaining...')
                    dsmtgroups[dfname]=dfgroups[dfname].callPositions(re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"gain.positions.xls")),width=width,distance=distance,edge=edge,fill_gap=0,fill_value=fill_value,pcut=1,height=logp,calculate_P_value=0,poscal=0,regions=peaks)
                    print('time elapsed:', time()-starttime,'seconds')
                for dfname in dfgroups:
                    print('\nposition calling for',dfname,'loss...')
                    dfgroups[dfname].foldChange(-1)
                    dsmtgroups2[dfname]=dfgroups[dfname].callPositions(re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"loss.positions.xls")),width=width,distance=distance,edge=edge,fill_gap=0,fill_value=fill_value,pcut=1,height=logp,calculate_P_value=0,poscal=0,regions=peaks)
                    dfgroups[dfname].foldChange(-1)
                    print('time elapsed:', time()-starttime,'seconds')
            '''
            else:
                for dfname in dfgroups:
                    print '\nretriving position value for',dfname,'gaining...'
                    dsmtgroups[dfname]=dfgroups[dfname].fillPositions(dic=refdic,file=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"gain.positions.xls")),width=width,distance=distance,edge=edge,pcut=1,      height=logp,  calculate_P_value=0,mode='w',title_line=1,poscal=0)
                    print 'time elapsed:', time()-starttime,'seconds'
                for dfname in dfgroups:
                    print '\nretriving position value for',dfname,'loss...'
                    dfgroups[dfname].foldChange(-1)
                    dsmtgroups[dfname]=dfgroups[dfname].fillPositions(dic=refdic,file=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"loss.positions.xls")),width=width,distance=distance,edge=edge,pcut=1,      height=logp,  calculate_P_value=0,mode='w',title_line=1,poscal=0)
                    dfgroups[dfname].foldChange(-1)
                    print 'time elapsed:', time()-starttime,'seconds'
            '''
    ###### step 12 end --- position calling --- ######



    ###### step 13 start --- map differential positions to nucleosome positions --- ######
    if len(dfgroups)>0:
        for dfname in dfgroups:
            #print '\nintegrative comparison for'+dfname
            groupnames=dfname.split(':')
            #print 'doing for',groupnames[0],'and',groupnames[1],'...'
            tgs=0
            for cr in dfgroups[dfname].data:tgs+=dfgroups[dfname].data[cr].size
            fdrsimu=min(100000,max(10000,div(tgs,1000)))
            if call_position==1:
                print('\nposition level integrative analysis for', dfname, '...')
                #if fdr==1:
                print('FDR simulation...')
                occFDRlist=occFDR(dwig=dfgroups[dfname],simu=fdrsimu,regions=peaks)
                fuzFDRlist=fuzFDR(cwig=pooledgroups[groupnames[1]],twig=pooledgroups[groupnames[0]],simu=fdrsimu,rd=rd,regions=peaks)
                #else:occFDRlist,fuzFDRlist=None,None
                print('analyzing...')
                if logp!=0:gainPositionFile,lossPositionFile=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"gain.positions.xls")),re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"loss.positions.xls"))
                else:gainPositionFile,lossPositionFile=None,None
                allPositionsInOneFile(controlPositionFile=os.path.join(opath,'pooled',groupnames[1]+addname+"positions.ref_adjust.xls"),\
                                  treatPositionFile=os.path.join(opath,'pooled',groupnames[0]+addname+"positions.ref_adjust.xls"),\
                gainPositionFile=gainPositionFile,lossPositionFile=lossPositionFile,\
                ofile=re.sub(':','-',os.path.join(opath,dfname+".positions.integrative.xls")),\
                cwig=pooledgroups[groupnames[1]],twig=pooledgroups[groupnames[0]],dwig=dfgroups[dfname],\
                dis=distance,rd=rd,test=test,fdrsimu=0,fdrRegions=peaks,occFDRlist=occFDRlist,fuzFDRlist=fuzFDRlist,fdr=fdr)
                print('time elapsed:', time()-starttime,'seconds')
            if call_peak==1:
                print('\npeak level integrative analysis for', dfname, '...')
                print('FDR simulation...')
                fdrlist=peakFDR(peakFile1=os.path.join(opath,'pooled',groupnames[1]+addname+"refpeaks.xls"),\
                            wg1=pooledgroups[groupnames[1]],wg2=pooledgroups[groupnames[0]],fdrsimu=fdrsimu,cut=height)
                print('analyzing...')
                region_differential(file1=os.path.join(opath,'pooled',groupnames[1]+addname+"refpeaks.xls"),\
                                    file2=os.path.join(opath,'pooled',groupnames[0]+addname+"refpeaks.xls"),\
                                    gainFile=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"local_gain.refpeaks.xls")),\
                                    lossFile=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"local_loss.refpeaks.xls")),\
                ofile=re.sub(':','-',os.path.join(opath,dfname+".peaks.integrative.xls")),widthFDRlist=fdrlist[0],\
                smtFDRlist=fdrlist[1],aucFDRlist=fdrlist[2],step=step,fdr=fdr)
                print('time elapsed:', time()-starttime,'seconds')
            if call_region==1:
                print('\nregion level integrative analysis for', dfname, '...')
                print('FDR simulation...')
                fdrlist=peakFDR(peakFile1=os.path.join(opath,'pooled',groupnames[1]+addname+"refregions.xls"),\
                            wg1=pooledgroups[groupnames[1]],wg2=pooledgroups[groupnames[0]],fdrsimu=fdrsimu,cut=height)
                print('analyzing...')
                region_differential(file1=os.path.join(opath,'pooled',groupnames[1]+addname+"refregions.xls"),\
                                    file2=os.path.join(opath,'pooled',groupnames[0]+addname+"refregions.xls"),\
                                    gainFile=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"local_gain.refregions.xls")),\
                                    lossFile=re.sub(':','-',os.path.join(opath,'diff',dfname+testname[test]+"local_loss.refregions.xls")),\
                ofile=re.sub(':','-',os.path.join(opath,dfname+".regions.integrative.xls")),widthFDRlist=fdrlist[0],\
                smtFDRlist=fdrlist[1],aucFDRlist=fdrlist[2],step=step,fdr=fdr)
                print('time elapsed:', time()-starttime,'seconds')
    ###### step 13 end --- map differential positions to nucleosome positions --- ######

    seconds=int(time()-starttime)
    hours=div(seconds,3600)
    minutes=div((seconds-hours*3600),60)
    seconds=seconds-hours*3600-minutes*60
    print('\ntotal time elapsed:',hours,'hours',minutes,'minutes',seconds,"seconds\n\njob done, cheers!\n\n")
def pathParser(tpath):
    groups={}
    pairs=tpath.split(',')
    if tpath==None:
        print('No input files detected\n')
        return [{},{}]
    else:
        for i in range(len(pairs)):
            pairs[i]=pairs[i].split(':')
            if len(pairs[i])>2:
                print('Wrong: each pair can not contain more than 2 groups, please also make sure that each group name does not contain the symbol "-"')
                return [{},{}]
            for j in range(len(pairs[i])):
                if pairs[i][j]=='None':continue
                if (not os.path.isfile(pairs[i][j])) and (not os.path.isdir(pairs[i][j])):
                    print('Wrong: file or directory',pairs[i][j],'does not exists')
                    return [{},{}]
                group=pairs[i][j]
                while group[-1]=='/':
                    group=group[:-1]
                # groupname=re.sub('/+','_',group) #old naming scheme having PATH_TO_FILE_SAMPLE
                groupname=group.split("/")[-1] #new naming scheme having just SAMPLE
                groupname=re.sub('.gz$','',groupname)
                while groupname[0] in ['.','_']:groupname=groupname[1:]
                if os.path.isfile(group):
                    if groupname[-6:]=='bowtie':groupname=groupname[:-7]
                    else:groupname=groupname[:-4]
                pairs[i][j]=groupname
                if groupname in groups:
                    if groups[groupname]!=group:
                        print('Wrong: different group (',groups[groupname],'and',group,') were found for same group name (',groupname,')')
                        return [{},{}]
                else:groups[groupname]=group
        return [pairs,groups]
def bgPathParser(tbg,groups=None):
    bggroups,subpairs={},{}
    if tbg!=None:
        bgpairs,bggroup=tbg.split(','),''
        for i in range(len(bgpairs)):
            bgpairs[i]=bgpairs[i].split(':')
            if len(bgpairs[i])>2:
                print('Wrong: each pair can not contain more than 2 groups, please also make sure that each group name does not contain the symbol "-"')
                return [{},{}]
            for j in range(len(bgpairs[i])):
                bggroup=bgpairs[i][j]
                if bggroup=='None':continue
                if (not os.path.isfile(bgpairs[i][j])) and (not os.path.isdir(bgpairs[i][j])):
                    print('Wrong: file or directory',bgpairs[i][j],'does not exists')
                    return [{},{}]
                while bggroup[-1]=='/':bggroup=bggroup[:-1]
                bggroupname=re.sub('/+','_',bggroup)
                bggroupname=re.sub('.gz$','',bggroupname)
                while bggroupname[0] in ['.','_']:bggroupname=bggroupname[1:]
                if os.path.isfile(bggroup):
                    if bggroupname[-6:]=='bowtie':bggroupname=bggroupname[:-7]
                    else:bggroupname=bggroupname[:-4]
                bgpairs[i][j]=bggroupname
                if (j==1) and (bggroup!='None'):
                    if bggroupname in bggroups:
                        if bggroups[bggroupname]!=bggroup:
                            print('Wrong: different group (',bggroups[bggroupname],'and',bggroup,') were found for same group name (',bggroupname,')')
                            return [{},{}]
                    else:bggroups[bggroupname]=bggroup
        #print bgpairs
        if len(bgpairs)==1 and len(bgpairs[0])==1:#all data set use the same background input
            if groups!=None:
                bggroups[bgpairs[0][0]]=bggroup
                for name in groups:subpairs[name]=bgpairs[0][0]
        else:
            for i in range(len(bgpairs)):
                if bgpairs[i][0] not in subpairs:subpairs[bgpairs[i][0]]=bgpairs[i][1]
                elif subpairs[bgpairs[i][0]]!=bgpairs[i][1]:
                    print('Wrong: different genomic background groups (', subpairs[bgpairs[i][0]],'and',bgpairs[i][1], ') were found for the same input group (',bgpairs[i][0],')\n')
                    return [{},{}]
        #if len(subpairs)!=len(groups):
        #    print 'count of groups for background subtraction does not equal to count of input groups\n'
        #    return [{},{}]
    #print subpairs
    return [bggroups,subpairs]

def scaleParser(amount,extend,groups):
    scalepairs={}
    if amount!=None:
        mpairs=amount.split(',')
        if len(mpairs)==1:#all data set normalized to the same amount
            temp=mpairs[0].split(':')
            if len(temp)==1:
                mpairs[0]=temp
                mpairs[0][0]=int(mpairs[0][0])*extend
                for name in groups:
                    scalepairs[name]=mpairs[0][0]
                return scalepairs
        for i in range(len(mpairs)):
            print(mpairs[i])
            mpairs[i]=mpairs[i].split(':')
            if len(mpairs[i])>2:
                print('\n!!!!! Wrong: format error, please make sure that name does not contain the symbol "-"!!!!!\n')
                return {}
            try:mpairs[i][1]=float(mpairs[i][1])
            except:
                print('\n!!!!! Wrong: amount is not readable for',mpairs[i][0],'!!!!!\n')
                return
            if (not os.path.isfile(mpairs[i][0])) and (not os.path.isdir(mpairs[i][0])):
                print('Wrong: file or directory',mpairs[i][0],'does not exists')
                return {}
            group=mpairs[i][0]
            while group[-1]=='/':
                group=group[:-1]
            # groupname=re.sub('/+','_',group) #old naming scheme having PATH_TO_FILE_SAMPLE
            groupname=group.split("/")[-1] #new naming scheme having just SAMPLE

            while groupname[0] in ['.','_']:groupname=groupname[1:]

            if os.path.isfile(group):
                if groupname[-2:]=='gz':groupname=groupname[:-3]
                if groupname[-6:]=='bowtie':groupname=groupname[:-7]
                else:groupname=groupname[:-4]
            mpairs[i][0]=groupname
            if mpairs[i][0] not in scalepairs:scalepairs[mpairs[i][0]]=mpairs[i][1]*extend
            elif scalepairs[mpairs[i][0]]!=mpairs[i][1]*extend:
                print('Wrong: different expectation of amounts (', scalepairs[mpairs[i][0]],'and',mpairs[i][1], ') were found for the same group (',mpairs[i][0],')\n')
                return {}
        if len(scalepairs)!=len(groups):
            print('count of groups for scaling does not equal to count of input groups\n')
            return {}
    return scalepairs

def region_differential(file1,file2,gainFile,lossFile,ofile=None,widthFDRlist=None,smtFDRlist=None,aucFDRlist=None,step=10,fdr=0):
    #r['options'](warn=-1)
    #prop=r('''function(vec){prop.test(matrix(vec,nrow=2))$p.value}''')
    f1,f2,gf,lf=open(file1).readlines(),open(file2).readlines(),open(gainFile).readlines(),open(lossFile).readlines()
    #print len(f1),len(f2),len(gf),len(lf)
    #f1.readline()
    #f2.readline()
    if ofile!=None:
        fo=open(ofile,'w')
        fo.write('chr\tstart\tend\tcenter\tcontrol_height\ttreat_height\theight_log2FC\theight_diff_log10Pval\theight_diff_FDR\tcontrol_width\ttreat_width\twidth_log2FC\twidth_diff_log10Pval\twidth_diff_FDR\tcontrol_total_signal\ttreat_total_signal\ttotal_signal_log2FC\ttotal_signal_diff_log10Pval\ttotal_signal_diff_FDR\tlocal_gain_log10Pval\tlocal_gain_FDR\tlocal_loss_log10Pval\tlocal_loss_FDR\tlocal_change_log10Pval\tlocal_change_FDR\n')
    i,lth=1,len(f1)
    lw,ls,lt=len(widthFDRlist),len(smtFDRlist),len(aucFDRlist)
    while i<lth:
        col1,col2,gcol,lcol=f1[i].split(),f2[i].split(),gf[i].split(),lf[i].split()
        #f1[i]=col1
        if col1[:3]!=col2[:3]:print('Wrong:',col1[:3],col2[:3])
        w1,w2=int(col1[4]),int(col2[4])
        w=int(col1[2])-int(col2[1])
        t1,t2=div(float(col1[5])*step,w),div(float(col2[5])*step,w)
        s1,s2=float(col1[6]),float(col2[6])#summit height
        if w1==w2:wdiff =0
        else:wdiff=log10PropTest([w1,w2,w2,w1])
        if t1==t2:tdiff=0
        else:tdiff=log10PropTest([t1,t2,t2,t1])
        if s1==s2:sdiff=0
        else:sdiff=float((ppois(unnumpy(max(s1,s2)+1),unnumpy(min(s1,s2)+1))/log(10)).split()[-1])
        gdiff,ldiff=float(gcol[6]),float(lcol[6])# (0-log10Pval) of gain and loss
        if fdr==1:
            wfdr = div(findRank(widthFDRlist,0-wdiff)*1.0,lw)
            tfdr = div(findRank(aucFDRlist,0-tdiff)*1.0,lt)
            sfdr = div(findRank(smtFDRlist,0-sdiff)*1.0,ls)
            gfdr = div(findRank(smtFDRlist,gdiff)*1.0,ls)
            lfdr = div(findRank(smtFDRlist,ldiff)*1.0,ls)
        else:wfdr,tfdr,sfdr,gfdr,lfdr='-','-','-','-','-'
        #print '\t'.join([str(wdiff),str(tdiff),str(sdiff)])
        if ofile!=None:
            try:fo.write('\t'.join(
                    col1[:4]+[str(s1),
                    str(s2),
                    str(div(log(div(max(s2,1),max(s1,1))),log(2))),
                    str(sdiff),
                    str(sfdr),
                    str(w1),
                    str(w2),
                    str(div(log(div(max(w2*1.0,1.0),max(w1,1))),log(2))),
                    str(wdiff),
                    str(wfdr),
                    str(t1),
                    str(t2),
                    str(div(log(div(max(t2,1),max(t1,1))),log(2))),
                    str(tdiff),
                    str(tfdr),
                    str(0-gdiff),
                    str(gfdr),
                    str(0-ldiff),
                    str(lfdr),
                    str(0-max(gdiff,ldiff)),
                    str(min(gfdr,lfdr))])+'\n')#+dcol[4:]+col1[4:-1]+col2[4:-1])+'\n')
            except: print(w2,w1,t2,t1,s2,s1)
        #print '\t'.join(col1[:4]+[str(wdiff),str(wfdr),str(tdiff),str(tfdr),str(sdiff),str(sfdr)])
        i+=1
def peakFDR(peakFile1,peakFile2=None,wg1=None,wg2=None,fdrsimu=1000000,cut=5):
    #r['options'](warn=-1)
    pk=open(peakFile1).readlines()
    if peakFile2!=None:pk+=open(peakFile2).readlines()[1:]
    i,lth=1,len(pk)
    step=wg1.step
    while i<lth:
        col=pk[i].split()
        col[1]=div(int(col[1]),step)
        col[2]=div(int(col[2]),step)
        pk[i]=col[:3]
        i+=1
    chrs={}
    for cr in wg1.data:
        if cr in wg2.data :chrs[cr]=min(wg1.data[cr].size,wg2.data[cr].size)-1
    i,lth=0,lth-1
    w,s,t=numpy.resize(numpy.array([0.0]),fdrsimu),numpy.resize(numpy.array([0.0]),fdrsimu),numpy.resize(numpy.array([0.0]),fdrsimu)
    while i<fdrsimu:
        if i%1000==0:print(i,'simulated')
        j=randint(1,lth)
        cr,rstart,rend,width=pk[j][0],pk[j][1],pk[j][2],pk[j][2]-pk[j][1]
        start=randint(rstart,rend)
        end=start+width
        if end>chrs[cr]-width:start,end=chrs[cr]-width,chrs[cr]
        w1=0
        w2=0
        s1=wg1.data[cr][int(start):int(end)].max()
        s2=wg2.data[cr][int(start):int(end)].max()
        t1=div(wg1.data[cr][int(start):int(end)].sum()*step,width)
        t2=div(wg2.data[cr][int(start):int(end)].sum()*step,width)
        k=start
        while k<end:
            if wg1.data[cr][int(k)]>=cut:w1+=step
            if wg2.data[cr][int(k)]>=cut:w2+=step
            k+=1
        if w1==w2:wdiff =0
        else:wdiff=log10PropTest([w1,w2,w2,w1])
        if t1==t2:tdiff=0
        else:tdiff=log10PropTest([t1,t2,t2,t1])
        if s1==s2:sdiff=0
        else:sdiff=float((ppois(unnumpy(max(s1,s2)+1),unnumpy(min(s1,s2)+1))/log(10)).split()[-1])
        w[i],s[i],t[i]=wdiff,sdiff,tdiff
        i+=1
        #print wdiff,sdiff,tdiff
    w.sort()
    s.sort()
    t.sort()
    return [0-w,0-s,0-t]


def log10PropTest(list=[]):
    '''
    Parameter:
        list: a list contain four number [count_in_A, count_not_in_A, count_in_B, count_not_in_B]
    Return:
        log scaled P value
    '''
    log10PropTest=r('''function (x, correct = TRUE)
    {
        x=matrix(x,nrow=2)
        l <- nrow(x)
        n <- rowSums(x)
        x <- x[, 1L]
        k <- length(x)
        ESTIMATE <- x/n
        correct <- as.logical(correct)
        YATES <- ifelse(correct && (k <= 2), 0.5, 0)
        DELTA <- ESTIMATE[1L] - ESTIMATE[2L]
        YATES <- min(YATES, abs(DELTA)/sum(1/n))
        p <- sum(x)/sum(n)
        PARAMETER <- k - 1
        x <- cbind(x, n - x)
        E <- cbind(n * p, n * (1 - p))
        STATISTIC <- sum((abs(x - E) - YATES)^2/E)
        names(STATISTIC) <- "X-squared"
        PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE,log.p = TRUE)
        return(PVAL/log(10))
    }
    ''')
    return float(str(log10PropTest(
        unnumpy(FloatVector(list)))).split()[-1])




def log10FisherTest(list=[]):
    '''
    Parameter:
        list: a list contain four number [count_in_A, count_not_in_B, overlap, total]
    Return:
        log scaled P value
    '''
    log10PropTest=r('''function (x)
    {
        return( phyper(x[3] - 1, x[1], x[4]-x[1], x[2], lower.tail = FALSE,log.p=TRUE)/log(10) )
    }

    ''')
    return float(str(log10PropTest(unnumpy(FloatVector(list)))).split()[-1])




def merge_peaks_by_head_tail_distance(peaks,distance):
    regions=peaks
    region_distance=distance
    for chr in regions :#merge neighboring regions
        ps=list(regions[chr].keys())
        ps.sort()
        lth=len(ps)-1
        i=0
        while i<lth:
            if ps[i+1]-regions[chr][ps[i]]<=region_distance:# or pks[chr][ps[i+1]]<=pks[chr][ps[i]]:
                if regions[chr][ps[i]]<regions[chr][ps[i+1]]:regions[chr][ps[i]]=regions[chr][ps[i+1]]
                regions[chr].pop(ps[i+1])
                ps[i+1]=ps[i]
            i+=1
    return regions
def readPeaks(file):
    regions={}
    lines=open(file).readlines()
    for line in lines[1:]:
        col=line.split()
        cr,start,end=col[0],int(col[1]),int(col[2])
        if cr not in regions:regions[cr]={}
        regions[cr][start]=end
    return regions

def all_wig_format(path):
    '''
    Description:
        whether all files are in wiggle format? return True is yes, else return False
    Parameters:
        path: a path to the directory of file(s) of sequence reads or occupancy data.

    '''
    #print 'checking whther all input are in wiggle format'
    tpaths=path.split(',')
    for tpath in tpaths:
        for path in tpath.split(':'):
            #print path
            if os.path.isdir(path):
                for infile in glob.glob( os.path.join(path, '*.bed') ):return False
                for infile in glob.glob( os.path.join(path, '*.bowtie') ):return False
                for infile in glob.glob( os.path.join(path, '*.bam') ):return False
                for infile in glob.glob( os.path.join(path, '*.sam') ):return False
                for infile in glob.glob( os.path.join(path, '*.bed.gz') ):return False
                for infile in glob.glob( os.path.join(path, '*.bowtie.gz') ):return False
                for infile in glob.glob( os.path.join(path, '*.bam.gz') ):return False
                for infile in glob.glob( os.path.join(path, '*.sam.gz') ):return False
            elif os.path.isfile(path):
                if path[-3:]=='bed':return False
                if path[-6:]=='bowtie':return False
                elif path[-3:]=='bam':return False
                elif path[-3:]=='sam':return False
                if path[-3:]=='bed.gz':return False
                if path[-6:]=='bowtie.gz':return False
                elif path[-3:]=='bam.gz':return False
                elif path[-3:]=='sam.gz':return False
    return True

def loadinput(path,fs=None,cut=1e-10,save=False,wgfmt='fixed',step=10,extend=100,mifrsz=100,mafrsz=200,paired=0,starttime=None):
    '''
    Description:
        load occupancy data in '.wig' format file, or calculate occupancy from sequencing reads, use sequencing reads in '.bed','.sam', and '.bam' format as input, generate occupancy data in wiggle format.

    Parameters:
        path: a path to the directory of file(s) of sequence reads or occupancy data.
        cut: the cutoff for removing clonal reads, could be P value larger than 0 and small than 1, or count as a positive integer.
        save: set to 'False' if don't need to save the occupancy data in wiggle files.
        nor: the method to normalize the occupancy if multiple replicates are provided in the input directory, 'F': fold change, 'Q':quantile normalization, 'S': resampling, 'N': no normalization to be done.
        wgfmt: the format of the ouput wiggle files, currently support 'fixed' only.
        step: the step size of the occupancy data
        fs: average size of fragments that are subject to sequencing and generate the reads, only for signgle-end reads. When this value is not given, a fs value will be infered by the program.
        extend: a interger value, each read will be extend to this length.
        mifrsz: the minimal estimated average fragment size
        mafrsz: the maximal estimated average fragment size
        paired: is the reads paired-end (set to 1) or single-end (set to 0)

    '''
    wigs={}
    if os.path.isdir(path):
        for infile in glob.glob( os.path.join(path, '*.bed') )+glob.glob( os.path.join(path, '*.bed.gz') ):
            rd=reads(infile,step=step,paired=paired,cut=cut,format='bed')
            fname=re.sub('/+','_',infile)
            fname=re.sub('.gz$','',fname)
            while fname[0] in ['.','_']:fname=fname[1:]
            if paired:wigs[fname[:-4]+'.wig']=rd.toWig(fs=0,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            else:wigs[fname[:-4]+'.wig']=rd.toWig(fs=fs,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            rd=''
            if save: wigs[fname[:-4]+'.wig'].save(fname[:-4]+'.wig',format=wgfmt,step=step)
            if starttime!=None: print('time elapsed:', time()-starttime,'seconds')
        for infile in glob.glob( os.path.join(path, '*.bowtie') )+glob.glob( os.path.join(path, '*.bowtie.gz') ):
            rd=reads(infile,step=step,paired=paired,cut=cut,format='bowtie')
            fname=re.sub('/+','_',infile)
            fname=re.sub('.gz$','',fname)
            while fname[0] in ['.','_']:fname=fname[1:]
            if paired:wigs[fname[:-7]+'.wig']=rd.toWig(fs=0,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            else:wigs[fname[:-7]+'.wig']=rd.toWig(fs=fs,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            rd=''
            if save: wigs[fname[:-7]+'.wig'].save(fname[:-7]+'.wig',format=wgfmt,step=step)
            if starttime!=None: print('time elapsed:', time()-starttime,'seconds')
        for infile in glob.glob( os.path.join(path, '*.bam') ):
            rd=reads(infile,step=step,paired=paired,cut=cut,format='bam')
            fname=re.sub('/+','_',infile)
            while fname[0] in ['.','_']:fname=fname[1:]
            if paired:wigs[fname[:-4]+'.wig']=rd.toWig(fs=0,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            else:wigs[fname[:-4]+'.wig']=rd.toWig(fs=fs,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            rd=''
            if save: wigs[fname[:-4]+'.wig'].save(fname[:-4]+'.wig',format=wgfmt,step=step)
            if starttime!=None: print('time elapsed:', time()-starttime,'seconds')
        for infile in glob.glob( os.path.join(path, '*.sam') ):
            rd=reads(infile,step=step,paired=paired,cut=cut,format='sam')
            fname=re.sub('/+','_',infile)
            while fname[0] in ['.','_']:fname=fname[1:]
            if paired:wigs[fname[:-4]+'.wig']=rd.toWig(fs=0,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            else:wigs[fname[:-4]+'.wig']=rd.toWig(fs=fs,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            rd=''
            if save: wigs[fname[:-4]+'.wig'].save(fname[:-4]+'.wig',format=wgfmt,step=step)
            if starttime!=None: print('time elapsed:', time()-starttime,'seconds')
        for infile in glob.glob( os.path.join(path, '*.wig') )+glob.glob( os.path.join(path, '*.wig.gz')):
            fname=re.sub('/+','_',infile)
            fname=re.sub('.gz$','',fname)
            while fname[0] in ['.','_']:fname=fname[1:]
            wigs[fname]=Wig(infile,step=step)
            if starttime!=None: print('time elapsed:', time()-starttime,'seconds')
    elif os.path.isfile(path):
        if path[-3:]=='bed' or path[-6:]=='bed.gz':
            rd=reads(path,step=step,paired=paired,cut=cut,format='bed')
            fname=re.sub('/+','_',path)
            fname=re.sub('.gz$','',fname)
            while fname[0] in ['.','_']:fname=fname[1:]
            if paired:wigs[fname[:-4]+'.wig']=rd.toWig(fs=0,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            else:wigs[fname[:-4]+'.wig']=rd.toWig(fs=fs,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            if save: wigs[fname[:-4]+'.wig'].save(fname[:-4]+'.wig',format=wgfmt,step=step)
            rd=''
            if starttime!=None: print('time elapsed:', time()-starttime,'seconds')
        if path[-6:]=='bowtie' or path[-9:]=='bowtie.gz':
            rd=reads(path,step=step,paired=paired,cut=cut,format='bowtie')
            fname=re.sub('/+','_',path)
            fname=re.sub('.gz$','',fname)
            while fname[0] in ['.','_']:fname=fname[1:]
            if paired:wigs[fname[:-7]+'.wig']=rd.toWig(fs=0,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            else:wigs[fname[:-7]+'.wig']=rd.toWig(fs=fs,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            if save: wigs[fname[:-7]+'.wig'].save(fname[:-7]+'.wig',format=wgfmt,step=step)
            rd=''
            if starttime!=None: print('time elapsed:', time()-starttime,'seconds')
        elif path[-3:]=='bam':
            rd=reads(path,step=step,paired=paired,cut=cut,format='bam')
            fname=re.sub('/+','_',path)
            while fname[0] in ['.','_']:fname=fname[1:]
            if paired:wigs[fname[:-4]+'.wig']=rd.toWig(fs=0,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            else:wigs[fname[:-4]+'.wig']=rd.toWig(fs=fs,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            if save: wigs[fname[:-4]+'.wig'].save(fname[:-4]+'.wig',format=wgfmt,step=step)
            rd=''
            if starttime!=None: print('time elapsed:', time()-starttime,'seconds')
        elif path[-3:]=='sam':
            rd=reads(path,step=step,paired=paired,cut=cut,format='sam')
            fname=re.sub('/+','_',path)
            while fname[0] in ['.','_']:fname=fname[1:]
            if paired:wigs[fname[:-4]+'.wig']=rd.toWig(fs=0,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            else:wigs[fname[:-4]+'.wig']=rd.toWig(fs=fs,extend=extend,mifrsz=mifrsz,mafrsz=mafrsz)#,fdrwigs=fdrwigs,fdrflag=fdrflag)
            if save: wigs[fname[:-4]+'.wig'].save(fname[:-4]+'.wig',format=wgfmt,step=step)
            rd=''
            if starttime!=None: print('time elapsed:', time()-starttime,'seconds')
        elif path[-3:]=='wig' or path[-6:]=='wig.gz':
            fname=re.sub('/+','_',path)
            fname=re.sub('.gz$','',fname)
            while fname[0] in ['.','_']:fname=fname[1:]
            wigs[fname]=Wig(path,step=step)
            if starttime!=None: print('time elapsed:', time()-starttime,'seconds')

    out=Wigs()
    out.data=wigs
    wigs=out
    out=''
    #wigs.nor(nor=nor)
    #if save and nor!='N':
    #    for k in wigs.keys(): wigs.get(k).save(k[:-3]+str(nor)+"nor.wig",format=wgfmt,step=step)
    return wigs



def neighborPosition(pos,dic,dis):
    for d in range(dis):
        if pos+d in dic:return d
        elif pos-d in dic:return 0-d
    return dis
def combinePositions(controlPositionFile,treatPositionFile,gainPositionFile=None,lossPositionFile=None,dis=75):
    cpd={}
    nc=0
    for line in open(controlPositionFile).readlines()[1:]:
        col=line.split()
        if col[0] not in cpd:cpd[col[0]]={}
        cpd[col[0]][int(col[3])]=-1
        nc+=1

    tpd={}
    nt=0
    for line in open(treatPositionFile).readlines()[1:]:
        col=line.split()
        if col[0] not in tpd:tpd[col[0]]={}
        tpd[col[0]][int(col[3])]=-1
        nt+=1

    c2t={}
    t2c={}
    while nc>0 or nt>0:
        #print nc,nt
        tcpd=deepcopy(cpd)
        ttpd=deepcopy(tpd)
        for cr in tcpd:
            if cr not in ttpd:ttpd[cr]={}
            for pos in tcpd[cr]:
                d=neighborPosition(pos,ttpd[cr],dis)
                if abs(d)<dis:tcpd[cr][pos]=pos+d
        for cr in ttpd:
            if cr not in tcpd:tcpd[cr]={}
            for pos in ttpd[cr]:
                d=neighborPosition(pos,tcpd[cr],dis)
                if abs(d)<dis:ttpd[cr][pos]=pos+d
        for cr in tcpd:
            if cr not in c2t:c2t[cr]={}
            for cpos in tcpd[cr]:
                tpos=tcpd[cr][cpos]
                if tpos in ttpd[cr]:
                    if ttpd[cr][tpos]==cpos:
                        c2t[cr][cpos]=tpos
                        cpd[cr].pop(cpos)
                        nc-=1
                else:
                    c2t[cr][cpos]=cpos
                    cpd[cr].pop(cpos)
                    nc-=1
        for cr in ttpd:
            if cr not in c2t:c2t[cr]={}
            for tpos in ttpd[cr]:
                cpos=ttpd[cr][tpos]
                if cpos in tcpd[cr]:
                    if tcpd[cr][cpos]==tpos:
                        c2t[cr][cpos]=tpos
                        tpd[cr].pop(tpos)
                        nt-=1
                else:
                    c2t[cr][tpos]=tpos
                    tpd[cr].pop(tpos)
                    nt-=1

    t2c={}
    for cr in c2t:
        t2c[cr]={}
        for cpos in c2t[cr]:
            tpos=c2t[cr][cpos]
            t2c[cr][tpos]=cpos
            c2t[cr][cpos]=[tpos]
    for file in [gainPositionFile,lossPositionFile]:
        if file==None:continue
        for line in open(file).readlines()[1:]:
            col=line.split()
            cr,pos=col[0],int(col[3])
            cdis=neighborPosition(pos,c2t[cr],dis)
            tdis=neighborPosition(pos,t2c[cr],dis)
            if abs(tdis)<dis or abs(cdis)<dis:
                if abs(tdis)<abs(cdis):cpos=t2c[cr][pos+tdis]
                else:cpos=pos+cdis
                c2t[cr][cpos].append(pos)
            '''
            else:#add by kaifu on Mar 10, 2013
                if not c2t.has_key(cr):c2t[cr]={}
                c2t[cr][pos]=[pos,pos]
            '''
    return c2t

def allPositionsInOneFile(ofile='result.xls',controlPositionFile=None,treatPositionFile=None,gainPositionFile=None,lossPositionFile=None,\
                          cwig=None,twig=None,dwig=None,dis=75,rd=None,test='P',fdrsimu=0,fdrRegions=None,occFDRlist=None,fuzFDRlist=None,fdr=0):
    if cwig.step!=twig.step:
        print("step values are different in the control and treatment wiggle data, work can not be done")
        return 0
    else:step=cwig.step
    print('combining nucleosome and differential positions...')
    c2tDic=combinePositions(controlPositionFile,treatPositionFile,gainPositionFile,lossPositionFile,dis)
    '''
    cpd={}
    for line in open(controlPositionFile).readlines()[1:]:
        col=line.split()
        if not cpd.has_key(col[0]):cpd[col[0]]={}
        cpd[col[0]][int(col[3])/step]=-1

    tpd={}
    for line in open(treatPositionFile).readlines()[1:]:
        col=line.split()
        if not tpd.has_key(col[0]):tpd[col[0]]={}
        tpd[col[0]][int(col[3])/step]=-1
    dpd={}
    for file in [gainPositionFile,lossPositionFile]:
        if file==None:continue
        for line in open(file).readlines()[1:]:
            col=line.split()
            if not dpd.has_key(col[0]):dpd[col[0]]={}
            dpd[col[0]][int(col[3])/step]=-1
    '''
    cpd,tpd,dpd={},{},{}
    for line in open(controlPositionFile).readlines()[1:]:
        col=line.split()
        if col[0] not in cpd:cpd[col[0]],tpd[col[0]],dpd[col[0]]={},{},{}
        cpd[col[0]][div(int(col[3]),step)]=-1
    for line in open(treatPositionFile).readlines()[1:]:
        col=line.split()
        if col[0] not in tpd:cpd[col[0]],tpd[col[0]],dpd[col[0]]={},{},{}
        tpd[col[0]][div(int(col[3]),step)]=-1
    for file in [gainPositionFile,lossPositionFile]:
        if file==None:continue
        for line in open(file).readlines()[1:]:
            col=line.split()
            if col[0] not in dpd:cpd[col[0]],tpd[col[0]],dpd[col[0]]={},{},{}
            dpd[col[0]][div(int(col[3]),step)]=-1

    print('calculating differential values for positions...')
    #rd=rd/step
    sumc=cwig.sum()
    sumt=twig.sum()
    outf=open(ofile,'w')
    outf.write('\t'.join(['chr','start','end','center','control_smt_loca','treat_smt_loca','diff_smt_loca','treat2control_dis','control_smt_val','treat_smt_val','smt_log2FC','smt_diff_log10pval','smt_diff_FDR',\
                          'control_point_val','treat_point_val','point_log2FC','point_diff_log10Pval','point_diff_FDR','control_fuzziness_score','treat_fuzziness_score','fuzziness_log2FC','fuzziness_diff_log10pval','fuzziness_diff_FDR\n']))
    nuc,olines,dr1,dr2,dr=-1,[],[],[],[]
    for cr in c2tDic:
        if cr in twig.data or cr in cwig.data:
            if cr not in twig.data:twig.data[cr]=deepcopy(cwig.data[cr])*0
            if cr not in cwig.data:cwig.data[cr]=deepcopy(twig.data[cr])*0
            if cr not in dwig.data:dwig.data[cr]=deepcopy(twig.data[cr])*0
        else:continue
        print(cr)
        cposes=list(c2tDic[cr].keys())
        cposes.sort()
        for cpos in cposes:
            p1,p2=cpos,c2tDic[cr][cpos][0] # cpos is the control nucleosome position position, c2tDic[cr][cpos][0] is the treatment nucleosome position position
            if div(p1,step)>=cwig.data[cr].size:continue #add by kaifu on Mar 10, 2013
            if div(p2,step)>=twig.data[cr].size:continue #add by kaifu on Mar 10, 2013
            temp=log10fuztest(pc=p1,pt=p2,cr=cr,cwig=cwig,twig=twig,rd=rd)
            dp,sd1,sd2=temp#,float(r.sd(FloatVector(temp[1]))[0]),float(r.sd(FloatVector(temp[2]))[0])
            nuc+=1

            if test=='P':
                #pdb.set_trace()
                if cwig.data[cr][div(p1, step)] > twig.data[cr][div(p2, step)]:
                    dp1 = float((0-(ppois(unnumpy(cwig.data[cr][div(p1, step)]),
                                          unnumpy(max(twig.data[cr][div(p2, step)], 1)))/log(10))))
                else:
                    dp1 = float((0-(ppois(unnumpy(twig.data[cr][div(p2, step)]),
                                          unnumpy(max(cwig.data[cr][div(p1, step)], 1)))/log(10))))

            if len(c2tDic[cr][cpos])<2:# no differential position is assigned to c2tDic[cr][cpos] (cpos is the control nucleosome position position, c2tDic[cr][cpos][0] is the treatment nucleosome position position)
                minp,maxp,maxv=min(p1,p2,dwig.data[cr].size*step-step,cwig.data[cr].size*step-step,twig.data[cr].size*step-step),max(p1,p2),0 # '*step' is add by kaifu on Mar 6, 2013
                midp=div((minp+maxp),2)
                minp,maxp=max(0,midp-dis),min(midp+dis,dwig.data[cr].size*step-step,cwig.data[cr].size*step-step,twig.data[cr].size*step-step) # '*step' is add by kaifu on Mar 6, 2013
                p=div((minp+maxp),(2*step)) #replace '2' by '2*step' by kaifu on Mar 6, 2013
                tp=div((minp+maxp),(2*step))
                tdis=div((maxp-minp),(2*step))
                while tdis>=0:
                    if abs(maxv)<abs(dwig.data[cr][int(tp+tdis)]):
                        maxv=dwig.data[cr][tp+tdis]
                        p=tp+tdis
                    if abs(maxv)<abs(dwig.data[cr][int(tp-tdis)]):
                        maxv=dwig.data[cr][tp-tdis]
                        p=tp-tdis
                    tdis-=1
                c2tDic[cr][cpos].append(int(p*step))
            p=c2tDic[cr][cpos][1]
            dp2=abs(dwig.data[cr][div(p,step)])
            if len(c2tDic[cr][cpos][1:])>2:
                for tp in c2tDic[cr][cpos][1:]:
                    tdp2=abs(dwig.data[cr][div(tp,step)]) #corrected by Kaifu Chen on Feb 13,2013, previously it was: tdp2=dwig.data[cr][tp/step]
                    if tdp2>dp2:
                        p=tp
                        dp2=tdp2
            dr2.append(dp2)
            dr1.append(dp1)
            dr.append(dp)
            #if not dpd.has_key(cr):strp='-'
            #elif not dpd[cr].has_key(p/step):strp='-'
            #else:
            strp=str(p)
            if cr not in cpd:strp1='-'
            elif div(p1,step) not in cpd[cr]:strp1='-'
            else:strp1=str(p1)
            if cr not in tpd:strp2='-'
            elif div(p2,step) not in tpd[cr]:strp2='-'
            else:strp2=str(p2)
            if strp1=='-' and strp2=='-':strp1,strp2=strp,strp
            elif strp1=='-':strp1==strp2
            elif strp2=='-':strp2=strp1
            middle=div((p1+p2),2)
            olines.append([cr, \
                str(middle-dis), \
                str(middle+dis), \
                str(middle), \
                str(strp1), \
                str(strp2), \
                str(strp), \
                str(abs(p2-p1)), \
                str(cwig.data[cr][div(p1,step)]), \
                str(twig.data[cr][div(p2,step)]), \
                str(div(log(div(max(twig.data[cr][div(p2,step)],1),max(cwig.data[cr][div(p2,step)],1))),log(2))),\
                str(0-dp1), '-', \
                str(cwig.data[cr][div(p,step)]), \
                str(twig.data[cr][div(p,step)]), \
                str(div(log(div(max(1,twig.data[cr][div(p,step)]),max(cwig.data[cr][div(p,step)],1))),log(2))), \
                str(0-dp2),'-', \
                str(sd1), \
                str(sd2), \
                str(div(log(div(max(1,sd2),max(sd1,1))),log(2))), \
                str(dp),'-',nuc])
    if fdr==1:
        if fdrsimu==0:
            tgs=0
            for cr in dwig.data:tgs+=dwig.data[cr].size
            fdrsimu=min(100000,max(10000,div(tgs,1000)))
        print('calculating occupancy differential FDR...')
        if type(occFDRlist).__module__=='builtins' and occFDRlist==None:
            occFDRlist=occFDR(dwig=dwig,simu=fdrsimu,regions=fdrRegions)
        for i in range(len(dr1)):dr1[i]=div(findRank(occFDRlist,dr1[i])*1.0,fdrsimu)
        for i in range(len(dr2)):dr2[i]=div(findRank(occFDRlist,dr2[i])*1.0,fdrsimu)

        print('calculating fuzziness differential FDR...')
        if type(fuzFDRlist).__module__=='builtins' and fuzFDRlist==None:
            fuzFDRlist=fuzFDR(cwig=cwig,twig=twig,simu=fdrsimu,rd=rd,regions=fdrRegions)
        for i in range(len(dr)):dr[i]=div(findRank(fuzFDRlist,0-dr[i])*1.0,fdrsimu)

    for line in olines:
        if fdr==1:line[12],line[17],line[-2]=str(dr1[line[-1]]),str(dr2[line[-1]]),str(dr[line[-1]])
        outf.write('\t'.join(line[:-1])+'\n')
    outf.close()
def fuzFDR(cwig,twig,simu=10000,rd=None,regions=None):
    gs=0
    if regions==None:
        regions={}
        for cr in cwig.data:
            regions[cr]={}
            regions[cr][0]=(cwig.data[cr].size-1)*cwig.step
            gs+=regions[cr][0]
    else:
        for cr in regions:
            for start in regions[cr]:
                gs+=regions[cr][start]-start

    #if simu==0:simu=min(100000,max(gs/100,1000))
    vec=numpy.array([0.0])
    vec.resize(simu,refcheck=0)
    id=0
    for cr in regions:
        for start in regions[cr]:
            count=div(simu*(regions[cr][start]-start),gs)
            rstart,rend=start-50,regions[cr][start]+50
            if rstart<rd+cwig.step:rstart=rd+cwig.step
            if cr not in cwig.data:continue
            if cr not in twig.data:continue
            if rend>cwig.data[cr].size*cwig.step-rd-cwig.step:rend=cwig.data[cr].size*cwig.step-rd-cwig.step
            if rend>twig.data[cr].size*twig.step-rd-twig.step:rend=twig.data[cr].size*twig.step-rd-twig.step
            if rend<=rstart:continue
            i=0
            while i<count:
                p=randint(rstart,rend)
                tempv=log10fuztest(pc=p,pt=p,cr=cr,cwig=cwig,twig=twig,rd=rd)
                #if len(tempv)<3:continue
                #if tempv[0]!=0:fdrv=log(tempv[0])/log(10)
                #else:fdrv=-350
                #if fdrv<-350:fdrv=-350
                vec[id]=tempv[0]
                i+=1
                id+=1
                if id%1000==0:print(id,'simulated...')
    print(id, 'simulated.')
    vec.sort()
    return 0-vec

def log10fuztest(pc,pt,cr,cwig,twig=None,rd=None):
    #ftest=r('''function(x,y){return(var.test(x,y)$p.value)}''')
    # removed this line below since rpy is not efficient
    #pf=r('''function(s,df1,df2){return(pf(s, df1, df2,log.p = TRUE)/log(10))}''')
    step=cwig.step
    bv,bc=0,div(rd*2,step)
    for d in range(-rd,rd+step,step):bv+=d*d
    bvc=div(bv,bc)

    if twig==None:
        if pc>=(cwig.data[cr].size-1)*rd:return [0,sqrt(bvc)]
        if pc<rd: return [0,sqrt(bvc)]
        v,c=var(p=pc,cr=cr,wig=cwig,step=step,rd=rd,bv=bv,bc=bc)
        if v>=bvc:return [0,sqrt(v)]
        else:
            p=pf(unnumpy(div(v,bvc)),
                unnumpy(c),unnumpy(bc))
            return[p,sqrt(v)]
    else:
        if cr not in twig.data:return(fuztest(pc=pc,pt=pt,cr=cr,cwig=cwig,twig=None,rd=rd))
        if cr not in cwig.data:
            temp=fuztest(pc=pc,pt=pt,cr=cr,cwig=twig,twig=None,rd=rd)
            reurn[temp[0],temp[2],temp[1]]
        if pc>=(cwig.data[cr].size-1)*rd:return [0,sqrt(bvc),sqrt(bvc)]
        if pt>=(twig.data[cr].size-1)*rd:return [0,sqrt(bvc),sqrt(bvc)]
        if pc<rd: return [0,sqrt(bvc),sqrt(bvc)]
        if pt<rd: return [0,sqrt(bvc),sqrt(bvc)]
        vc,cc,=var(p=pc,cr=cr,wig=cwig,step=step,rd=rd,bv=bv,bc=bc)
        vt,ct=var(p=pt,cr=cr,wig=twig,step=step,rd=rd,bv=bv,bc=bc)
        if vc<vt:p=pf(unnumpy(div(vc,vt)),unnumpy(cc),unnumpy(ct))
        else:p=pf(unnumpy(div(vt,vc)),
            unnumpy(ct),unnumpy(cc))
        return[p,sqrt(vc),sqrt(vt)]

def var(p,cr,wig,step,rd,bv,bc):
    start,end=0-rd,rd+step
    tp1,v1,c1=div(p,step),bv,bc
    mi=0#wig.data[cr][(tp1+start/step):(tp1+end/step)].min()
    for d in range(start,end,step):
        vd=d*d
        try:
            v1+=vd*(wig.data[cr][tp1+div(d,step)]-mi)
            c1+=wig.data[cr][tp1+div(d,step)]-mi
        except:
            print(p,d,wig.data[cr].size)
            #return [v1/c1,c1]
    return [div(v1,c1),c1]

def occFDR(dwig,simu=1000000,regions=None):
    gs=0
    if regions==None:
        regions={}
        for cr in dwig.data:
            regions[cr]={}
            regions[cr][0]=(dwig.data[cr].size-1)*dwig.step
            gs+=regions[cr][0]
    else:
        for cr in regions:
            for start in regions[cr]:
                gs+=regions[cr][start]-start

    #if simu==0:simu=min(100000,max(gs/100,1000))
    vec=numpy.array([0.0])
    vec.resize(simu,refcheck=0)
    id=0
    for cr in regions:
        for start in regions[cr]:
            count=div(simu*(regions[cr][start]-start),gs)
            rlth=div((regions[cr][start]-start),2)
            rstart,rend=start-rlth,regions[cr][start]+rlth
            if rstart<dwig.step:rstart=dwig.step
            if rend>dwig.data[cr].size*dwig.step-dwig.step:rend=dwig.data[cr].size*dwig.step-dwig.step
            if rend<=rstart:continue
            i=0
            while i<count:
                p=randint(rstart,rend)
                try: vec[id]=0-abs(dwig.data[cr][div(p,dwig.step)])
                except: print('wrong',id,vec.size,p,dwig.data[cr].size)
                id+=1
                i+=1
    vec.sort()
    return 0-vec
def findRank(vec,v):
    '''
    find the rank of v in a list vec
    Note: the vec must have been ranked in decreasing order
    '''
    v=float(v)
    terminal=len(vec)-1
    start,end=0,terminal+1
    find=True
    while find:
        p=div((start+end),2)
        if p>=terminal:return terminal
        if p<=0: return 0
        if vec[p-1]>=v and vec[p]<=v:
            return p
        elif vec[p]>v:start=p
        elif vec[p]<v:end=p

def refPositions(positionFiles,distance=100):
    dic={}
    for file in positionFiles:
        print('reading from',file,'...')
        fi=open(file)
        fi.readline()
        for line in fi:
            col=line.split()
            cr,pos,fuz=col[0],int(col[3]),float(col[5])
            if cr not in dic:dic[cr]={}
            if pos not in dic[cr]:dic[cr][pos]=fuz
            elif dic[cr][pos]<fuz:dic[cr][pos]=fuz
        fi.close()
    for cr in dic:
        ks=list(dic[cr].keys())
        ks.sort()
        vs=numpy.array([0.0])
        vs.resize(len(ks),refcheck=0)
        i=0
        while i <vs.size:
            vs[i]=dic[cr][ks[i]]
            i+=1
        ks=numpy.array(ks)
        dic[cr]={}
        dic[cr]['p']=ks
        dic[cr]['v']=vs



    print('mering...')
    ctime=time()
    tnum=0
    onum=0
    #mergeto={}
    for cr in dic:
        #mergeto[cr]={}
        print(cr,":")
        ps=dic[cr]['p']#summits positions
        vs=dic[cr]['v']# values
        onum+=ps.size# original number of summits
        print(ps.size,"summits, merging...")
        merge=1
        while(merge>0):
            merge=0
            nps=numpy.array([0])
            nvs=numpy.array([0.0])
            #ps=dic.keys()
            #ps.sort()
            lth=ps.size-2
            nps.resize(lth+2,refcheck=0)
            nvs.resize(lth+2,refcheck=0)
            if lth<0:continue
            i=0
            ni=0
            while i<lth:
                td=ps[i+1]-ps[i]
                if td>=distance:
                    nps[ni],nvs[ni]=ps[i],vs[i]
                    ni+=1
                else:
                    merge+=1
                    td2=ps[i+2]-ps[i+1]
                    if td2<td:
                        nps[ni],nvs[ni]=ps[i],vs[i]
                        ni+=1
                    else:
                        if vs[i]<vs[i+1]:
                            vs[i+1]=vs[i]
                            #mergeto[cr][ps[i+1]]=ps[i]
                            ps[i+1]=ps[i]
                        elif vs[i]==vs[i+1]:
                            pos=div((ps[i]+ps[i+1]),2)
                            ps[i+1]=pos
                            #mergeto[cr][ps[i+1]]=pos
                            vs[i+1]=div((vs[i]+vs[i+1]),2)#wg.data[cr][pos/wg.step] ###### added by Kaifu Chen Jul 10,2012 ######
                        #else:mergeto[cr][ps[i]]=ps[i+1]
                i+=1
            if (ps[-1]-ps[-2])>=distance:
                nps[ni],nps[ni+1],nvs[ni],nvs[ni+1]=ps[-2],ps[-1],vs[-2],vs[-1]
                ni+=2
            else:
                if vs[-2]<vs[-1]:
                    nps[ni],nvs[ni]=ps[-2],vs[-2]
                elif vs[-2]==vs[-1]:
                    nps[ni],nvs[ni]=div((ps[-2]+ps[-1]),2),div((vs[-2]+vs[-1]),2)#wg.data[cr][((ps[-2]+ps[-1])/2)/wg.step] ###### added by Kaifu Chen Jul 10,2012 ######
                else:nps[ni],nvs[ni]=ps[-1],vs[-1]
                ni+=1
                merge+=1
            ps=nps[:ni]
            vs=nvs[:ni]
        print(ps.size, 'left')
        dic[cr]['p']=ps
        dic[cr]['v']=vs
    return dic

def readPositions(file):
    dic=[]
    fi=open(file)
    fi.readline()
    for line in fi:
        col=line.split()
        cr,pos=col[0],col[1]
        if cr not in dic:dic[cr]={'p':[]}
        dic[cr]['p'].append(pos)
        #dic[cr]['v'].append(val)
    for cr in dic:
        dic[cr]['p']=numpy.array(dic[cr]['p'])
        #dic[cr]['v']=numpy.array(dic[cr]['v'])
    return dic
def positionAdjust(dic,infile,outfile,wg,distance=100,fcut=None,hdiff=0.1):
    lines=open(infile).readlines()
    if outfile!=None:fo=open(outfile,'w')
    if outfile!=None:fo.write(lines[0])
    tdic={}
    td,btd,step=div(distance,wg.step),div(distance,wg.step)+1,wg.step
    for line in lines[1:]:
        col=line.split()
        cr,pos,fuz=col[0],int(col[3]),float(col[5])
        if cr not in tdic:tdic[cr]={}
        if fuz<=fcut:tdic[cr][pos]=line
        else:
            tp=div(pos,step)
            h = max(wg.data[cr][int(tp)],1)
            mi1 =wg.data[cr][int(max(tp-td,0)):int(tp)].min()
            mi2 = wg.data[cr][int(tp):int(min(tp+btd,wg.data[cr].size))].min()
            mi=max(mi1,mi2)
            if div((h-mi),h)>=hdiff and div((h-mi),h)>=hdiff*1.5:tdic[cr][pos]=line#retuire occupancy to be hdiff fold lower at one linker relative to position center, and 1.5*hdiff lower at the other linker.

        '''
        if mergeto[cr].has_key(pos):
            tpos=mergeto[cr][pos]
            m+=1
            while mergeto[cr].has_key(tpos):tpos=mergeto[cr][tpos]
            if not tdic[cr].has_key(tpos):tdic[cr][tpos]=[]
            if fuz>fcut:
                if len(tdic[cr][tpos])==0:tdic[cr][tpos].append([cr,'-','-',str(tpos),'-',str(fuz)])
            else:
                if len(tdic[cr][tpos])==0:tdic[cr][tpos].append(col)
                elif tdic[cr][tpos][0][1]=='-':tdic[cr][tpos][0]=[cr,'-','-',str(tpos),'-',str(fuz)]
                else:tdic[cr][tpos].append(col)
        else:
            c+=1
            tdic[cr][pos]=[col]
        '''
    for cr in dic:
        if cr not in tdic:tdic[cr]={}
        for pos in dic[cr]['p']:
            if pos not in tdic[cr]:
                d=neighborPosition(pos=pos,dic=tdic[cr],dis=distance)
                if abs(d)>=distance:tdic[cr][pos]='\t'.join([cr,'-','-',str(pos),'-','-'])+'\n'
    for cr in tdic:
        ps=list(tdic[cr].keys())
        ps.sort()
        if outfile!=None:
            for p in ps:
                for line in tdic[cr][p]:fo.write(line)
        tdic[cr]={}
        tdic[cr]['p']=numpy.array(ps)
    if outfile!=None:fo.close()
    return tdic

def scnMerge(wfs,width=0,distance=250,pheight=1,height=5,linkfold=0,linkLogP=0,binSize=1000,wsize=500,wstep=0):
    wd={}
    for file in wfs:wd[file]=Wig(file)
    pkg={}

    for name in wd:pkg[name]=wd[name].callRegions(ofile=None,width=width,distance=distance,pheight=pheight,height=height,calculate_P_value=1,mode='w',title_line=1,pos_only=False,fold=0,suppress=True)
    pks={}
    for name in pkg:
        tpks=pkg[name]
        for cr in tpks:
            if cr not in pks:pks[cr]={}
            for start in tpks[cr]:
                if start not in pks[cr]:pks[cr][start]=tpks[cr][start]
                elif pks[cr][start]<tpks[cr][start]:pks[cr][start]=tpks[cr][start]
    for name in wd:
        wd[name].fillRegions(regions=pks,file=name[:-3]+'peaks.xls',pheight=pheight,height=height,width=width,calculate_P_value=1,pos_only=False)

def scn(tpath,name='result',pdis=3000,step=1,mapq=30,clipSize=3,inter=True,intra=True,saveWig=False,width=0,distance=250,pheight=1,height=5,zscore=3,linkfold=0,linkLogP=0,binSize=1000,wsize=500,wstep=0,exclude_low_percent=0,exclude_high_percent=1,bnum=100000):
    '''
    parameter:
    norto: may be 'mappable','unique','trans'
    '''
    from lib import translocationReads
    from lib import translocationLinks
    if not os.path.exists(name):os.mkdir(name)
    pairs, groups=pathParser(tpath)
    allWigs=deepcopy(groups)
    bindic,readsCount={},{'unmappable':0,'non_unique':0,'unique_inter':0,'unique_intra':0,'unique_other':0,'unique':0,'mappable':0,'trans':0,}#readsCount
    for groupname in groups:
        print('\n',groupname)
        groupfilename,path=os.path.join(name,groupname),groups[groupname]
        open(groupfilename+'.sam','w').close()
        readsCount[groupname]={}
        if os.path.isfile(path):groups[groupname],allWigs[groupname]=translocationReads(path,bindic=bindic,binSize=binSize,outReadsFile=groupfilename+'.sam',outmode='a',pdis=pdis,step=step,mapq=mapq,clipSize=clipSize,inter=inter,intra=intra,readsCount=readsCount[groupname])
        else:
            groups[groupname],allWigs[groupname]=Wig(step=step),Wig(step=step)
            for file in glob.glob(os.path.join(path,'*sam')):
                print('\n'+file)
                temp=translocationReads(file,bindic=bindic,binSize=binSize,outReadsFile=groupfilename+'.sam',outmode='a',pdis=pdis,step=step,mapq=mapq,clipSize=clipSize,inter=inter,intra=intra,readsCount=readsCount[groupname])
                groups[groupname].add(temp[0])
                allWigs[groupname].add(temp[1])
        print('\nreads in group',groupname)
        for k in ['all','mappable','unmappable','non_unique','unique_inter','unique_intra','unique_other','unique','trans']:#readsCount[groupname]:
            print(k,readsCount[groupname][k])
        if saveWig:
            print('')
            groups[groupname].save(groupfilename+'.trans.wig')
            allWigs[groupname].save(groupfilename+'.all.wig')
    uwigs=Wigs(step=step)
    uwigs.data=allWigs
    sampling_total=uwigs.samplingTotal(exclude_low_percent=exclude_low_percent,exclude_high_percent=exclude_high_percent,bnum=bnum)
    tsum=0
    for groupname in groups:tsum+=sampling_total[groupname]
    taverage=div(tsum*1.0,len(list(sampling_total.keys())))
    for groupname in groups:
        tsum=uwigs.data[groupname].sum()
        uwigs.data[groupname].foldChange(div(taverage,sampling_total[groupname]))
        print('\nnormalize',groupname,'wiggle data of all all reads from',tsum,'to',uwigs.data[groupname].sum())
        if saveWig:uwigs.data[groupname].save(os.path.join(name,groupname)+'.all.nor.wig')

    for groupname in groups:
        print('\ncalling for',groupname,'...')
        tsum=groups[groupname].sum()
        groups[groupname].foldChange(div(taverage,sampling_total[groupname]))
        print('normalize from',tsum,'to',groups[groupname].sum())
        if saveWig:groups[groupname].save(os.path.join(name,groupname)+'.trans.nor.wig')
        #The value of groups[groupname] will be setted as the peaks data.
        groups[groupname]=groups[groupname].callRegions(ofile=os.path.join(name,groupname)+'.trans.sites.xls',width=width,distance=distance,pheight=pheight,height=height,calculate_P_value=0,mode='w',title_line=1,pos_only=False,fold=0,suppress=True)

    print('\nmerging peaks by head-tail distance...')
    pks={}
    for tpks in groups:
        for cr in groups[tpks]:
            if cr not in pks:pks[cr]={}
            for start in groups[tpks][cr]:
                if start not in pks[cr]:pks[cr][start]=groups[tpks][cr][start]
                elif pks[cr][start]<groups[tpks][cr][start]:pks[cr][start]=groups[tpks][cr][start]
    pks=merge_peaks_by_head_tail_distance(pks,distance=distance)
    for groupname in groups:
        print('\nlinking',groupname,'...')
        #The value of groups[groupname] will be setted as total count link reads.
        #groups[groupname]=
        translocationLinks(peaks=pks,samFile=os.path.join(name,groupname+'.sam'),linkfile=os.path.join(name,groupname+'.trans.links.xls'),bindic=bindic,fold=0,logP=0,binSize=binSize,wsize=wsize,wstep=wstep)

    if len(pairs)>0:
        for pair in pairs:
            print('\ncomparing',pair,'...')
            #print 'normalize by',norto,'reads'
            diffwig=uwigs.data[pair[0]].ppois(uwigs.data[pair[1]])
            diffwig.save(os.path.join(name,pair[0]+'-'+pair[1]+'.all.cnv.wig'))
            diffwig.callRegions(ofile=os.path.join(name,pair[0]+'-'+pair[1]+'.all.cnv.sites.gain.xls'),width=width,distance=distance,pheight=pheight,height=height,calculate_P_value=0,mode='w',title_line=1,pos_only=True,fold=0,suppress=True)
            diffwig.foldChange(-1)
            diffwig.callRegions(ofile=os.path.join(name,pair[0]+'-'+pair[1]+'.all.cnv.sites.loss.xls'),width=width,distance=distance,pheight=pheight,height=height,calculate_P_value=0,mode='w',title_line=1,pos_only=True,fold=0,suppress=True)
            tn1,tn2,sf1,sf2,lf1,lf2,lf=sampling_total[pair[0]],sampling_total[pair[1]],os.path.join(name,pair[0]+'.trans.sam'),os.path.join(name,pair[1]+'.trans.sam'),os.path.join(name,pair[0]+'.trans.links.xls'),os.path.join(name,pair[1]+'.trans.links.xls'),os.path.join(name,'-'.join(pair)+'.trans.links.xls')
            scnCompare(tn1=tn1,tn2=tn2,sf1=sf1,sf2=sf2,lf1=lf1,lf2=lf2,lf=lf,bindic=bindic,width=width,distance=distance,pheight=pheight,height=height,linkfold=linkfold,linkLogP=linkLogP,binSize=binSize,wsize=wsize,wstep=wstep)

    print('\nall job done, cheers!\n\n')

def scnCompare(tn1,tn2,sf1,sf2,lf1,lf2,lf,pvalue=1e-3,bindic={},width=0,distance=250,pheight=1,height=5,zscore=3,linkfold=0,linkLogP=0,binSize=1000,wsize=500,wstep=0):
    pvalue=div(log(pvalue),log(10))
    tn=div((tn1+tn2),2.0)
    print('normalize reads count from',tn1,'and',tn2,'to',tn)
    nf1,nf2=div(tn,tn1), div(tn,tn2)
    f1,f2,f=open(lf1),open(lf2),open(lf,'w'),
    t1=f1.readline().split()
    t2=f2.readline().split()
    f.write('\t'.join(t1[:5]+t1[6:11]+['observationA','observationB','expectedA','expectedB','log10PA','log10PB','normalized_observationA','normalized_observationB','log2FC','log10Pdiff'])+'\n')
    for l1 in f1:
        l2=f2.readline()
        col1,col2=l1.split(),l2.split()
        wrong=False
        for i in [0,1,2,3,4,6,7,8,9,10]:
            if col1[i]!=col2[i]:wrong=True
        if wrong:
            print('wrong line:',line[:-1])
            continue
        if col1[0]=='-' or col1[6]=='-':continue
        o1,o2,e1,e2,p1,p2=int(col1[12]),int(col2[12]),float(col1[13]),float(col2[13]),float(col1[14]),float(col2[14])#,float(col2[15])
        if o1+o2==0:continue
        pv=log10PropTest(list=[o1,tn1,o2,tn2])
        no1,no2=int(o1*nf1+0.5),int(o2*nf2+0.5)
        #try:
        logFC=div(log(div(max(1,no1)*1.0,max(1,no2))),log(2))
        if pv<pvalue:f.write('\t'.join(col1[:5]+col1[6:11]+[str(o1),str(o2),str(e1),str(e2),str(p1),str(p2),str(no1),str(no2),str(logFC),str(pv)])+'\n')
        #except:print no1,no2

if __name__ == "__main__":
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0) # This allow DANPOS to print each message on screen immediately.
