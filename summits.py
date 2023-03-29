#!/usr/bin/env python
from time import time
from copy import deepcopy
import numpy

# since this code was originally written in
# python2 which has integer division
# we much do division differently if both
# dividing numbers are ints


def div(a, b):
    if isinstance(a, int) and isinstance(b, int):
        return(a//b)
    else:
        return(a/b)


class Summits:
    def __init__(self):
        # self.data[chr]['feature']=numpy_array, 'feature' can be 'p':summit position, 'v': summit occupancy value, 's': positioning score, or 'ppos':positioning p value
        self.data = {}

    def fetchValueFromWig(self, wig):
        '''
        Description:
            retrieve occupancy value from a Wig class instance

        Parameter:
            wig: an Wig class instance

        Value:
            None
        '''
        smt = self
        for cr in smt.data:
            if 'v' not in smt.data[cr]:
                smt.data[cr]['v'] = numpy.array([0.0])
                smt.data[cr]['v'].resize(smt.data[cr]['p'].size, refcheck=0)
            ps = smt.data[cr]['p']
            # w=wig.data[cr]
            lth = ps.size
            i = 0
            while(i < lth):
                tp = div(ps[i], wig.step)
                if cr in wig.data:
                    if tp > wig.data[cr].size:
                        # issue: https://groups.google.com/g/danpos/c/dPxLBr_g6Bc
                        wig.data[cr].resize(int(tp)+1, refcheck=0)
                    smt.data[cr]['v'][i] = wig.data[cr][int(tp)]
                else:
                    smt.data[cr]['v'][i] = 1
                i += 1

    def positioning(self, wig, rd=None):
        '''
        Description:
            Calculate nucleosome (summit) positioning score and P value.

        Parameter:
            wig: a Wig class instance containing nucleosome occupancy data
            rd: the half size of the region flanking each summit used to calculate nucleosome fuzziness

        Value:
            None
        '''
        from rpy2.robjects import r, FloatVector
        from random import shuffle
        from functions import log10fuztest
        # ct=r('''function(x){return(chisq.test(x)$p.value)}''')
        # xx=numpy.array(range(rd*2))
        # ftest=r('''function(x,y){return(var.test(x,y)$p.value)}''')
        smt = self
        step = wig.step
        # rd=rd/step+1
        # tarray=numpy.array(0.0)
        # tarray.resize(rd*2,refcheck=0)
        for cr in smt.data:
            print(cr)
            ps = smt.data[cr]['p']
            s, ppos = numpy.array([0.0]*rd*2), numpy.array([0.0]*rd*2)
            lth = ps.size
            s.resize(lth, refcheck=0)
            ppos.resize(lth, refcheck=0)
            if cr in wig.data:  # add by kaifu on Jan 27, 2014
                w = wig.data[cr]
                i = 0
                while(i < lth):
                    temp = log10fuztest(
                        pc=ps[i], pt=ps[i], cr=cr, cwig=wig, twig=None, rd=rd)
                    ppos[i], s[i] = temp
                    i += 1
            smt.data[cr]['s'] = s
            smt.data[cr]['ppos'] = ppos

    def merge(self, wg, distance=110):
        '''
        Description:
            Merge neighboring summits whos distance from each other shorter than a specified distance.

        Parameter:
            wig: a Wig class instance containing nucleosome occupancy data.
            distance: merge neighboring summits whos distance from each other shorter than this value.

        Value:
            None.
        '''

        step = wg.step
        ctime = time()
        tnum = 0
        onum = 0
        for cr in self.data:
            print(cr, ":", end=' ')
            ps = self.data[cr]['p']  # summits positions
            vs = self.data[cr]['v']  # summits values
            onum += ps.size  # original number of summits
            print(ps.size, "summits, merging ...")
            merge = 1
            while(merge > 0):
                merge = 0
                nps = numpy.array([0])
                nvs = numpy.array([0.0])
                # ps=dic.keys()
                # ps.sort()
                lth = ps.size-2
                nps.resize(lth+2, refcheck=0)
                nvs.resize(lth+2, refcheck=0)
                if lth < 0:
                    continue
                i = 0
                ni = 0
                while i < lth:
                    td = ps[i+1]-ps[i]
                    if td >= distance:
                        nps[ni], nvs[ni] = ps[i], vs[i]
                        ni += 1
                    else:
                        merge += 1
                        td2 = ps[i+2]-ps[i+1]
                        if td2 < td:
                            nps[ni], nvs[ni] = ps[i], vs[i]
                            ni += 1
                        else:
                            if vs[i] > vs[i+1]:
                                vs[i+1] = vs[i]
                                ps[i+1] = ps[i]
                            elif vs[i] == vs[i+1]:
                                pos = div((ps[i]+ps[i+1]), 2)
                                ps[i+1] = pos
                                # print(div(pos,wg.step))
                                # added by Kaifu Chen Jul 10,2012 ######
                                vs[i+1] = wg.data[cr][int(div(pos, wg.step))]
                    i += 1
                if (ps[-1]-ps[-2]) >= distance:
                    nps[ni], nps[ni+1], nvs[ni], nvs[ni +
                                                     1] = ps[-2], ps[-1], vs[-2], vs[-1]
                    ni += 2
                else:
                    if vs[-2] > vs[-1]:
                        nps[ni], nvs[ni] = ps[-2], vs[-2]
                    elif vs[-2] == vs[-1]:
                        ###### added by Kaifu Chen Jul 10,2012 ######
                        nps[ni] = div((ps[-2]+ps[-1]), 2)
                        nvs[ni] = wg.data[cr][int(
                            div(div((ps[-2]+ps[-1]), 2), wg.step))]
                    else:
                        nps[ni], nvs[ni] = ps[-1], vs[-1]
                    ni += 1
                    merge += 1
                ps = nps[:ni]
                vs = nvs[:ni]
            print(ps.size, 'left')
            self.data[cr]['p'] = ps
            self.data[cr]['v'] = vs
        return True

    def fillgap(self, wg, height=5, distance=110):
        '''
        Description:
            Insert a summit between two neighboring summits whose distance from each other larger than a specified range.

        Parameter:
            hight: require the nucleosome occupancy of the inserted summit to be higher than or equal than this value.
            distance: require the distance between neighboring summits to be larger than distance*2.5 and smaller than distance*3.5 to allow insertion.

        Value:
            None
        '''

        step = wg.step
        midis, madis = distance*2.5, distance*3.5
        ctime = time()
        for cr in self.data:
            id = []
            np = []
            nv = []
            print(cr, ":", end=' ')
            ps = self.data[cr]['p']  # summits positions
            vs = self.data[cr]['v']  # summits values
            lth = ps.size-1
            i = 0
            while(i < lth):
                gs = ps[i+1]-ps[i]
                if gs > midis and gs < madis:
                    p = int(div((ps[i+1]+ps[i])*0.5, step))
                    # print gs, wg.data[cr][p],wg.data[cr][ps[i]/step],wg.data[cr][ps[i+1]/step]
                    if wg.data[cr][int(p)] >= height:
                        id.append(i+1)
                        np.append(p*step)
                        nv.append(wg.data[cr][int(p)])
                i += 1
            alth = len(id)  # total count of insert smts
            print(alth)
            nlth = lth+1+alth  # the final smts count
            for k in self.data[cr]:
                self.data[cr][k] = deepcopy(self.data[cr][k])
                self.data[cr][k].resize(nlth, refcheck=0)
            j = 0
            while j < alth:
                # the list has been extend j step backward in the previous j cycling
                i = id[j]+j
                # the reason to do this is, there might be some other lists in addition to the 'p' and 'v' lists
                for k in self.data[cr]:
                    self.data[cr][k][i+1:nlth] = self.data[cr][k][i:(nlth-1)]
                    self.data[cr][k][i] = 0
                self.data[cr]['p'][i] = np[j]
                self.data[cr]['v'][i] = nv[j]
                j += 1
        return True


if __name__ == "__main__":
    print('')
