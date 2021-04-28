import os
import stat
import sys
import logging
import subprocess as sp
import ROOT
import itertools
import glob
import getpass
import array
from Utils.ColoredPrintout import colors
from collections import defaultdict
import getopt # command line parser
import argparse
import CombineHarvester.CombineTools.plotting as plot #Combine plotting utils
from EFTPlotter import BuildScan

def test():

    f = open("exclusions.txt", "w")

    counter_SMexclusion = 0
    counter_SMexclusion_secondMin = 0
    counter_secondMin_preferred = 0
    counter_total_toys=0
    counter_tmp = 0

    for i in range(1,50+1):
    # for i in range(1,2):
        print('-- i', i)

        filename = "higgsCombine.1ktoys_v3_cpt_part{}.POINTS.0.39.MultiDimFit.mH120.*root".format(i)
        #print('== Opening file', filename)
        #print(glob.glob(filename))
        filename = glob.glob(filename)[0]
        rootFile = ROOT.TFile.Open(filename)
        print('== Opened file', filename)
        limitTree = rootFile.Get('limit')

        #for itoy in range(1,10+1):
        #for itoy in range(1,20+1):
        for itoy in range(1,10+1):
            print('-- itoy', itoy)
            counter_total_toys+= 1

            #-- Use CombineTool utils (see: https://github.com/cms-analysis/CombineHarvester/blob/master/CombineTools/python/plotting.py)
            graph = plot.TGraphFromTree(limitTree, "cpt", '2*deltaNLL', 'quantileExpected > -1.5 && iToy=={}'.format(itoy))

            yvals = [3.84] #1sigma, 95%CL intervals
            #func, crossings, val, val_2sig, cross_1sig, cross_2sig, other_1sig, other_2sig = BuildScan(graph, ROOT.kBlack, yvals)
            main_scan = BuildScan(graph, ROOT.kBlack, yvals)

            minz = min(graph.GetY())

            if minz < 0: #Always indicates second min...?
                counter_secondMin_preferred+= 1

            if graph.Eval(0) - minz >= 4:
                counter_SMexclusion+= 1
                print('SM EXCLUDED: toy', itoy, 'file', filename)
                f.write('SM EXCLUDED: plot scan1D_1ktoys_cpt_part{}_POINTS_0_49_toy{}.png\n'.format(i,itoy))
                if minz<-4:
                # if minz<0:
                    counter_SMexclusion_secondMin+= 1 #Due to second deeper well
            
            if graph.Eval(0) - minz >= 2: counter_tmp+= 1

            '''
            sm_rel = 0
            min_nll = 0
            for ientry in range(limitTree.GetEntries()): #Loop just once, and compare entries 0 and 1
                # print('=== Entry', ientry, '===')

                limitTree.GetEntry(ientry)
                cpt = limitTree.GetLeaf("cpt").GetValue(0)
                iToy = limitTree.GetLeaf("iToy").GetValue(0)
                deltaNLL = limitTree.GetLeaf("deltaNLL").GetValue(0)
                quantileExpected = limitTree.GetLeaf("quantileExpected").GetValue(0)
                # print('cpt', cpt)
                # print('deltaNLL', deltaNLL)
                # print('quantileExpected', quantileExpected)
                if iToy != itoy: continue
                if deltaNLL < min_nll: min_nll = deltaNLL
                if cpt==-0.5: sm_rel = 2*deltaNLL
                # if cpt==-0.5 or cpt==+0.5: sm_rel=deltaNLL
                # elif quantileExpected==-1: sm_rel-= 2*deltaNLL
                # else: continue


                # for ibr in range(len(list_br)):
                #     br_name = list_br.At(ibr).GetName()
                #     # print(br_name, ' = ', limitTree.GetLeaf(br_name).GetValue(ientry))
                #     limitTree.GetEntry(0)
                #     val0 = limitTree.GetLeaf(br_name).GetValue(0)

            sm_rel-= 2*min_nll
            print('sm_rel', sm_rel)
            print('min_nll', min_nll)

            if sm_rel>=3:
                counter_SMexclusion+= 1
                print('EXCLUDED: toy', itoy, 'file', filename)
            '''

    print('== counter_SMexclusion:', counter_SMexclusion)
    f.write('== counter_SMexclusion: {}\n'.format(counter_SMexclusion))

    print('== counter_SMexclusion:', counter_SMexclusion_secondMin)
    # f.write('== counter_SMexclusion_secondMin: {}\n'.format(counter_SMexclusion_secondMin))

    print('== counter_secondMin_preferred:', counter_secondMin_preferred)

    f.write('== counter_total_toys: {}\n'.format(counter_total_toys))

    f.close()

    print('== counter_tmp:', counter_tmp)

    return

if __name__ == "__main__":

    test()
