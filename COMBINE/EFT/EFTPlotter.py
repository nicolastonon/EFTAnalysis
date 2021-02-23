#-- Plot fit results
#-- Adapted from: https://github.com/cms-govner/EFTFit
'''
*** Critical values of chi2 distribution wih N dofs: (cf. PDG18 Table 39.2)
* N = 1
- 68.27% (1sigma) -> 1.00
- 95% -> 3.84
- 95.45% (2sigma) -> 4.00
- 99% -> 6.63
- 99.73% (3sigma) -> 9.00

* N = 2
- 68.27% (1sigma) -> 2.30
- 95% -> 5.99
- 95.45% (2sigma) -> 6.18
- 99% -> 9.21
- 99.73% (3sigma) -> 11.83
'''

import ROOT
from ROOT import TCanvas, TGraph, gStyle, TMath
import logging
import os
import sys
import numpy
import itertools
import subprocess as sp
from ContourHelper import *
from Utils.ColoredPrintout import colors
import getopt # command line parser
import argparse
import math
import numpy as np
from functools import partial
from settings import opts #Custom dictionnary of settings
import CombineHarvester.CombineTools.plotting as plot #Combine plotting utils
import ctypes


 #    # ###### #      #####  ###### #####
 #    # #      #      #    # #      #    #
 ###### #####  #      #    # #####  #    #
 #    # #      #      #####  #      #####
 #    # #      #      #      #      #   #
 #    # ###### ###### #      ###### #    #

#-- Get the X positions at which the log-likelihood function intersects a given y-line (only keep first and last intersections. There may be more if the function has several minima)
def Get_Intersection_X(graph, y_line):

    xmin, xmax, y_tmp = ROOT.Double(0), ROOT.Double(0), ROOT.Double(0) #Necessary to pass by reference in GetPoint()
    graph.GetPoint(0, xmin, y_tmp) #First point
    graph.GetPoint(graph.GetN()-1, xmax, y_tmp) #Last point
    #print('xmin', xmin, 'xmax', xmax)

    list_X_intersects = []
    step = 0.001 #Desired precision in X (if too coarse, intersection lines/fill areas are displaced)
    y_previousPoint = 99999
    searchingNewMinimum = True #Once an intersection has been found, must wait for the curve to 'invert' before we can find another one
    for x in np.arange(xmin, xmax, step):
        y_thisPoint = graph.Eval(x)
        # print('x', x, 'y', y_thisPoint, 'diff', abs(y_thisPoint-y_line), 'comp', abs(y_previousPoint-y_line), 'searchingNewMinimum', searchingNewMinimum)
        if searchingNewMinimum and abs(y_thisPoint-y_line) > abs(y_previousPoint-y_line): #Check if minimum was found at previous point
            # print('!!! Found intersection', x-step)
            list_X_intersects.append(x-step/2.) #Found intersection, taken as mid-point between 2 graph points
            searchingNewMinimum = False
        elif searchingNewMinimum is False and abs(y_thisPoint-y_line) < abs(y_previousPoint-y_line): searchingNewMinimum = True #Curve has inverted, can look for next intersection
        else: #Check next points
            y_previousPoint = y_thisPoint

    '''
    #-- Obsolete
    step = 0.001 #Desired precision in X (if too coarse, intersection lines/fill areas are displaced)
    for x in np.arange(xmin, xmax, step):
        # print('x', x)
        y_tmp = graph.Eval(x)
        #print('x', x, 'y', y_tmp, 'diff', abs(y_tmp-y_line))
        if abs(y_tmp-y_line) < step*10: #Arbitrary threshold to define intersection
            if len(list_X_intersects)>0 and abs(x-list_X_intersects[len(list_X_intersects)-1])<0.1: continue #Avoid returning several contiguous x positions in a row
            list_X_intersects.append(x)
    '''

    if(len(list_X_intersects)==0):
        print(colors.fg.orange + '[Get_Intersection_X] Warning: no intersection found with zz function. This indicates that the fit rootfile may be empty/incorrect, or simply that the precision is not sufficient to exclude any point at the given threshold ' + str(y_line) + ' (<-> no intersection)' + colors.reset)
        return []

    if(len(list_X_intersects)>=2): list_X_intersects = [list_X_intersects[idx] for idx in [0,len(list_X_intersects)-1]] #remove?

    #print(list_X_intersects)
    return list_X_intersects


def Get_Parameter_LegName(name):

    # EFT operators
    if name == 'ctz': return 'C_{tZ} / #Lambda^{2} [TeV^{-2}]'
    elif name == 'ctw': return 'C_{tW} / #Lambda^{2} [TeV^{-2}]'
    elif name == 'cpq3': return 'C^{3}_{#varphiQ} / #Lambda^{2} [TeV^{-2}]'
    elif name == 'cpqm': return 'C^{-}_{#varphiQ} / #Lambda^{2} [TeV^{-2}]'
    elif name == 'cpt': return 'C_{#varphit} / #Lambda^{2} [TeV^{-2}]'

    #if name == 'ctz': return 'C_{tZ} (#Lambda/TeV)^{2}'
    #elif name == 'ctw': return 'C_{tW} (#Lambda/TeV)^{2}'
    #elif name == 'cpq3': return 'C^{3}_{#varphiQ} (#Lambda/TeV)^{2}'
    #elif name == 'cpqm': return 'C^{-}_{#varphiQ} (#Lambda/TeV)^{2}'
    #elif name == 'cpt': return 'C_{#varphit} (#Lambda/TeV)^{2}'

    # SM signal strengths
    elif name == 'r_tzq': return '#mu(tZq)'
    elif name == 'r_ttz': return '#mu(ttZ)'

    return name


#See: https://github.com/cms-analysis/CombineHarvester/blob/ada429b880fd475d5faf3f4e0780957f8d08df66/CombineTools/python/plotting.py#L1032
def FindCrossingsWithSpline(graph, func, yval):
    crossings = []
    intervals = []
    current = None
    for i in xrange(graph.GetN() - 1):
        if (graph.GetY()[i] - yval) * (graph.GetY()[i + 1] - yval) < 0.:
            cross = func.GetX(yval, graph.GetX()[i], graph.GetX()[i + 1])
            if (graph.GetY()[i] - yval) > 0. and current is None:
                current = {
                    'lo': cross,
                    'hi': graph.GetX()[graph.GetN() - 1],
                    'valid_lo': True,
                    'valid_hi': False
                }
            if (graph.GetY()[i] - yval) < 0. and current is None:
                current = {
                    'lo': graph.GetX()[0],
                    'hi': cross,
                    'valid_lo': False,
                    'valid_hi': True
                }
                intervals.append(current)
                current = None
            if (graph.GetY()[i] - yval) < 0. and current is not None:
                current['hi'] = cross
                current['valid_hi'] = True
                intervals.append(current)
                current = None
            # print 'Crossing between: (%f, %f) -> (%f, %f) at %f' %
            # (graph.GetX()[i], graph.GetY()[i], graph.GetX()[i+1],
            # graph.GetY()[i+1], cross)
            crossings.append(cross)
    if current is not None:
        intervals.append(current)
    if len(intervals) == 0:
        current = {
            'lo': graph.GetX()[0],
            'hi': graph.GetX()[graph.GetN() - 1],
            'valid_lo': False,
            'valid_hi': False
        }
        intervals.append(current)
    print 'intervals', intervals
    return intervals
    # return crossings


#See: https://github.com/cms-analysis/CombineHarvester/blob/master/CombineTools/scripts/plot1DScan.py#L35
def BuildScan(graph, color, yvals):
    bestfit = None
    for i in xrange(graph.GetN()):
        #print('graph.GetY()[i]', graph.GetY()[i])
        if graph.GetY()[i] == 0.:
            bestfit = graph.GetX()[i]
    graph.SetMarkerColor(color)
    spline = ROOT.TSpline3("spline3", graph)
    NAMECOUNTER = 0
    func = ROOT.TF1('splinefn'+str(NAMECOUNTER), partial(Eval, spline), graph.GetX()[0], graph.GetX()[graph.GetN() - 1], 1)
    NAMECOUNTER += 1
    func.SetLineColor(color)
    func.SetLineWidth(3)
    assert(bestfit is not None)
    crossings = {}
    cross_1sig = None
    cross_2sig = None
    other_1sig = []
    other_2sig = []
    val = None
    val_2sig = None
    for yval in yvals:
        crossings[yval] = plot.FindCrossingsWithSpline(graph, func, yval)
        for cr in crossings[yval]:
            cr["contains_bf"] = cr["lo"] <= bestfit and cr["hi"] >= bestfit
    for cr in crossings[yvals[0]]:
        if cr['contains_bf']:
            val = (bestfit, cr['hi'] - bestfit, cr['lo'] - bestfit)
            cross_1sig = cr
        else:
            other_1sig.append(cr)
    if len(yvals) > 1:
        for cr in crossings[yvals[1]]:
            if cr['contains_bf']:
                val_2sig = (bestfit, cr['hi'] - bestfit, cr['lo'] - bestfit)
                cross_2sig = cr
            else:
                other_2sig.append(cr)
    else:
        val_2sig = (0., 0., 0.)
        cross_2sig = cross_1sig

    #print('crossings', crossings)
    #print('val', val)
    #print('val_2sig', val_2sig)
    #print('cross_1sig', cross_1sig)
    #print('cross_2sig', cross_2sig)
    #print('other_1sig', other_1sig)
    #print('other_2sig', other_2sig)

    #return func, crossings, val, val_2sig, cross_1sig, cross_2sig, other_1sig, other_2sig
    return {
            "graph"     : graph,
            "spline"    : spline,
            "func"      : func,
            "crossings" : crossings,
            "val"       : val,
            "val_2sig": val_2sig,
            "cross_1sig" : cross_1sig,
            "cross_2sig" : cross_2sig,
            "other_1sig" : other_1sig,
            "other_2sig" : other_2sig
        }


def Eval(obj, x, params):
    return obj.Eval(x[0])


  ####  ##### #   # #      # #    #  ####
 #        #    # #  #      # ##   # #    #
  ####    #     #   #      # # #  # #
      #   #     #   #      # #  # # #  ###
 #    #   #     #   #      # #   ## #    #
  ####    #     #   ###### # #    #  ####

def Load_Canvas_Style():

    gStyle.SetCanvasBorderMode(0)
    gStyle.SetCanvasColor(0)
    gStyle.SetCanvasDefH(600)
    gStyle.SetCanvasDefW(600)
    gStyle.SetCanvasDefX(0)
    gStyle.SetCanvasDefY(0)
    gStyle.SetPadBorderMode(0)
    gStyle.SetPadColor(0)
    gStyle.SetPadGridX(0)
    gStyle.SetPadGridY(0)
    gStyle.SetGridColor(0)
    gStyle.SetGridStyle(3)
    gStyle.SetGridWidth(1)
    gStyle.SetFrameBorderMode(0)
    gStyle.SetFrameBorderSize(1)
    gStyle.SetFrameFillColor(0)
    gStyle.SetFrameFillStyle(0)
    gStyle.SetFrameLineColor(1)
    gStyle.SetFrameLineStyle(1)
    gStyle.SetFrameLineWidth(1)
    gStyle.SetHistLineColor(1)
    gStyle.SetHistLineStyle(0)
    gStyle.SetHistLineWidth(1)
    gStyle.SetEndErrorSize(2)
    gStyle.SetOptFit(1011)
    gStyle.SetFitFormat("5.4g")
    gStyle.SetFuncColor(2)
    gStyle.SetFuncStyle(1)
    gStyle.SetFuncWidth(1)
    gStyle.SetOptDate(0)
    gStyle.SetOptFile(0)
    gStyle.SetOptStat(0)
    gStyle.SetStatColor(0)
    gStyle.SetStatFont(42)
    gStyle.SetStatFontSize(0.04)
    gStyle.SetStatTextColor(1)
    gStyle.SetStatFormat("6.4g")
    gStyle.SetStatBorderSize(1)
    gStyle.SetStatH(0.1)
    gStyle.SetStatW(0.15)
    gStyle.SetPadTopMargin(0.07)
    gStyle.SetPadBottomMargin(0.13)
    gStyle.SetPadLeftMargin(0.16)
    gStyle.SetPadRightMargin(0.03)
    # gStyle.SetOptTitle(0)
    gStyle.SetOptTitle(1)
    gStyle.SetTitleFont(42)
    gStyle.SetTitleColor(1)
    gStyle.SetTitleTextColor(1)
    gStyle.SetTitleFillColor(10)
    gStyle.SetTitleFontSize(0.05)
    gStyle.SetTitleColor(1, "XYZ")
    gStyle.SetTitleFont(42, "XYZ")
    gStyle.SetTitleSize(0.06, "XYZ")
    gStyle.SetTitleXOffset(0.9)
    gStyle.SetTitleYOffset(1.25)
    gStyle.SetLabelColor(1, "XYZ")
    gStyle.SetLabelFont(42, "XYZ")
    gStyle.SetLabelOffset(0.007, "XYZ")
    gStyle.SetLabelSize(0.05, "XYZ")
    gStyle.SetAxisColor(1, "XYZ")
    gStyle.SetStripDecimals(1)
    gStyle.SetTickLength(0.03, "XYZ")
    gStyle.SetNdivisions(510, "XYZ")
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetOptLogx(0)
    gStyle.SetOptLogy(0)
    gStyle.SetOptLogz(0)
    gStyle.SetPaperSize(20.,20.)


######## ######## ######## ########  ##        #######  ########
##       ##          ##    ##     ## ##       ##     ##    ##
##       ##          ##    ##     ## ##       ##     ##    ##
######   ######      ##    ########  ##       ##     ##    ##
##       ##          ##    ##        ##       ##     ##    ##
##       ##          ##    ##        ##       ##     ##    ##
######## ##          ##    ##        ########  #######     ##

class EFTPlot(object):

 # #    # # #####
 # ##   # #   #
 # # #  # #   #
 # #  # # #   #
 # #   ## #   #
 # #    # #   #

    def __init__(self, opts, mode='1D'):

        self.SM_mus = opts["SM_mus"]
        self.wcs = opts["wcs"]
        self.wcs_pairs = opts["wcs_pairs"]
        self.wc_ranges = opts["wc_ranges"]
        self.SMmu_ranges = opts["SMmu_ranges"]
        self.wc_ranges_scan1D = opts['wc_ranges_scan1D']
        self.wc_ranges_scan2D = opts['wc_ranges_scan2D']

        self.logger = logging.getLogger(__name__)
        self.ContourHelper = ContourHelper()
        self.histosFileName = 'Histos.root'

        #-- Define operator names for plotting
        # self.texdic = {'ctp': '#it{c}_{t#varphi}/#Lambda^{2}', 'ctG': '#it{c}_{tG}/#Lambda^{2}', 'cbW': '#it{c}_{bW}/#Lambda^{2}', 'cptb': '#it{c}_{#varphitb}/#Lambda^{2}', 'cQl3': '#it{c}^{3(#it{l})}_{Ql}/#Lambda^{2}', 'cQlM': '#it{c}^{-(#it{l})}_{Ql}/#Lambda^{2}', 'cQe': '#it{c}^{(#it{l})}_{Qe}/#Lambda^{2}', 'ctl': '#it{c}^{(#it{l})}_{tl}/#Lambda^{2}', 'cte': '#it{c}^{(#it{l})}_{te}/#Lambda^{2}', 'ctlS': '#it{c}^{S(#it{l})}_{t}/#Lambda^{2}', 'ctlT': '#it{c}^{T(#it{l})}_{t}/#Lambda^{2}'}
        self.texdic = {'ctw': '#it{c}_{tW}/#Lambda^{2}', 'ctz': '#it{c}_{tZ}/#Lambda^{2}', 'cpqm': '#it{c}^{-}_{#varphiQ}/#Lambda^{2}', 'cpq3': '#it{c}^{3(#it{l})}_{#varphiQ}/#Lambda^{2}', 'cpt': '#it{c}_{#varphit}/#Lambda^{2}'}
        self.texdicfrac = {'ctw': '#frac{#it{c}_{tW}}{#Lambda^{2}}', 'ctz': '#frac{#it{c}_{tZ}}{#Lambda^{2}}', 'cpqm': '#frac{#it{c}^{-}_{#varphiQ}}{#Lambda^{2}}', 'cpq3': '#frac{#it{c}^{3(#it{l})}_{#varphiQ}}{#Lambda^{2}}', 'cpt': '#frac{#it{c}_{#varphit}}{#Lambda^{2}}'}
        self.texdicrev = {v: k for k,v in self.texdic.items()}

        left = 0.17
        top = 0.91
        self.CMS_text = ROOT.TLatex(left, top, "CMS")
        self.CMS_text.SetNDC()
        self.CMS_text.SetTextColor(ROOT.kBlack)
        self.CMS_text.SetTextFont(61)
        self.CMS_text.SetTextAlign(11)
        self.CMS_text.SetTextSize(0.06)

        self.extraText = ROOT.TLatex(left+0.11, top, "Preliminary")
        self.extraText.SetNDC()
        self.extraText.SetTextFont(52)
        self.extraText.SetTextSize(0.05)

        self.lumiText = ROOT.TLatex(0.95, top, "137 fb^{-1} (13 TeV)")
        if mode == '2D': self.lumiText = ROOT.TLatex(0.82, top, "137 fb^{-1} (13 TeV)") #Need more space on right side for z-axis
        self.lumiText.SetNDC()
        self.lumiText.SetTextFont(42)
        self.lumiText.SetTextAlign(31)
        self.lumiText.SetTextSize(0.04)

        # Style
        # //--------------------------------------------
        self.linestyle = 7 #7=small dashes, 9=long dashes
        self.fillstyle = 1001 # 1001=solid,3001/2/3=dense,3004=strips

        # Logger
        # //--------------------------------------------
        log_file = 'plotter.log'

        FORMAT1 = '%(message)s'
        FORMAT2 = '[%(levelname)s] %(message)s'
        FORMAT3 = '[%(levelname)s][%(name)s] %(message)s'

        frmt1 = logging.Formatter(FORMAT1)
        frmt2 = logging.Formatter(FORMAT2)
        frmt3 = logging.Formatter(FORMAT3)

        logging.basicConfig(
            level=logging.DEBUG,
            format=FORMAT2,
            filename=log_file,
            filemode='w'
        )

        # Configure logging to also output to stdout
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(frmt2)
        logging.getLogger('').addHandler(console)

        return

    def ResetHistoFile(self, name=''):
        ROOT.TFile('Histos{}.root'.format(name),'RECREATE')
        self.histosFileName = 'Histos{}.root'.format(name)


 #     # #       #          ######                           #   ######
 ##    # #       #          #     # #       ####  #####     ##   #     #
 # #   # #       #          #     # #      #    #   #      # #   #     #
 #  #  # #       #          ######  #      #    #   #        #   #     #
 #   # # #       #          #       #      #    #   #        #   #     #
 #    ## #       #          #       #      #    #   #        #   #     #
 #     # ####### #######    #       ######  ####    #      ##### ######


    def Plot_NLLscan_1D(self, param='', log=False, paper=False, filepath='', name=''):
        '''
        Plot the NLL function versus a single POI (1D scan).
        NB: this is the prefered function, which uses CombineTool helper functions for convenience.

        wc: parameter to plot

        NB: 'NLL' = negative profiled log-likelihood function, as read in TFile.
        '''

        superimpose_second_scan = True #True <-> will look for second rootfile (hardcoded path) to overlay

        if filepath == '': filepath = './higgsCombine{}.MultiDimFit.mH120.root'.format(name)

        if not param:
            logging.error("No param specified!")
            return
        if not os.path.exists(filepath):
            logging.error("File " + filepath + " does not exist!".format(name))
            return

        if name == '': name = param

        logging.info(colors.fg.lightblue + "Enter function Plot_NLLscan_1D()" + colors.reset + '\n')
        print(colors.fg.orange + 'Reading file: ' + filepath + colors.reset)
        print(colors.fg.orange + 'Param: ' + param + colors.reset)

        ROOT.gROOT.SetBatch(True)
        c = ROOT.TCanvas('','',1000,800)
        #pads = plot.OnePad() #Helper func
        c.SetGrid(1)
        c.SetTopMargin(0.1)
        l = c.GetLeftMargin()

        colorLine_68 = 1 #68% band (12=gray, ...)
        colorLine_95 = 1 #95% band (ROOT.kAzure-8, ...)
        colorFill_68 = ROOT.kMagenta-1
        colorFill_95 = ROOT.kAzure-8

        yvals = [1., 3.84] #1sigma, 95%CL intervals

        #-- Get scan TTree
        rootFile = ROOT.TFile.Open(filepath)
        limitTree = rootFile.Get('limit')
        if limitTree.GetEntries()<2: print(colors.fg.red + 'ERROR: <2 entries in file {}, are you sure you performed a grid scan ?'.format(filepath) + colors.reset)

        #-- Get x-y coordinates for TGraph representing NLL function

        #-- Use CombineTool utils (see: https://github.com/cms-analysis/CombineHarvester/blob/master/CombineTools/python/plotting.py)
        graph = plot.TGraphFromTree(limitTree, param, '2*deltaNLL', 'quantileExpected > -1.5')
        #print(graph.GetN())
        graph.Sort()
        plot.RemoveGraphXDuplicates(graph)
        plot.RemoveGraphYAbove(graph, 7.)
        graph.SetLineColor(ROOT.kBlack)
        graph.SetLineWidth(4)

        #-- Overlay second scan (expected)
        if superimpose_second_scan == True:
            filepath = './higgsCombine.{}Exp.MultiDimFit.mH120.root'.format(param) #HARDCODED naming convention
            if not os.path.exists(filepath):
                logging.error("File " + filepath + " does not exist ! Can't overlay second TGraph...".format(name))
                superimpose_second_scan = False

        if superimpose_second_scan == True:
            rootFile = ROOT.TFile.Open(filepath)
            limitTree = rootFile.Get('limit')
            if limitTree.GetEntries()<2: print(colors.fg.red + 'ERROR: <2 entries in file {}, are you sure you performed a grid scan ?'.format(filepath) + colors.reset)
            graph2 = plot.TGraphFromTree(limitTree, param, '2*deltaNLL', 'quantileExpected > -1.5')
            graph2.Sort()
            plot.RemoveGraphXDuplicates(graph2)
            plot.RemoveGraphYAbove(graph2, 7.)
            graph2.SetLineColor(ROOT.kGray)
            graph2.SetLineWidth(3)

        #-- Xmin, xmax
        if SM:
            if superimpose_second_scan == False:
                xmin = max([graph.GetXaxis().GetXmin(), self.SMmu_ranges[param][0]], key=abs)
                xmax = max([graph.GetXaxis().GetXmax(), self.SMmu_ranges[param][1]], key=abs)
            else:
                xmin = max([graph.GetXaxis().GetXmin(), graph2.GetXaxis().GetXmin()], key=abs)
                xmax = max([graph.GetXaxis().GetXmax(), graph2.GetXaxis().GetXmax()], key=abs)
        else:
            if superimpose_second_scan == False:
                xmin = max([graph.GetXaxis().GetXmin(), self.wc_ranges_scan1D[param][0]], key=abs)
                xmax = max([graph.GetXaxis().GetXmax(), self.wc_ranges_scan1D[param][1]], key=abs)
            else:
                xmin = max([graph.GetXaxis().GetXmin(), graph2.GetXaxis().GetXmin()], key=abs)
                xmax = max([graph.GetXaxis().GetXmax(), graph2.GetXaxis().GetXmax()], key=abs)
        # print('xmin', xmin, 'xmax', xmax)

        #-- Create and draw filled TGraph (shaded area below NLL function representing 2sigmas area)
        #NB: moved here because needs to be drawn first, else hides other elements
        main_scan = BuildScan(graph, ROOT.kBlack, yvals)
        #main_scan['func'].Draw("same") #Can also draw the corresponding smooth function

        crossings = main_scan['crossings'][yvals[0]][0]
        fgraph68 = None
        if crossings['valid_lo'] and crossings['valid_hi']:
            fgraph68 = graph.Clone()
            ipt = 0
            while ipt<fgraph68.GetN():
                x, y = ROOT.Double(0), ROOT.Double(0) #Necessary to pass by reference in GetPoint()
                fgraph68.GetPoint(ipt, x, y)
                #print('x,y=', x, y)
                if x < crossings['lo'] or x > crossings['hi']:
                    # print(ipt, x, y)
                    fgraph68.RemovePoint(ipt)
                else: ipt+= 1
        # Trick: manually set first and last points to get the filled area properly defined
        if fgraph68 is not None: #Only draw filled area if there are multiple intersections (else, need to care about direction)
            fgraph68.GetPoint(0, x, y) #Trick: can't do InsertPointBefore(0) -> Save first point, modify it, insert it next
            fgraph68.SetPoint(0, crossings['lo'], 0) #Add first point at y=0
            fgraph68.InsertPointBefore(1, crossings['lo'], yvals[0]) #Add second point at y=Y
            fgraph68.InsertPointBefore(2, x, y) #cf. trick
            fgraph68.SetPoint(fgraph68.GetN(), crossings['hi'], yvals[0]) #Add first-to-last point at y=Y
            fgraph68.SetPoint(fgraph68.GetN(), crossings['hi'], 0) #Add last point at y=0
            fgraph68.SetFillColorAlpha(colorFill_68, 0.60) #Gray
            fgraph68.SetFillStyle(self.fillstyle)
            fgraph68.SetLineColor(ROOT.kBlack); fgraph68.SetLineWidth(0)

        #-- Idem for 95%CL graph
        crossings = main_scan['crossings'][yvals[1]][0]
        fgraph95 = None
        if crossings['valid_lo'] and crossings['valid_hi']:
            fgraph95 = graph.Clone()
            ipt = 0
            while ipt<fgraph95.GetN():
                x, y = ROOT.Double(0), ROOT.Double(0) #Necessary to pass by reference in GetPoint()
                fgraph95.GetPoint(ipt, x, y)
                #print('x,y=', x, y)
                if x < crossings['lo'] or x > crossings['hi']:
                    #print(ipt, x, y, 'REMOVED', crossings['lo'])
                    fgraph95.RemovePoint(ipt)
                else: ipt+= 1

            list_ipt_overlap68 = [] #List the indices of the points in fgraph95 which overlap with fgraph68, so that we can remove them in the plot
            for ipt in range(0, fgraph95.GetN()):
                x, y = ROOT.Double(0), ROOT.Double(0) #Necessary to pass by reference in GetPoint()
                fgraph95.GetPoint(ipt, x, y)
                # print('ipt', ipt, 'x', x, 'y', y)
                if x > main_scan['crossings'][yvals[0]][0]['lo'] and x < main_scan['crossings'][yvals[0]][0]['hi']:
                    fgraph95.SetPoint(ipt, x, 0) #Set overlap points to y=0
                    list_ipt_overlap68.append(ipt)
            #Add points to adjust the band perfectly on the crossing points #NB: the weird indices (+3...) are due to the fact that we keep adding point to the graph !
            fgraph95.InsertPointBefore(list_ipt_overlap68[0], main_scan['crossings'][yvals[0]][0]['lo'], 0)
            fgraph95.InsertPointBefore(list_ipt_overlap68[0], main_scan['crossings'][yvals[0]][0]['lo'], yvals[0])
            fgraph95.InsertPointBefore(list_ipt_overlap68[len(list_ipt_overlap68)-1]+3, main_scan['crossings'][yvals[0]][0]['hi'], yvals[0])
            fgraph95.InsertPointBefore(list_ipt_overlap68[len(list_ipt_overlap68)-1]+3, main_scan['crossings'][yvals[0]][0]['hi'], 0)

        # Trick: manually set first and last points to get the filled area properly defined
        if fgraph95 is not None: #Only draw filled area if there are multiple intersections (else, need to care about direction)
            fgraph95.GetPoint(0, x, y) #Trick: can't do InsertPointBefore(0) -> Save first point, modify it, insert it next
            fgraph95.SetPoint(0, crossings['lo'], 0) #Add first point at y=0
            fgraph95.InsertPointBefore(1, crossings['lo'], yvals[1]) #Add second point at y=Y
            fgraph95.InsertPointBefore(2, x, y) #cf. trick
            fgraph95.SetPoint(fgraph95.GetN(), crossings['hi'], yvals[1]) #Add first-to-last point at y=Y
            fgraph95.SetPoint(fgraph95.GetN(), crossings['hi'], 0) #Add last point at y=0
            # fgraph95.SetFillColorAlpha(ROOT.kAzure-7, 0.30) #Blue
            fgraph95.SetFillColorAlpha(colorFill_95, 0.50) #Blue
            fgraph95.SetFillStyle(self.fillstyle)
            fgraph95.SetLineColor(ROOT.kBlack); fgraph95.SetLineWidth(0)

            #-- Printouts
            #print(fgraph95.GetN())
            #for ipt in range(0, fgraph95.GetN()):
            #    x, y = ROOT.Double(0), ROOT.Double(0) #Necessary to pass by reference in GetPoint()
            #    fgraph95.GetPoint(ipt, x, y)
            #    print(ipt, x, y)

        #-- Axes ranges
        graph.SetMinimum(0.) #Display 0 label (else, line not displayed properly in .eps)
        # graph.SetMinimum(0.001) #Don't display 0 label
        graph.SetMaximum(10.) #Arbitrary

        #-- Draw graphs
        graph.Draw("AL") #Create axes #A:axes, P: markers, L:line
        graph.GetXaxis().SetLimits(xmin,xmax)
        fgraph95.Draw("F same")
        fgraph68.Draw("F same")
        graph.Draw("L same") #Redraw line on top

        line = ROOT.TLine()
        line.SetLineStyle(self.linestyle)
        line.SetLineWidth(2)

        for yval in yvals: #For each CL value
            if yval == 1: line.SetLineColor(colorLine_68)
            elif yval == 3.84: line.SetLineColor(colorLine_95)

            #Full horizontal line
            #plot.DrawHorizontalLine(pads[0], line, yval) #Helper func
            # line.DrawLine(graph.GetXaxis().GetXmin(), yval, graph.GetXaxis().GetXmax(), yval)

            #Split horizontal lines
            line.DrawLine(xmin, yval, main_scan['crossings'][yval][0]['lo'], yval)
            line.DrawLine(main_scan['crossings'][yval][0]['hi'], yval, xmax, yval)

            #Vertical lines
            for cr in main_scan['crossings'][yval]: #For each CL value
                if yval == 1: line.SetLineColor(colorLine_68)
                elif yval == 3.84: line.SetLineColor(colorLine_68)

                if cr['valid_lo']: line.DrawLine(cr['lo'], 0, cr['lo'], yval) #Low bound
                if cr['valid_hi']: line.DrawLine(cr['hi'], 0, cr['hi'], yval) #High bound

        if superimpose_second_scan == True: graph2.Draw("L same") #Draw second graph last ('overlay')

        #-- Legend
        leg = ROOT.TLegend(0.37,0.71,0.75,0.86)
        if superimpose_second_scan == True: leg = ROOT.TLegend(0.28,0.67,0.88,0.86) #More space
        leg.SetTextSize(0.04)
        leg.SetBorderSize(0) #Remove legend border

        if superimpose_second_scan == True:
            leg.AddEntry(graph, "Profiled log-likelihood (Observed)", "L")
            leg.AddEntry(graph2, "Profiled log-likelihood (Expected)", "L")
        else: leg.AddEntry(graph, "Profiled log-likelihood", "L")
        crossings = main_scan['crossings'][yvals[0]][0]
        #print('crossings', crossings)
        if crossings['valid_lo'] and crossings['valid_hi'] and fgraph68 is not None:
            leg.AddEntry(fgraph68, "68% CL [{:.2f}, {:.2f}]".format(crossings['lo'],crossings['hi']), "F")
        crossings = main_scan['crossings'][yvals[1]][0]
        if crossings['valid_lo'] and crossings['valid_hi'] and fgraph95 is not None:
            leg.AddEntry(fgraph95, "95% CL [{:.2f}, {:.2f}]".format(crossings['lo'],crossings['hi']), "F")
            print(colors.fg.orange + "95% CL interval: [" + str(crossings['lo']) + ", " + str(crossings['hi']) + "]" + colors.reset)

        #-- Look for additional other (disconnected) intervals
        if(len(main_scan['crossings'][yvals[1]])>1):
            crossings = main_scan['crossings'][yvals[1]][1]
            if crossings['valid_lo'] and crossings['valid_hi'] and fgraph95 is not None:
                leg.AddEntry(fgraph95, "95% CL [{:.2f}, {:.2f}]".format(crossings['lo'],crossings['hi']), "F")
                print(colors.fg.orange + "95% CL interval: [" + str(crossings['lo']) + ", " + str(crossings['hi']) + "]" + colors.reset)
        if(len(main_scan['crossings'][yvals[1]])>2):
            crossings = main_scan['crossings'][yvals[1]][2]
            if crossings['valid_lo'] and crossings['valid_hi'] and fgraph95 is not None:
                leg.AddEntry(fgraph95, "95% CL [{:.2f}, {:.2f}]".format(crossings['lo'],crossings['hi']), "F")
                print(colors.fg.orange + "95% CL interval: [" + str(crossings['lo']) + ", " + str(crossings['hi']) + "]" + colors.reset)

        leg.Draw("same")

        #-- Labels
        graph.SetTitle("")
        # graph.SetMarkerStyle(26) #Change markers from invisible dots to nice triangles
        graph.SetMarkerStyle(8)
        graph.SetMarkerSize(1)
        # graph.GetXaxis().SetTitle(param)
        graph.GetXaxis().SetTitle(Get_Parameter_LegName(param))
        graph.GetYaxis().SetTitle("-2 #Delta log(L)")
        # graph.GetYaxis().SetTitle("{} 2#DeltaNLL".format(param))

        self.CMS_text.Draw('same')
        if paper==False: self.extraText.Draw('same')
        self.lumiText.Draw('same')

        #-- Save
        outname = 'scan1D{}'.format(name.replace('.','_'))
        if log:
            graph.SetMinimum(0.1) #Crash if zero
            graph.SetLogz()
            outname+= '_log'

        if paper: outname+= '_paper.eps' #.pdf
        else: outname+= '.png'

        c.Print(outname)

        return


 #    #   ##   #    # #    #   ##   #
 ##  ##  #  #  ##   # #    #  #  #  #
 # ## # #    # # #  # #    # #    # #
 #    # ###### #  # # #    # ###### #
 #    # #    # #   ## #    # #    # #
 #    # #    # #    #  ####  #    # ######

    def Plot1DManualNLLScan(self, param, paper=False, name=''):
        '''
        Plot a NLL scan from multiple fixed points (NLL scanned manually with algo=fixed rather than automatically with algo=grid)
        '''

        if not param:
            logging.error("No param specified!")
            return

        logging.info(colors.fg.lightblue + "\nEnter function Plot1DManualNLLScan" + colors.reset)

        if name == '': name = param

        ROOT.gROOT.SetBatch(True)
        c = ROOT.TCanvas('','',1000,800)
        #c = ROOT.TCanvas('','',2000,1600)
        c.SetGrid(1)
        c.SetTopMargin(0.1)
        l = c.GetLeftMargin()

        #WC_values = [-4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 3, 4]
        #WC_values = [-4, -2, -1, 0, 1, 2, 4]
        #WC_values = [-2, -1.5, -1, 0, 1, 2]
    	#WC_values = [-1.5, -1.2, -1, -0.8, -0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5, 0.8, 1, 1.2, 1.5]
    	WC_values = [-1.5, -1, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.5]

        list_filepaths = []
        for val in WC_values:
            filepath = './higgsCombine_'+param+'_'+str(val)+'.MultiDimFit.mH120.root'
            if not os.path.exists(filepath):
                logging.error("File " + filepath + " does not exist!")
                return
            else:
                list_filepaths.append(filepath)

        #-- Get x-y coordinates for TGraph representing NLL function
        NLL_values = []
        for filepath in list_filepaths:
            print('filepath', filepath)

            #-- Get scan TTree
            rootFile = ROOT.TFile.Open(filepath)

            limitTree = rootFile.Get('limit')
            nentries = limitTree.GetEntries()
            limitTree.GetEntry(nentries-1) #This is where we read the relevant entry #Not sure why, there is sometimes 1 or 2 entries

            print('deltaNLL', limitTree.GetLeaf('deltaNLL').GetValue(0)) #For non-arrays, always read leaf element 0
            print('nll0', limitTree.GetLeaf('nll0').GetValue(0))
            print('nll', limitTree.GetLeaf('nll').GetValue(0))

            #NLL_tmp = limitTree.GetLeaf('deltaNLL').GetValue(0) + limitTree.GetLeaf('nll0').GetValue(0) + limitTree.GetLeaf('nll').GetValue(0) #nll+nll0+deltaNLL
            NLL_tmp = limitTree.GetLeaf('deltaNLL').GetValue(0) #only deltaNLL
            NLL_values.append(2 * NLL_tmp)

        graph = ROOT.TGraph(len(WC_values), numpy.asarray(WC_values, dtype=np.double), numpy.asarray(NLL_values))

        for ipt in range(graph.GetN()):
            x, y = ROOT.Double(0), ROOT.Double(0) #Necessary to pass by reference in GetPoint()
            graph.GetPoint(ipt, x, y)
            print(ipt, x, y)

        f = ROOT.TF1("f", "[2] * x * x + [1] * x + [0]") #Quadratic polynomial
        # f = ROOT.TF1("f", "[4] * x * x * x * x + [3] * x * x * x + [2] * x * x + [1] * x + [0]")

        min = TMath.MinElement(graph.GetN(),graph.GetY());
        #print(min)

        NLL_values = [val-min for val in NLL_values] #Substract ymin
        #print(NLL_values)
        graph = ROOT.TGraph(len(WC_values), numpy.asarray(WC_values, dtype=np.double), numpy.asarray(NLL_values)) #Update graph

        #graph.Fit(f)

        # for ipt in range(graph.GetN()):
        #     x, y = ROOT.Double(0), ROOT.Double(0) #Necessary to pass by reference in GetPoint()
        #     graph.GetPoint(ipt, x, y)
        #     print(ipt, x, y)

        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(1.3)
        graph.SetMarkerColor(ROOT.kBlue)
        graph.Draw("APL") #A:axes, P: markers, L:line

        #-- Axis ranges
        xmin = graph.GetXaxis().GetXmin()
        xmax = graph.GetXaxis().GetXmax()
        graph.SetMinimum(0.001) #Don't display 0 label
        graph.SetMaximum(10.) #Arbitrary
        xmin = self.wc_ranges[param][0]
        xmax = self.wc_ranges[param][1]
        graph.GetXaxis().SetRangeUser(xmin,xmax)

        yval_68, list_X_intersects68, fgraph68, lines_X_intersects68, yval_95, list_X_intersects95, fgraph95, lines_X_intersects95 = Get_nSigmaBand_Graphs(graph)
        fgraph68.Draw("F same")
        fgraph95.Draw("F same")

        #-- Axis ranges
        #xmin = graph.GetXaxis().GetXmin()
        #xmax = graph.GetXaxis().GetXmax()
        #graph.SetMinimum(0.001) #Don't display 0 label
        graph.SetMaximum(10.) #Arbitrary

        #-- Legend
        leg = ROOT.TLegend(0.37,0.70,0.75,0.86)
        leg.SetTextSize(0.04)
        leg.AddEntry(graph, "Profiled log-likelihood", "L")
        if len(list_X_intersects68)>=2: leg.AddEntry(fgraph68, "68% CL [{:.2f}, {:.2f}]".format(list_X_intersects68[0],list_X_intersects68[1]), "F")
        if len(list_X_intersects95)>=2: leg.AddEntry(fgraph95, "95% CL [{:.2f}, {:.2f}]".format(list_X_intersects95[0],list_X_intersects95[1]), "F")
        leg.Draw("same")
        if len(list_X_intersects95)>1: print(colors.fg.orange + "95% interval: [" + str(list_X_intersects95[0]) + ", " + str(list_X_intersects95[1]) + "]" + colors.reset)

        #-- Labels
        graph.SetTitle("")
        graph.SetMarkerStyle(8)
        graph.SetMarkerSize(1)
        graph.GetXaxis().SetTitle(Get_Parameter_LegName(param))
        graph.GetYaxis().SetTitle("-2 #Delta log(L)")

        cmsText = "CMS";
        latex = ROOT.TLatex()
        latex.SetNDC();
        latex.SetTextColor(ROOT.kBlack);
        latex.SetTextFont(61);
        latex.SetTextAlign(11);
        latex.SetTextSize(0.06);
        latex.DrawLatex(l+0.01, 0.92, cmsText)

        extraText = "Preliminary simulation";
        latex.SetTextFont(52);
        latex.SetTextSize(0.04);
        if paper==False: latex.DrawLatex(l+0.12, 0.92, extraText)

        latex.SetTextFont(42);
        latex.SetTextAlign(31);
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.96, 0.92, "41.5 fb^{-1} (13 TeV)");

        outname = 'scan1D_manual{}'.format(name.replace('.','_'))
        if paper: outname+= '_paper.eps' #.pdf
        else: outname+= '.png'

        c.Print(outname)

        return


 #     # #       #          ######                          #####  ######
 ##    # #       #          #     # #       ####  #####    #     # #     #
 # #   # #       #          #     # #      #    #   #            # #     #
 #  #  # #       #          ######  #      #    #   #       #####  #     #
 #   # # #       #          #       #      #    #   #      #       #     #
 #    ## #       #          #       #      #    #   #      #       #     #
 #     # ####### #######    #       ######  ####    #      ####### ######

    def Plot_NLLscan_2D(self, mode='EFT', params=[], ceiling=1, log=False, paper=False, filepath='', name=''):
        '''
        Plot the NLL function versus 2 POIs (2D scan).

        params: POIs (signal strengths of 2 processes of interest)
        mode: 'SM' or 'EFT'
        ceiling: maximum NLL value
        '''

        ROOT.gStyle.SetPalette(57) #kBird (default, blue to yellow)

        # ROOT.gStyle.SetPalette(ROOT.kDeepSea) #Dark blues
        # ROOT.gStyle.SetPalette(ROOT.kOcean) #Dark greens/blues
        # ROOT.gStyle.SetPalette(ROOT.kDarkBodyRadiator) #Dark to light
        # ROOT.gStyle.SetPalette(ROOT.kBlackBody) #Light red to blue
        # ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)
        # ROOT.gStyle.SetPalette(ROOT.kLightTemperature)

        drawContours = True #True <-> superimpose contours (default)
        valCont = [2.30, 5.99] #1sigma and 95% contours

        if filepath == '': filepath = './higgsCombine.'+mode+'.MultiDimFit.mH120.root'

        if len(params) != 2:
            logging.error("Function 'Plot_NLLscan_2D' requires exactly two parameters!")
            return
        if not os.path.exists(filepath):
            logging.error("File {} does not exist!".format(filepath))
            return

        logging.info(colors.fg.lightblue + "Enter function Plot_NLLscan_2D()\n" + colors.reset)

        if name == '': name = param

        ROOT.gROOT.SetBatch(True)
        c = ROOT.TCanvas('','',1000,800)
        c.SetGrid(1)
        c.SetTopMargin(0.1)
        c.SetRightMargin(0.17);
        l = c.GetLeftMargin()

        c.SetTopMargin(0.15) #Leave space for legend

        hname = 'scan2D'
        hname+= name.replace('.','_') if name not in ['','EFT','SM'] else params[0]+'_'+params[1]
        if log: hname += "_log"

        #-- Open file and draw 2D histogram
        xvar = params[0]; yvar = params[1]
        logging.info(colors.fg.orange + "Open file : " + colors.reset)
        logging.info(colors.fg.orange + "xvar={} / yvar={}".format(xvar,yvar) + colors.reset)
        logging.info(filepath)
        rootFile = ROOT.TFile.Open(filepath)
        limitTree = rootFile.Get('limit')
        minZ = limitTree.GetMinimum('deltaNLL')
        xmin = limitTree.GetMinimum(xvar)
        xmax = limitTree.GetMaximum(xvar)
        ymin = limitTree.GetMinimum(yvar)
        ymax = limitTree.GetMaximum(yvar)
        # xmin = -40
        # xmax = 40
        # ymin = -40
        # ymax = 40

        # maxZ = 20 #Max z-axis threshold
        maxZ = 1000 #Max z-axis threshold #If want whole plot to be colored
        nbins = 150

    	#-- Auto-binning
    	nbins = (int) (math.sqrt(limitTree.GetEntries())) #If there are N^2 points, will use N bins in each axis
    	# print('nbins', nbins)

        #-- WARNING: zvar:yvar:xvar (cf. https://root.cern.ch/doc/master/classTTree.html... although it says different?)
        expr = '2*(deltaNLL-{}):{}:{}>>{}({},{},{},{},{},{})'.format(minZ, yvar, xvar, hname, nbins, xmin, xmax, nbins, ymin, ymax)
        limitTree.Draw(expr, '2*deltaNLL<{}'.format(maxZ), 'prof colz') #Option 'prof' generates a TProfile

        hist = c.GetPrimitive(hname)
        # hist.Draw('colz') #Clean draw (not necessary)

        #-- Trick: can happen that just few points fail, an corresponding bins set 0 contents
        #--> set bin content to same value as neighbouring bin (usually safe because happen far away from minimum)
        #NB: for some reason, can't use SetBinContent on 'hist' object directly; so duplicate it in new, modificable TH2F object
        newhist = ROOT.TH2F('', '', nbins, xmin, xmax, nbins, ymin, ymax)
        for xbin in range(1,hist.GetNbinsX()+1):
            # print('-- xbin', xbin)
            for ybin in range(1,hist.GetNbinsY()+1):
                # print('- ybin', ybin)
                newhist.SetBinContent(xbin, ybin, hist.GetBinContent(xbin,ybin))
                if newhist.GetBinContent(xbin,ybin)==0:
                    if ybin>1:
                        # print('newhist.GetBinContent(xbin,ybin)', newhist.GetBinContent(xbin,ybin))
                        newhist.SetBinContent(xbin,ybin, newhist.GetBinContent(xbin,ybin-1))
                        # print('--> newhist.GetBinContent(xbin,ybin)', newhist.GetBinContent(xbin,ybin))
                    elif xbin>1:
                        newhist.SetBinContent(xbin,ybin, newhist.GetBinContent(xbin-1,ybin))

        hist = newhist #Reuse default name
        hist.Draw('colz') #Clean draw

        if drawContours:
            histsForCont = hist.Clone() #Clone histo
            c_contlist = ((ctypes.c_double)*(len(valCont)))(*valCont) #Define contour thresholds ('levels')
            histsForCont.SetContour(len(c_contlist), c_contlist) #Set the contours
            histsForCont.Draw("cont z list") #Get the contour in memory
            c.Update()
            conts = ROOT.gROOT.GetListOfSpecials().FindObject("contours") #Retrieve the object
            cont_p1 = conts.At(0).Clone()
            cont_p2 = conts.At(1).Clone()

            hist.Draw('colz') #Fresh draw to overwrite previous one

            for conts in [cont_p1]:
                for cont in conts:
                    #cont.SetLineColor(ROOT.kBlue-6)
                    cont.SetLineColor(ROOT.kBlack)
                    cont.SetLineWidth(3)
                    cont.SetLineStyle(self.linestyle)
                    cont.Draw("L same") #Draw second contour
                    h_cont68 = cont.Clone()
            for conts in [cont_p2]:
                for cont in conts:
                    cont.SetLineColor(ROOT.kBlack)
                    cont.SetLineWidth(3)
                    cont.Draw("L same") #Draw first contour
                    h_cont95 = cont.Clone()

        #-- Marker for SM point
        marker_SM = ROOT.TMarker()
        marker_SM.SetMarkerSize(3.)
        marker_SM.SetMarkerStyle(34)
        marker_SM.SetMarkerColor(ROOT.kBlack)
        marker_SM.DrawMarker(0,0) #SM=(0,0)

        #-- Draw best fit point from grid scan (for observed fit only)
        if 'exp' not in name and 'Exp' not in name:
            for ievt in range(limitTree.GetEntries()):
                limitTree.GetEntry(ievt)
                if limitTree.GetLeaf("quantileExpected").GetValue(ievt)==-1: #Best fit point has quantileExpected==-1
                    bestfit_x = limitTree.GetLeaf(params[0]).GetValue(ievt) #Read x
                    bestfit_y = limitTree.GetLeaf(params[1]).GetValue(ievt) #Read y
                    # print(bestfit_x, bestfit_y)
                    break #The minimum may be included multiple times in the file (if stitched scans together) --> only read once
            best_fit = ROOT.TMarker()
            best_fit.SetMarkerSize(2)
            best_fit.SetMarkerStyle(47) #34 full +, 47 full cross, 2 slim +, 5 slim cross, 28 empty +, 46 empty cross, 4 empty circle, 8 full circle
            best_fit.SetMarkerColor(ROOT.kRed+1)
            if bestfit_x != 0. or bestfit_y != 0.: best_fit.DrawMarker(bestfit_x, bestfit_y)

        if log: c.SetLogz()
        hist.GetXaxis().SetTitle(Get_Parameter_LegName(params[0]))
        hist.GetYaxis().SetTitle(Get_Parameter_LegName(params[1]))
        hist.GetZaxis().SetTitle("-2 #Delta log(L)")
        hist.GetYaxis().SetTitleOffset(0.9)
        hist.SetMaximum(20.) #If want to cut z-axis below

        hist.SetTitle('')
        hist.SetStats(0)

        if paper:
            left = 0.17
            top = 0.91-0.04
            self.CMS_text = ROOT.TLatex(left, top, "CMS")
            self.CMS_text.SetNDC()
            self.CMS_text.SetTextColor(ROOT.kBlack)
            self.CMS_text.SetTextFont(61)
            self.CMS_text.SetTextAlign(11)
            self.CMS_text.SetTextSize(0.06)

            self.extraText = ROOT.TLatex(left+0.11, top, "Preliminary")
            self.extraText.SetNDC()
            self.extraText.SetTextFont(52)
            self.extraText.SetTextSize(0.05)

            self.lumiText = ROOT.TLatex(left+0.12, top, "137 fb^{-1} (13 TeV)")
            # self.lumiText = ROOT.TLatex(0.95, top, "137 fb^{-1} (13 TeV)")
            if mode == '2D': self.lumiText = ROOT.TLatex(0.82, top, "137 fb^{-1} (13 TeV)") #Need more space on right side for z-axis
            self.lumiText.SetNDC()
            self.lumiText.SetTextFont(42)
            self.lumiText.SetTextAlign(11)
            # self.lumiText.SetTextAlign(31)
            self.lumiText.SetTextSize(0.04)

        else: #Preliminary
            left = 0.17
            top = 0.91+0.02
            self.CMS_text = ROOT.TLatex(left, top, "CMS")
            self.CMS_text.SetNDC()
            self.CMS_text.SetTextColor(ROOT.kBlack)
            self.CMS_text.SetTextFont(61)
            self.CMS_text.SetTextAlign(11)
            self.CMS_text.SetTextSize(0.06)

            self.extraText = ROOT.TLatex(left+0.11, top, "Preliminary")
            self.extraText.SetNDC()
            self.extraText.SetTextFont(52)
            self.extraText.SetTextSize(0.05)

            self.lumiText = ROOT.TLatex(left+0.01, top-0.06, "137 fb^{-1} (13 TeV)")
            # self.lumiText = ROOT.TLatex(0.95, top, "137 fb^{-1} (13 TeV)")
            if mode == '2D': self.lumiText = ROOT.TLatex(0.82, top, "137 fb^{-1} (13 TeV)") #Need more space on right side for z-axis
            self.lumiText.SetNDC()
            self.lumiText.SetTextFont(42)
            self.lumiText.SetTextAlign(11)
            # self.lumiText.SetTextAlign(31)
            self.lumiText.SetTextSize(0.04)

        #-- For paper

        self.CMS_text.Draw('same')
        if paper==False: self.extraText.Draw('same')
        self.lumiText.Draw('same')

        #-- Legend
        # legend = ROOT.TLegend(0.53,0.70,0.80,0.87)
        legend = ROOT.TLegend(0.60,0.86,0.83,0.99)
        legend.AddEntry(h_cont68, "68% CL",'L')
        legend.AddEntry(h_cont95, "95% CL",'L')
        legend.AddEntry(marker_SM, "SM",'p')
        legend.AddEntry(best_fit, "Best fit",'p')
        legend.SetBorderSize(0)
        legend.SetTextSize(0.04)
        legend.SetNColumns(1)
        legend.Draw('same')

        #-- Save plot
        outname = hname
        if paper: outname+= '_paper.eps' #.pdf
        else: outname+= '.png'

        c.Print(outname)

        #-- Save to root file
        if not log:
            outfile = ROOT.TFile(self.histosFileName,'UPDATE')
            hist.Write()
            outfile.Close()


  ####  #    # ###### #####  #        ##   #   #    #####  #       ####  #####  ####
 #    # #    # #      #    # #       #  #   # #     #    # #      #    #   #   #
 #    # #    # #####  #    # #      #    #   #      #    # #      #    #   #    ####
 #    # #    # #      #####  #      ######   #      #####  #      #    #   #        #
 #    #  #  #  #      #   #  #      #    #   #      #      #      #    #   #   #    #
  ####    ##   ###### #    # ###### #    #   #      #      ######  ####    #    ####

    #-- Not used yet
    def OverlayLLPlot1DEFT(self, name1='.test', name2='.test', wc='', log=False):
        if not wc:
            logging.error("No wc specified!")
            return
        if not os.path.exists('./higgsCombine{}.MultiDimFit.root'.format(name1)):
            logging.error("File higgsCombine{}.MultiDimFit.root does not exist!".format(name1))
            return
        if not os.path.exists('./higgsCombine{}.MultiDimFit.root'.format(name2)):
            logging.error("File higgsCombine{}.MultiDimFit.root does not exist!".format(name2))
            return

        ROOT.gROOT.SetBatch(True)

        canvas = ROOT.TCanvas('canvas', 'canvas', 700, 530)
        p1 = ROOT.TPad('p1', 'p1', 0, 0.05, 1.0, 1.0)
        p1.Draw()
        p1.cd()

        # Get scan trees
        rootFile1 = ROOT.TFile.Open('./higgsCombine{}.MultiDimFit.root'.format(name1))
        limitTree1 = rootFile1.Get('limit')

        rootFile2 = ROOT.TFile.Open('./higgsCombine{}.MultiDimFit.root'.format(name2))
        limitTree2 = rootFile2.Get('limit')

        # Get coordinates for TGraphs
        graph1wcs = []
        graph2wcs = []
        graph1nlls = []
        graph2nlls = []
        for entry in range(limitTree1.GetEntries()):
            limitTree1.GetEntry(entry)
            graph1wcs.append(limitTree1.GetLeaf(wc).GetValue(0))
            graph1nlls.append(2*limitTree1.GetLeaf('deltaNLL').GetValue(0))
        for entry in range(limitTree2.GetEntries()):
            limitTree2.GetEntry(entry)
            graph2wcs.append(limitTree2.GetLeaf(wc).GetValue(0))
            graph2nlls.append(2*limitTree2.GetLeaf('deltaNLL').GetValue(0))

        # Rezero the y axis and make the tgraphs
        graph1nlls = [val-min(graph1nlls) for val in graph1nlls]
        graph2nlls = [val-min(graph2nlls) for val in graph2nlls]
        graph1 = ROOT.TGraph(len(graph1wcs),numpy.asarray(graph1wcs),numpy.asarray(graph1nlls))
        graph2 = ROOT.TGraph(len(graph2wcs),numpy.asarray(graph2wcs),numpy.asarray(graph2nlls))
        del graph1nlls,graph2nlls,graph1wcs,graph2wcs

        # Combine into TMultiGraph
        multigraph = ROOT.TMultiGraph()
        multigraph.Add(graph1)
        multigraph.Add(graph2)
        multigraph.Draw("AP")
        multigraph.GetXaxis().SetLabelSize(0.05)
        multigraph.GetYaxis().SetLabelSize(0.05)
        multigraph.GetXaxis().SetTitleSize(0.05)
        multigraph.GetXaxis().SetTitleOffset(0.8)
        #multigraph.GetXaxis().SetNdivisions(7)

        # Squeeze X down to whatever range captures the float points
        xmin = self.wc_ranges[wc][1]
        xmax = self.wc_ranges[wc][0]
        for idx in range(graph1.GetN()):
            if graph1.GetY()[idx] < 10 and graph1.GetX()[idx] < xmin:
                xmin = graph1.GetX()[idx]
            if graph1.GetY()[idx] < 10 and graph1.GetX()[idx] > xmax:
                xmax = graph1.GetX()[idx]
        multigraph.GetXaxis().SetRangeUser(xmin,xmax)
        multigraph.GetYaxis().SetRangeUser(-0.1,10)

        #Change markers from invisible dots to nice triangles
        graph1.SetMarkerColor(1)
        graph1.SetMarkerStyle(26)
        graph1.SetMarkerSize(1)

        graph2.SetMarkerColor(2)
        graph2.SetMarkerStyle(26)
        graph2.SetMarkerSize(1)

        #Add 1-sigma and 2-sigma lines. (Vertical lines were too hard, sadly)
        canvas.SetGrid(1)
        p1.SetGrid(1)

        line68 = ROOT.TLine(xmin,1,xmax,1)
        line68.Draw('same')
        line68.SetLineColor(ROOT.kYellow+1)
        line68.SetLineWidth(3)
        line68.SetLineStyle(self.linestyle)

        line95 = ROOT.TLine(xmin,4,xmax,4)
        line95.Draw('same')
        line95.SetLineColor(ROOT.kCyan-2)
        line95.SetLineWidth(3)
        line95.SetLineStyle(self.linestyle)

        # Labels
        Title = ROOT.TLatex(0.5, 0.95, "{} 2#DeltaNLL".format(wc))
        Title.SetNDC(1)
        Title.SetTextAlign(20)
        #Title.Draw('same')
        #multigraph.GetXaxis().SetTitle(wc)
        #multigraph.GetXaxis().SetTitle(self.texdic[wc.rstrip('i')])
        XTitle = ROOT.TLatex(0.85, 0.01, self.texdic[wc.rstrip('i')])
        XTitle.SetNDC(1)
        XTitle.SetTextAlign(20)
        XTitle.SetTextFont(42)
        canvas.cd()
        XTitle.Draw('same')

        # CMS-required text
        self.CMS_text.Draw('same')
        self.lumiText.Draw('same')

        # Legend
        legend = ROOT.TLegend(0.1,0.85,0.45,0.945)
        legend.AddEntry(graph1,"Others Profiled (2#sigma)",'p')
        legend.AddEntry(graph2,"Others Fixed to SM (2#sigma)",'p')
        legend.SetTextSize(0.035)
        legend.SetNColumns(1)
        legend.Draw('same')

        #Check log option, then save as image
        if log:
            multigraph.SetMinimum(0.1)
            multigraph.SetLogz()
            canvas.Print('Overlay{}1DNLL_log.png'.format(wc),'png')
        else:
            canvas.Print('Overlay{}1DNLL.png'.format(wc),'png')

        rootFile1.Close()
        rootFile2.Close()

        return

    def OverlayZoomLLPlot1DEFT(self, name1='.test', name2='.test', wc='', log=False):
        if not wc:
            logging.error("No wc specified!")
            return
        if not os.path.exists('./higgsCombine{}.MultiDimFit.root'.format(name1)):
            logging.error("File higgsCombine{}.MultiDimFit.root does not exist!".format(name1))
            return
        if not os.path.exists('./higgsCombine{}.MultiDimFit.root'.format(name2)):
            logging.error("File higgsCombine{}.MultiDimFit.root does not exist!".format(name2))
            return

        ROOT.gROOT.SetBatch(True)

        canvas = ROOT.TCanvas()

        # Get scan trees
        rootFile1 = ROOT.TFile.Open('./higgsCombine{}.MultiDimFit.root'.format(name1))
        limitTree1 = rootFile1.Get('limit')

        rootFile2 = ROOT.TFile.Open('./higgsCombine{}.MultiDimFit.root'.format(name2))
        limitTree2 = rootFile2.Get('limit')

        # Get coordinates for TGraphs
        graph1wcs = []
        graph2wcs = []
        graph1nlls = []
        graph2nlls = []
        for entry in range(limitTree1.GetEntries()):
            limitTree1.GetEntry(entry)
            graph1wcs.append(limitTree1.GetLeaf(wc).GetValue(0))
            graph1nlls.append(2*limitTree1.GetLeaf('deltaNLL').GetValue(0))
        for entry in range(limitTree2.GetEntries()):
            limitTree2.GetEntry(entry)
            graph2wcs.append(limitTree2.GetLeaf(wc).GetValue(0))
            graph2nlls.append(2*limitTree2.GetLeaf('deltaNLL').GetValue(0))

        # Rezero the y axis and make the tgraphs
        graph1nlls = [val-min(graph1nlls) for val in graph1nlls]
        graph2nlls = [val-min(graph2nlls) for val in graph2nlls]
        graph1 = ROOT.TGraph(len(graph1wcs),numpy.asarray(graph1wcs),numpy.asarray(graph1nlls))
        graph2 = ROOT.TGraph(len(graph2wcs),numpy.asarray(graph2wcs),numpy.asarray(graph2nlls))
        del graph1nlls,graph2nlls,graph1wcs,graph2wcs

        # Combine into TMultiGraph
        multigraph = ROOT.TMultiGraph()
        multigraph.Add(graph1)
        multigraph.Add(graph2)
        multigraph.Draw("AP")

        # Squeeze X down to 20 pts around 0.
        width = self.wc_ranges[wc][1]-self.wc_ranges[wc][0]
        xmin = -float(width)/50
        xmax = float(width)/50
        ymax = max(graph1.Eval(xmin),graph1.Eval(xmax),graph2.Eval(xmin),graph2.Eval(xmax))
        ymin = -ymax/10
        multigraph.GetXaxis().SetRangeUser(xmin, xmax)
        multigraph.GetYaxis().SetRangeUser(ymin, ymax)

        #Change markers from invisible dots to nice triangles
        graph1.SetMarkerColor(1)
        graph1.SetMarkerStyle(26)
        graph1.SetMarkerSize(1)

        graph2.SetMarkerColor(2)
        graph2.SetMarkerStyle(26)
        graph2.SetMarkerSize(1)

        #Add 1-sigma and 2-sigma lines. (Vertical lines were too hard, sadly)
        canvas.SetGrid(1)

        line68 = ROOT.TLine(xmin,1,xmax,1)
        line68.Draw('same')
        line68.SetLineColor(ROOT.kYellow+1)
        line68.SetLineWidth(3)
        line68.SetLineStyle(self.linestyle)

        line95 = ROOT.TLine(xmin,4,xmax,4)
        line95.Draw('same')
        line95.SetLineColor(ROOT.kCyan-2)
        line95.SetLineWidth(3)
        line95.SetLineStyle(self.linestyle)

        # Labels
        Title = ROOT.TLatex(0.5, 0.95, "{} 2#DeltaNLL".format(wc))
        Title.SetNDC(1)
        Title.SetTextAlign(20)
        Title.Draw('same')
        multigraph.GetXaxis().SetTitle(wc)

        # CMS-required text
        self.CMS_text.Draw('same')
        self.lumiText.Draw('same')

        #Check log option, then save as image
        if log:
            multigraph.SetMinimum(0.1)
            multigraph.SetLogz()
            canvas.Print('OverlayZoom{}1DNLL_log.png'.format(wc),'png')
        else:
            canvas.Print('OverlayZoom{}1DNLL.png'.format(wc),'png')

        rootFile1.Close()
        rootFile2.Close()


  ####   ####  #    # #####  ####  #    # #####   ####
 #    # #    # ##   #   #   #    # #    # #    # #
 #      #    # # #  #   #   #    # #    # #    #  ####
 #      #    # #  # #   #   #    # #    # #####       #
 #    # #    # #   ##   #   #    # #    # #   #  #    #
  ####   ####  #    #   #    ####   ####  #    #  ####

    def ContourPlotEFT(self, mode='EFT', params=[], paper=False, filepath='', name=''):
        '''
        Make 2D contour plots.
        '''

        if len(params)!=2:
            logging.error("Function 'ContourPlot' requires exactly two params!")
            return

        if filepath == '': filepath = './higgsCombine{}.MultiDimFit.mH120.root'.format(mode)

        if not os.path.exists(filepath):
            logging.error("File {} does not exist!".format(filepath))
            return

        ranges = self.wc_ranges #Default: EFT
        if mode=='SM': ranges = self.SMmu_ranges #SM

        if name == '': name = params[0] + '_' + params[1]

        #-- Get Grid scan and copy to h_contour
        #-- params[0] is x-axis variable, params[1] is y-axis variable
        gridFile = ROOT.TFile.Open(filepath)
        tree = gridFile.Get('limit')
        minZ = tree.GetMinimum('deltaNLL')

        ROOT.gROOT.SetBatch(True)
        c = ROOT.TCanvas('c','c',1000,800)
        c.SetGrid(1)
        c.SetTopMargin(0.1)
        l = c.GetLeftMargin()

        nbins = 150
        xvar = params[0]; yvar = params[1]
        # xmin = ranges[xvar][0]; xmax = ranges[xvar][1]
        # ymin = ranges[yvar][0]; ymax = ranges[yvar][1]
        xmin = -10; xmax = 10 #Why range impact plot ?
        ymin = -10; ymax = 10

        c.cd()
        expr = '2*(deltaNLL-{}):{}:{}>>grid({},{},{},{},{},{})'.format(minZ, xvar, yvar, nbins, xmin, xmax, nbins, ymin, ymax)
        tree.Draw(expr, '', 'prof colz')

        original = ROOT.TProfile2D(c.GetPrimitive('grid'))

        h_contour = ROOT.TProfile2D('h_contour', 'h_contour', nbins, xmin, xmax, nbins, ymin, ymax)
        h_contour = original.Clone('h_contour')

        #-- Adjust scale so that the best bin has content 0
        best2DeltaNLL = original.GetMinimum()
        for xbin in range(original.GetNbinsX()):
            xcoord = original.GetXaxis().GetBinCenter(xbin)
            for ybin in range(original.GetNbinsY()):
                ycoord = original.GetYaxis().GetBinCenter(ybin)
                if original.GetBinContent(1+xbin,1+ybin)==0:
                    h_contour.Fill(xcoord,ycoord,1000)
                if original.GetBinContent(1+xbin,1+ybin)!=0:
                    h_contour.Fill(xcoord,ycoord,original.GetBinContent(1+xbin,1+ybin)-best2DeltaNLL)
                #h_contour.SetBinContent(1+xbin,1+ybin,original.GetBinContent(1+xbin,1+ybin)-best2DeltaNLL)

        #-- Exclude data outside of the contours
        #h_contour.SetMaximum(11.83)
        #h_contour.SetContour(200)
        #h_contour.GetZaxis().SetRangeUser(0,21);
        h_contour.GetXaxis().SetRange(1,h_contour.GetNbinsX()-3)
        h_contour.GetYaxis().SetRange(1,h_contour.GetNbinsY()-3)

        #-- Set Contours
        c68 = self.ContourHelper.GetContour(h_contour,2.30)
        # print('c68', c68)
        # print('c68[0]', c68[0])
        # print('c68[0].GetN()', c68[0].GetN())
        # for ipt in range(c68[0].GetN()):
        #     x, y = ROOT.Double(0), ROOT.Double(0) #Necessary to pass by reference in GetPoint()
        #     c68[0].GetPoint(ipt, x, y)
        #     print(x, y)

        c95 = self.ContourHelper.GetContour(h_contour,6.18)
        c997 = self.ContourHelper.GetContour(h_contour,11.83)
        c681D = self.ContourHelper.GetContour(h_contour,1.00)
        c951D = self.ContourHelper.GetContour(h_contour,4.00)
        c9971D = self.ContourHelper.GetContour(h_contour,9.00)
        self.ContourHelper.styleMultiGraph(c68,ROOT.kYellow+1,3,1)
        self.ContourHelper.styleMultiGraph(c95,ROOT.kCyan-2,3,1)
        self.ContourHelper.styleMultiGraph(c997,ROOT.kBlue-2,3,1)

        #place holders for the legend, since TLine is weird
        hc68 = ROOT.TH1F('c68', 'c68', 1, 0, 1)
        hc95 = ROOT.TH1F('c95', 'c68', 1, 0, 1)
        hc997 = ROOT.TH1F('c997', 'c68', 1, 0, 1)
        hc68.SetLineColor(ROOT.kYellow+1)
        hc95.SetLineColor(ROOT.kCyan-2)
        hc997.SetLineColor(ROOT.kBlue-2)
        self.ContourHelper.styleMultiGraph(c681D,ROOT.kYellow+1,1,3)
        self.ContourHelper.styleMultiGraph(c951D,ROOT.kCyan-2,1,3)
        self.ContourHelper.styleMultiGraph(c9971D,ROOT.kBlue-2,1,3)

        # Marker for SM point
        marker_1 = ROOT.TMarker()
        marker_1.SetMarkerSize(2.0)
        marker_1.SetMarkerColor(97)
        marker_1.SetMarkerStyle(33)
        marker_2 = ROOT.TMarker()
        marker_2.SetMarkerSize(1.2)
        marker_2.SetMarkerColor(89)
        marker_2.SetMarkerStyle(33)
        hSM = ROOT.TH1F('SM', 'SM', 1, 0, 1)
        hSM.SetMarkerStyle(33)
        hSM.SetMarkerColor(97)

        # h_contour.SetStats(0)
        #h_contour.SetTitle("Significance Contours")
        # h_contour.SetTitle("")
        # h_contour.GetYaxis().SetTitle(self.texdic[params[0].rstrip('i')])
        # h_contour.GetXaxis().SetTitle(self.texdic[params[1].rstrip('i')])

        #-- Labels
        h_contour.SetTitle("")
        h_contour.SetMarkerStyle(8)
        h_contour.SetMarkerSize(1)
        h_contour.GetXaxis().SetTitle(Get_Parameter_LegName(params[0]))
        h_contour.GetYaxis().SetTitle(Get_Parameter_LegName(params[1]))

        # h_contour.GetXaxis().SetTitleOffset(1.1)
        # h_contour.GetXaxis().SetTitleSize(0.04)
        # h_contour.GetXaxis().SetLabelSize(0.04)
        # h_contour.GetYaxis().SetTitleOffset(1.1)
        # h_contour.GetYaxis().SetTitleSize(0.04)
        # h_contour.GetXaxis().SetLabelSize(0.04)
        #h_contour.GetYaxis().SetNdivisions(7)

        h_contour.Draw('AXIS')
        c68.Draw('L SAME')
        # c68[0].Draw('hist SAME')
        c95.Draw('L SAME')
        c997.Draw('L SAME')
        #C681D.Draw('L SAME')
        #C951D.Draw('L SAME')
        #C9971D.Draw('L SAME')
        marker_1.DrawMarker(0,0)
        marker_2.DrawMarker(0,0)

        #c = [2.3, 6.18, 11.83]
        #original.SetContourLevel(0, c[0])
        #original.SetContourLevel(1, c[1])
        #original.SetContourLevel(2, c[2])
        #import numpy
        #original.SetContour(3, numpy.array(c))
        #original.SetMaximum(3)
        #ROOT.gStyle.SetOptStat(0)
        #original.Draw("cont1z")
        #marker_1.DrawMarker(0,0)
        #marker_2.DrawMarker(0,0)

        legend = ROOT.TLegend(0.1,0.9,0.45,0.945)
        legend.AddEntry(hc68,"1#sigma",'l')
        legend.AddEntry(hc95,"2#sigma",'l')
        legend.AddEntry(hc997,"3#sigma",'l')
        legend.AddEntry(hSM,"SM value",'p')
        legend.SetTextSize(0.035)
        legend.SetNColumns(4)
        # legend.Draw('same')

        self.CMS_text.Draw('same')
        if paper==False: self.extraText.Draw('same')
        self.lumiText.Draw('same')
        c.SetGrid()

        outname = '{}_2Dcontour'.format(name.replace('.','_'))
        if paper: outname+= '_paper'
        outname+= '.png'

        c.Print(outname)

        # Save contour to histogram file
        outfile = ROOT.TFile(self.histosFileName,'UPDATE')
        h_contour.Write()
        outfile.Close()

        return


  ####   ####  #####  #####  ###### #        ##   ##### #  ####  #    #  ####
 #    # #    # #    # #    # #      #       #  #    #   # #    # ##   # #
 #      #    # #    # #    # #####  #      #    #   #   # #    # # #  #  ####
 #      #    # #####  #####  #      #      ######   #   # #    # #  # #      #
 #    # #    # #   #  #   #  #      #      #    #   #   # #    # #   ## #    #
  ####   ####  #    # #    # ###### ###### #    #   #   #  ####  #    #  ####

    def PlotCorrelationMatrix(self, name='', allNuisances=True, SM=True):
        '''
        Plot correlation matrix.

        allNuisances: True <-> plot full correlation matrix (NB: labels invisible if lots of nuisances); else, only plot correlations between POIs (useful only for n-D fits)
        '''

        logging.info(colors.fg.lightblue + "Enter function PlotCorrelationMatrix()\n" + colors.reset)

        filepath = './multidimfit{}.root'.format(name)
        if not os.path.exists(filepath):
            logging.error("File " + filepath + " does not exist!".format(name))
            return

        ROOT.gROOT.SetBatch(True)
        canvas = ROOT.TCanvas()
        canvas.Clear()
        ROOT.gStyle.SetPaintTextFormat('.2f')

        #-- Get rooFitResult object
        logging.info("Reading corrMatrix from file: {}".format(filepath))
        rooFitFile = ROOT.TFile.Open(filepath)
        rooFitResult = rooFitFile.Get('fit_mdf')

        #-- Get correlation matrix
        rooFitResult.correlationHist().Draw('colz')
        matrix = canvas.GetPrimitive('correlation_matrix')
        matrix.SetName("corrMatrix")

        params = opts['wcs']
        if SM==True: params = opts['SM_mus']
        npar = len(params)

        if allNuisances == False: #Will restrict plot only to POIs #NB: only makes sense for nD fits, not 1D
            BinsX, BinsY = [], []
            for idx,label in enumerate(matrix.GetXaxis().GetLabels()):
                if label in params: BinsX.append(1+idx)
            for idy,label in enumerate(matrix.GetYaxis().GetLabels()):
                if label in params: BinsY.append(1+idy)
            newmatrix = ROOT.TH2D("Correlation Matrix","Correlation Matrix",npar,0,npar,npar,0,npar)
            newmatrix.SetTitle("Correlation Matrix")
            for idx,binx in enumerate(BinsX):
                for idy,biny in enumerate(BinsY):
                    newmatrix.SetBinContent(1+idx,npar-idy,matrix.GetBinContent(binx,matrix.GetNbinsY()-biny+1))
                newmatrix.GetXaxis().SetBinLabel(1+idx,matrix.GetXaxis().GetBinLabel(binx))
            for idy,biny in enumerate(BinsY):
                newmatrix.GetYaxis().SetBinLabel(npar-idy,matrix.GetYaxis().GetBinLabel(matrix.GetNbinsY()-biny+1))

            #-- Change plot format
            newmatrix.SetMaximum(1)
            newmatrix.SetMinimum(-1.)
            newmatrix.SetStats(0)
            newmatrix.SetName("corrMatrixSM")
            newmatrix.SetTitle("Correlation Matrix")
            newmatrix.Draw('colz')
            newmatrix.Draw('same text') #Also draw values in bins (only possible if few bins)

        else: #Draw default matrix
            matrix.Draw('colz')

        # Save the plot
        outname = 'correlationMatrix{}.png'.format(name.replace('.','_'))
        canvas.Print(outname)

        #-- Save the plot to the histogram file
        # outfile = ROOT.TFile(self.histosFileName,'UPDATE')
        # newmatrix.Write()
        # outfile.Close()

        return


                         #####  ######  #     #  #####
  ####  #####  # #####  #     # #     # #  #  # #     #
 #    # #    # # #    #       # #     # #  #  # #
 #      #    # # #    #  #####  #     # #  #  # #
 #  ### #####  # #    # #       #     # #  #  # #
 #    # #   #  # #    # #       #     # #  #  # #     #
  ####  #    # # #####  ####### ######   ## ##   #####

    '''
    def grid2DWC(self, name='', wc='', fits=[]):

        ROOT.gROOT.SetBatch(True)
        canvas = ROOT.TCanvas()

        # Get limit tree
        fit_file = ROOT.TFile.Open('./higgsCombine{}.{}.MultiDimFit.root'.format(name,wc))
        limit_tree = fit_file.Get('limit')

        def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
            return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

        # Fined event with minimum NLL
        min_val = 0
        wc_values = []
        for entry in range(limit_tree.GetEntries()):
            limit_tree.GetEntry(entry)
            #min_values.append((entry,2*limit_tree.GetLeaf('deltaNLL').GetValue(0)))
            if limit_tree.GetLeaf(wc).GetValue(0) == fits[wc]: min_val = entry

        # Load other 15 WCs in the minimum entry
        best_vals = []
        for w in self.wcs:
            limit_tree.GetEntry(min_val)
            if w == wc:
                best_vals.append([w, str(limit_tree.GetLeaf(w).GetValue(0))])
            else:
                best_vals.append([w, str(limit_tree.GetLeaf('trackedParam_' + w).GetValue(0))])

        # Close files
        fit_file.Close()
        return best_vals

    def batchGrid2DWC(self, name=''):
        best = []

        fits_float = self.getIntervalFits(name)
        fits = {lst[0] : lst[1] for lst in fits_float}

        for wc in self.wcs:
            best.append(self.grid2DWC(name, wc, fits))

        j=0
        for i in range(0, len(self.wcs)):
            if i ==0: print '  '.join(self.wcs)
            lst = [str(round(float(w[1]), 3)) for w in best[i]]
            print best[i][j][0], ' '.join(lst)
            j = j + 1
    '''

 # #    # ##### ###### #####  #    #   ##   #       ####
 # ##   #   #   #      #    # #    #  #  #  #      #
 # # #  #   #   #####  #    # #    # #    # #       ####
 # #  # #   #   #      #####  #    # ###### #           #
 # #   ##   #   #      #   #   #  #  #    # #      #    #
 # #    #   #   ###### #    #   ##   #    # ######  ####

    '''
    def getIntervalFits(self, basename='.EFT.SM.Float', params=[], siginterval=2):
        ### Return a table of parameters, their best fits, and their uncertainties ###
        ### Use 1D scans instead of regular MultiDimFit ###
        if not params:
            params = self.wcs


        ROOT.gROOT.SetBatch(True)

        fit_array = [] # List of [WC, WC value of minimum, [2sig lowedges], [2sig highedges]]

        for param in params:

            # Get scan TTree
            logging.debug("Obtaining result of scan: higgsCombine{}.{}.MultiDimFit.root".format(basename,param))
            fit_file = ROOT.TFile.Open('./higgsCombine{}.{}.MultiDimFit.root'.format(basename,param))
            limit_tree = fit_file.Get('limit')

            # Extract points
            wc_values = []
            nll_values = []
            for entry in range(limit_tree.GetEntries()):
                limit_tree.GetEntry(entry)
                wc_values.append(limit_tree.GetLeaf(param).GetValue(0))
                nll_values.append(2*limit_tree.GetLeaf('deltaNLL').GetValue(0))

            # Rezero deltanll values
            bestNLL = min(nll_values)
            logging.debug("Best nll value is {}".format(bestNLL))
            logging.debug("nll_values:")
            logging.debug(nll_values)
            nll_values = [oldValue-bestNLL for oldValue in nll_values]

            # Sort values just in case
            coords = zip(wc_values,nll_values)
            coords.sort(key = lambda t: t[0])
            wc_values, nll_values = zip(*coords)
            wc_values = numpy.asarray(wc_values)
            nll_values = numpy.asarray(nll_values)

            # Prep a `spline to get the exact crossings of the 1,2 sigma levels
            graph = ROOT.TGraph()
            graph = ROOT.TGraph(len(coords), wc_values, nll_values)
            spline = ROOT.TSpline3("spline3", graph)

            #f1 = ROOT.TF1('f1','poln')
            #graph.Fit('poln','N')
            #fitfunc = graph.GetFunction('poln')

            # Extract 2-sig certainty intervals and save WC value of minumum
            lowedges=[]
            l1sigma=[]
            highedges=[]
            h1sigma=[]
            true_minimum = -1000
            prevnll = 1000
            prevnll1 = 1000
            for idx,coord in enumerate(coords):
                wc,nll = coord[0],coord[1]
                # Did we cross a low edge?
                if prevnll>4 and 4>nll:
                    #cross = fitfunc.GetX(4, graph.GetX()[idx-1], graph.GetX()[idx])
                    interval = prevnll-nll
                    linPctInterp = (prevnll-4)/interval
                    cross = graph.GetX()[idx-1]+(graph.GetX()[idx]-graph.GetX()[idx-1])*linPctInterp
                    lowedges.append(cross)
                # Did we cross a high edge?
                if prevnll<4 and 4<nll:
                    #cross = fitfunc.GetX(4, graph.GetX()[idx-1], graph.GetX()[idx])
                    interval = nll-prevnll
                    linPctInterp = (4-prevnll)/interval
                    cross = graph.GetX()[idx-1]+(graph.GetX()[idx]-graph.GetX()[idx-1])*linPctInterp
                    highedges.append(cross)
                # Is this the best fit?
                if prevnll>1 and 1>nll:
                    #cross = fitfunc.GetX(2, graph.GetX()[idx-1], graph.GetX()[idx])
                    interval = prevnll-nll
                    linPctInterp = (prevnll-2)/interval
                    cross = graph.GetX()[idx-1]+(graph.GetX()[idx]-graph.GetX()[idx-1])*linPctInterp
                    l1sigma.append(cross)
                # Did we cross a high edge?
                if prevnll<1 and 1<nll:
                    #cross = fitfunc.GetX(2, graph.GetX()[idx-1], graph.GetX()[idx])
                    interval = nll-prevnll
                    linPctInterp = (1-prevnll)/interval
                    cross = graph.GetX()[idx-1]+(graph.GetX()[idx]-graph.GetX()[idx-1])*linPctInterp
                    h1sigma.append(cross)
                # Is this the best fit?
                if nll == min(nll_values):
                    true_minimum = wc
                # Continue
                prevnll = nll
            if not len(lowedges) == len(highedges):
                logging.error("Something is strange! Interval is missing endpoint!")
            if not len(l1sigma) == len(h1sigma):
                logging.error("Something is strange! Interval is missing endpoint!")
            ## uncomment for 2 decimal place printing for AN
            #true_minimum = '%.2f' % float(true_minimum)
            #lowedges = ['%.2f' % elem for elem in lowedges]
            #highedges = ['%.2f' % elem for elem in highedges]
            if siginterval==2: fit_array.append([param,true_minimum,lowedges,highedges])
            elif siginterval==1: fit_array.append([param,true_minimum,l1sigma,h1sigma])
            else: fit_array.append([param,true_minimum,lowedges,highedges])

        for line in fit_array:
            pline = line[:]
            if pline[0][-1] == 'i': pline[0] = pline[0][:-1]
            pline[0] = '\\' + pline[0] + '$/\\Lambda^{2}$'
            pline[0] = pline[0].replace('3','a')
            #print line
            one = pline[2]
            one = ['%.2f' % elem for elem in one]
            two = pline[3]
            two = ['%.2f' % elem for elem in two]
            s = pline[0] + ' & '
            if len(one)==2:
                one = ', '.join(one)
                two = ', '.join(two)
                s += '[' + str(one) + ']' + ' and [' + str(two) + ']'
            else:
                s += '[' + str(one[0]) + ', ' + str(two[0]) + ']'
            print s

        return fit_array
    '''


 #######  ######## ##     ## ######## ########   ######
##     ##    ##    ##     ## ##       ##     ## ##    ##
##     ##    ##    ##     ## ##       ##     ## ##
##     ##    ##    ######### ######   ########   ######
##     ##    ##    ##     ## ##       ##   ##         ##
##     ##    ##    ##     ## ##       ##    ##  ##    ##
 #######     ##    ##     ## ######## ##     ##  ######

def Get_nSigmaBand_Graphs(self, graph):
    #-- y coordinates corresponding to relevant NLL thresholds
    yval_68 = 1 #68% CL <-> 1 sigma <-> nll=1
    yval_95 = 3.84 #95% CL <-> ~2 sigma <-> nll=3.84
    # yval_2sigmas = 4 #~95.45% CL <-> 2 sigma <-> nll=4

    #-- Display vertical lines at each intersection of the NLL function with the 95% CL threshold
    list_X_intersects95 = Get_Intersection_X(graph, yval_95)
    lines_X_intersects95 = []
    if len(list_X_intersects95) > 0:
        for i,x in enumerate(list_X_intersects95):
            lines_X_intersects95.append(ROOT.TLine(x,0,x,yval_95))
            # lines_X_intersects95[i].Draw('same') #Draw last
            lines_X_intersects95[i].SetLineColor(ROOT.kAzure-1) #ROOT.kRed+1
            lines_X_intersects95[i].SetLineWidth(3)
            lines_X_intersects95[i].SetLineStyle(self.linestyle)

        # Create and draw filled TGraph (shaded area below NLL function representing 2sigmas area)
        # fgraph95 = ROOT.TGraph(graph.GetN(),graph.GetX(),graph.GetY())
        fgraph95 = graph.Clone()
        ipt = 0
        while ipt<fgraph95.GetN():
        # for ipt in range(fgraph95.GetN()):
            x, y = ROOT.Double(0), ROOT.Double(0) #Necessary to pass by reference in GetPoint()
            fgraph95.GetPoint(ipt, x, y)
            # print(x, y)
            if len(list_X_intersects95)>=2 and (x<=list_X_intersects95[0] or x>=list_X_intersects95[1]):
                # print(ipt, x, y)
                fgraph95.RemovePoint(ipt)
            else: ipt+= 1

        # Trick: manually set first and last points to get the filled area properly defined
        if len(list_X_intersects95)>=2: #Only draw filled area if there are multiple intersections (else, need to care about direction)
            fgraph95.SetPoint(0, list_X_intersects95[0], 0) #Add first point at y=0
            fgraph95.InsertPointBefore(1, list_X_intersects95[0], yval_95) #Add second point at y=Y
            fgraph95.SetPoint(fgraph95.GetN(), list_X_intersects95[1], yval_95) #Add first-to-last point at y=Y
            fgraph95.SetPoint(fgraph95.GetN(), list_X_intersects95[1], 0) #Add last point at y=0

            #-- Printouts
            # print(fgraph95.GetN())
            # for ipt in range(0, fgraph95.GetN()):
            #     x, y = ROOT.Double(0), ROOT.Double(0) #Necessary to pass by reference in GetPoint()
            #     fgraph95.GetPoint(ipt, x, y)
            #     print(ipt, x, y)

            fgraph95.SetFillColorAlpha(ROOT.kRed+1, 0.55)
            fgraph95.SetFillStyle(self.fillstyle)

    #-- Idem for 68% CL -- not drawn for now
    list_X_intersects68 = Get_Intersection_X(graph, yval_68)
    lines_X_intersects68 = []
    if len(list_X_intersects68) > 0:
        for i,x in enumerate(list_X_intersects68):
            lines_X_intersects68.append(ROOT.TLine(x,0,x,yval_68))
            lines_X_intersects68[i].SetLineColor(ROOT.kRed+1)
            lines_X_intersects68[i].SetLineWidth(3)
            lines_X_intersects68[i].SetLineStyle(self.linestyle)

        fgraph68 = ROOT.TGraph(graph.GetN(),graph.GetX(),graph.GetY())
        # fgraph68 = graph.Clone()
        ipt = 0
        while ipt<fgraph68.GetN():
            x, y = ROOT.Double(0), ROOT.Double(0) #Necessary to pass by reference in GetPoint()
            fgraph68.GetPoint(ipt, x, y)
            # print(x, y)
	#print(list_X_intersects68)
            if (len(list_X_intersects68)>0 and x<=list_X_intersects68[0]) or (len(list_X_intersects68)>1 and x>=list_X_intersects68[1]) or (x==0 and y==0):
                # print('ipt', ipt)
                fgraph68.RemovePoint(ipt)
            else: ipt+= 1

        if len(list_X_intersects68)>=2: #Only draw filled area if there are multiple intersections (else, need to care about direction)
            fgraph68.SetPoint(0, list_X_intersects68[0], 0) #Add first point at y=0
            fgraph68.InsertPointBefore(1, list_X_intersects68[0], yval_68) #Add second point at y=Y
            fgraph68.SetPoint(fgraph68.GetN(), list_X_intersects68[1], yval_68) #Add first-to-last point at y=Y
            fgraph68.SetPoint(fgraph68.GetN(), list_X_intersects68[1], 0) #Add last point at y=0

            #-- Printouts
            # print(fgraph68.GetN())
            # for ipt in range(0, fgraph68.GetN()):
            #     x, y = ROOT.Double(0), ROOT.Double(0) #Necessary to pass by reference in GetPoint()
            #     fgraph68.GetPoint(ipt, x, y)
            #     print(ipt, x, y)

            fgraph68.SetFillColor(12)
            fgraph68.SetFillStyle(self.fillstyle)

    return yval_68, list_X_intersects68, fgraph68, lines_X_intersects68, yval_95, list_X_intersects95, fgraph95, lines_X_intersects95











# //--------------------------------------------
# //--------------------------------------------

##     ##    ###    #### ##    ##
###   ###   ## ##    ##  ###   ##
#### ####  ##   ##   ##  ####  ##
## ### ## ##     ##  ##  ## ## ##
##     ## #########  ##  ##  ####
##     ## ##     ##  ##  ##   ###
##     ## ##     ## #### ##    ##

# //--------------------------------------------
# //--------------------------------------------

if __name__ == "__main__":

# User options
# //--------------------------------------------
    SM = False #True <-> consider SM scenario (rather than SMEFT)
    mode = '1D'
    POI = []
    paper = False
    filepath = ''
    name = '' #Suffix added to output filenames


# Set up the command line arguments
# //--------------------------------------------
    parser = argparse.ArgumentParser(description='Perform SM and EFT fits using custom Physics Model')
    parser.add_argument("-m", metavar="mode", help="Plotting mode: 1D, 2D, contour, manual")
    #parser.add_argument("-scan", metavar="scan", help="Scan type (1D, 2D, manual, ...)")
    parser.add_argument('-P','--POI', metavar="POI", nargs='+', help='Define POI(s)', required=False) #Takes >=0 args
    parser.add_argument("--sm", metavar="SM", help="Consider SM scenario (rather than SMEFT)", nargs='?', const=1)
    parser.add_argument("--paper", metavar="paper", help="Make plots for paper (remove Prelim. labels)", nargs='?', const=1)
    parser.add_argument("-f", metavar="file path", help="Path to the rootfile containing the object to plot")
    parser.add_argument("-name", metavar="name", help="add suffix to output filename")

    args = parser.parse_args()
    if len(args.POI) == 0:
    	print('ERROR: no POI selected !')
        exit(1)
    if len(args.POI) == 2:
    	print('2 POIs detected --> Make 2D plot')
        mode = '2D'
        if args.m == 'contour':
            print('Not supported !')
            exit(1)
    elif args.m in ['1D','2D','manual','contour']: mode = args.m
    elif args.m:
        print('ERROR ! Wrong [-m] arg !'); exit(1)
    if args.POI: POI = args.POI
    if args.sm: SM = True
    if args.paper: paper = True
    if args.f: filepath = args.f
    if args.name: name = args.name

    if name != '' and name[0] != '.': name = '.{}'.format(name) #Convention: add a '.' before custom name

    plotter = EFTPlot(opts, mode)
    Load_Canvas_Style()


# SM
# //--------------------------------------------
    if SM:
        if name == '': name = '.SM'

        if mode=='1D': plotter.Plot_NLLscan_1D(param=opts['SM_mu'], log=False, paper=paper, filepath=filepath, name=name)
        elif mode=='2D': plotter.Plot_NLLscan_2D(mode='SM', params=opts['SM_mus'], ceiling=100, log=False, paper=paper, filepath=filepath, name=name)
        elif mode=='contour': plotter.ContourPlotEFT(mode='SM', params=opts['SM_mus'], paper=paper, filepath=filepath, name=name)

# SMEFT
# //--------------------------------------------
    else:
        if name == '': name = '.EFT'

        if mode=='1D':
            param_tmp = POI if len(POI) == 1 else [opts['wc']]
            plotter.Plot_NLLscan_1D(param=param_tmp[0], log=False, paper=paper, filepath=filepath, name=name)
        elif mode=='2D':
            param_tmp = POI if len(POI) == 2 else [opts['wcs_pairs']]
            plotter.Plot_NLLscan_2D(mode='EFT', params=param_tmp, ceiling=100, log=False, paper=paper, filepath=filepath, name=name)
        elif mode=='contour':
            param_tmp = POI if len(POI) == 2 else [opts['wcs_pairs']]
            plotter.ContourPlotEFT(mode='EFT', params=param_tmp, paper=paper, filepath=filepath, name=name)
        elif mode == 'manual':
            param_tmp = POI if len(POI) == 1 else [opts['wc']]
            plotter.Plot1DManualNLLScan(param=param_tmp[0], paper=paper, name=name)

# OTHER
# //--------------------------------------------
    # plotter.PlotCorrelationMatrix(name=name, allNuisances=True, SM=SM)
