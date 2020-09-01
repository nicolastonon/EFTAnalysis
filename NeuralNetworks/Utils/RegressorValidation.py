#Specific validation functions for regressor NN

import ROOT
import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from Utils.ColoredPrintout import colors
from Utils.Helper import *

                                    #     #  #####
 #####  #####  ###### #####         #     # #     #    ##### #####  #    # ##### #    #
 #    # #    # #      #    #        #     # #            #   #    # #    #   #   #    #
 #    # #    # #####  #    #        #     #  #####       #   #    # #    #   #   ######
 #####  #####  #      #    # ###     #   #        #      #   #####  #    #   #   #    #
 #      #   #  #      #    # ###      # #   #     #      #   #   #  #    #   #   #    #
 #      #    # ###### #####  ###       #     #####       #   #    #  ####    #   #    #

# //--------------------------------------------
# //--------------------------------------------

def Plot_LR_Pred_vs_Truth(opts, list_features, list_labels, list_yTrain_allClasses, list_yTest_allClasses, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_truth_Test_allClasses, list_xTrain_allClasses, list_xTest_allClasses, weight_dir):
    """
    Make validation plots for regressor NN using test data. Compare predictions to true target values.
    """

    if opts["strategy"] not in ["regressor", "ROLR", "RASCAL"]: return #Useful for regressor validation only

    nodename='r'
    # if opts['regress_onLogr'] == True: nodename='log(r)'

    xmin=-0.; xmax=5
    ymin=-0.; ymax=5
    # xmin=-0.; xmax=50
    # ymin=-0.; ymax=50

    mycol = 'g'

# //--------------------------------------------
#DATA

    # print(list_yTest_allClasses[0].shape); print(list_yTest_allClasses[1].shape)
    if opts["nofOutputNodes"] == 1:
        truth_data = np.squeeze(np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1])))
    else:
        truth_data = np.squeeze(np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1]))[:,0])

    if len(list_predictions_test_allNodes_allClasses[0]) > 1: #Concatenate 'sig' and 'bkg'
        pred_data = np.squeeze(np.concatenate((list_predictions_test_allNodes_allClasses[0][0],list_predictions_test_allNodes_allClasses[0][1])))
        class_data = np.squeeze(np.concatenate((list_truth_Test_allClasses[0],list_truth_Test_allClasses[1])))
    else: #For simple regressor, there may be only 1 process class to consider
        pred_data = np.squeeze(list_predictions_test_allNodes_allClasses[0][0])
        class_data = np.squeeze(list_truth_Test_allClasses[0])

    if opts["comparVarIdx"] >= 0: #Also plot kinReco results
        truth_data = np.squeeze(np.concatenate((np.squeeze(list_yTest_allClasses[0]), list_xTest_allClasses[0][:,opts["comparVarIdx"]])))
        pred_data = np.squeeze(np.concatenate((list_predictions_test_allNodes_allClasses[0][0],list_predictions_test_allNodes_allClasses[0][0])))
        class_data = np.squeeze(np.concatenate((list_truth_Test_allClasses[0],list_truth_Test_allClasses[0]-1))) #Set different (arbitrary) class values for predictions / comparison variable

    truth_data, pred_data, class_data = unison_shuffled_copies(truth_data, pred_data, class_data) #Randomize events (else, only see latest events on scatteplot)

    # print(truth_data.shape); print(pred_data.shape); print(class_data.shape)
    # print(truth_data[:20]); print(pred_data[:20])

    # quantile = np.quantile(pred_data, 0.995)
    quantile = np.quantile(truth_data, 0.995) #Define axis range based on population of truth values
    xmax=abs(quantile)
    ymax=xmax #Need to have same axes range, so that diagonal represents perfect correlation

# //--------------------------------------------
#Scatterplot
#See https://seaborn.ppred_data.org/tutorial/distributions.html

    fig = plt.figure('splot', figsize=(10, 10))
    # timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
    # timer.add_callback(close_event)
    if opts["strategy"] in ["ROLR","RASCAL"]:
        plt.title('Predicted VS True '+nodename)
        plt.xlabel(r'True '+nodename+r'(x|$\theta_0,\theta_1$)', fontsize=15) # add 'r' in front <-> interpreted as raw string
        plt.ylabel(r'Learned '+nodename+r'(x|$\theta_0,\theta_1$)', fontsize=15) # add 'r' in front <-> interpreted as raw string #color='darkorange'
    else:
        plt.title('Predictions VS Truth')
        plt.xlabel('Truth', fontsize=15)
        plt.ylabel('Prediction', fontsize=15)

    nmax=2000 #Can't see if there are too many points
    splot = sns.scatterplot(x=truth_data[:nmax], y=pred_data[:nmax], hue=class_data[:nmax], style=class_data[:nmax])
    leg_handles = splot.get_legend_handles_labels()[0]
    splot.legend(leg_handles, ['EFT', 'SM']) #NB: EFT first because class label is 0
    if opts["comparVarIdx"] >= 0: splot.legend(leg_handles, [list_features[opts["comparVarIdx"]], 'NN'])
    ax = fig.gca()
    ax.set(xlim=(xmin, xmax))
    ax.set(ylim=(ymin, ymax))

    # Plot your initial diagonal line based on the starting
    # xlims and ylims.
    diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")

    # timer.start()
    # plt.show()

    plotname = weight_dir + 'ScatterPlotLR_PredvsTruth.png'
    fig.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR scatter plot as :", colors.reset, plotname)
    fig.clear(); plt.close('splot')

# //--------------------------------------------
#Regplot

    """
    fig = plt.figure('regplot', figsize=(10, 10))
    ax = fig.gca(); ax.set(xlim=(xmin, xmax)); ax.set(ylim=(xmin, xmax))
    regplot = sns.regplot(x=truth_data, y=pred_data, color=mycol)
    plotname = weight_dir + 'RegPlot_PredvsTruth.png'
    fig.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR reg plot as :", colors.reset, plotname)
    fig.clear(); plt.close('regplot')
    """

# //--------------------------------------------
#Jointplot

    #NB: 'jointplot' func creates its own figure (as do: FacetGrid, factorplot, lmplot, PairGrid, pairplot, JointGrid, jointplot.)
    jplot = sns.jointplot(x=truth_data, y=pred_data, color=mycol, alpha=0.5, xlim=(xmin,xmax), ylim=(ymin,ymax))
    # jplot.ax_marg_x.set_xlim(xmin, xmax); jplot.ax_marg_y.set_ylim(ymin, ymax)
    plotname = weight_dir + 'JointPlot_PredvsTruth.png'
    jplot.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR joint plot as :", colors.reset, plotname)

# //--------------------------------------------
#Jointplot HEX

    """
    hexplot = sns.jointplot(x=truth_data, y=pred_data, color=mycol, alpha=0.5, xlim=(xmin,xmax), ylim=(ymin,ymax), kind="hex")
    # hexplot.ax_marg_x.set_xlim(xmin, xmax); jplot.ax_marg_y.set_ylim(ymin, ymax)
    plotname = weight_dir + 'JointPlotHEX_PredvsTruth.png'
    hexplot.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR joint HEX plot as :", colors.reset, plotname)
    """

# //--------------------------------------------
#Jointplot KDE

    kdeplot = sns.jointplot(x=truth_data, y=pred_data, color=mycol, alpha=0.5, xlim=(xmin,xmax), ylim=(ymin,ymax), kind="kde")
    # kdeplot.ax_marg_x.set_xlim(xmin, xmax); jplot.ax_marg_y.set_ylim(ymin, ymax)
    plotname = weight_dir + 'JointPlotKDE_PredvsTruth.png'
    kdeplot.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR joint KDE plot as :", colors.reset, plotname)

# //--------------------------------------------
#Jointplot KDE+scatterpol

    """
    kdesplot = sns.jointplot(x=truth_data, y=pred_data, color=mycol, alpha=0.5, xlim=(xmin,xmax), ylim=(ymin,ymax), kind="kde")
    kdesplot.plot_joint(plt.scatter, c=mycol, s=30, linewidth=1, marker="+")
    kdesplot.ax_joint.collections[0].set_alpha(0)
    # kdesplot.ax_marg_x.set_xlim(xmin, xmax); jplot.ax_marg_y.set_ylim(ymin, ymax)
    # kdesplot.set_axis_labels("$X$", "$Y$");
    plotname = weight_dir + 'JointPlotKDESPLOT_PredvsTruth.png'
    kdesplot.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR joint KDE+SPLOT plot as :", colors.reset, plotname)
    """

# //--------------------------------------------
#Scatterplot
#See https://seaborn.ppred_data.org/tutorial/distributions.html

    """
    fig = plt.figure('splot2', figsize=(10, 10))
    plt.title('')
    plt.xlabel(list_features[0], fontsize=15)
    plt.ylabel(list_features[1], fontsize=15)
    # plt.zlabel('xxx', fontsize=15)
    x1data = np.concatenate((list_xTest_allClasses[0][:500,0],list_xTest_allClasses[1][:500,0]))
    x2data = np.concatenate((list_xTest_allClasses[0][:500,1],list_xTest_allClasses[1][:500,1]))
    splot2 = plt.scatter(x1data, x2data, c=truth_data[:len(x1data)])
    plt.colorbar()
    # ax = fig.gca()
    # ax.set(xlim=(xmin, xmax))
    # ax.set(ylim=(ymin, ymax))
    plotname = weight_dir + 'ScatterPlot2LR_PredvsTruth.png'
    fig.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR scatter plot 2 as :", colors.reset, plotname)
    fig.clear(); plt.close('splot2')
    """

# //--------------------------------------------
#Scatterplot

    if opts["parameterizedNN"] == True:
        fig = plt.figure('splot3', figsize=(10, 10))
        # plt.title('Predicted VS True '+nodename)
        plt.xlabel(r'$\theta$', fontsize=15) # add 'r' in front <-> interpreted as raw string
        plt.ylabel(r'True '+nodename+r'(x|$\theta_0,\theta_1$)', fontsize=15) # add 'r' in front <-> interpreted as raw string #color='darkorange'
        xdata = np.squeeze(list_xTest_allClasses[1][:1000,-1])
        if opts["nofOutputNodes"] == 1: ydata = np.squeeze(list_yTest_allClasses[1][:1000])
        else: ydata = np.squeeze(list_yTest_allClasses[1][:1000,0])
        splot3 = sns.scatterplot(x=xdata, y=ydata)
        # ax = fig.gca()
        # ax.set(xlim=(xmin, xmax))
        # ax.set(ylim=(ymin, ymax))
        plotname = weight_dir + 'ScatterPlot3LR_PredvsTruth.png'
        fig.savefig(plotname)
        print(colors.fg.lightgrey, "\nSaved LR scatter plot3 as :", colors.reset, plotname)
        fig.clear(); plt.close('splot3')

    return


# //--------------------------------------------
# //--------------------------------------------

 #####  #    # #      #
 #    # #    # #      #
 #    # #    # #      #
 #####  #    # #      #
 #      #    # #      #
 #       ####  ###### ######

def Make_Pull_Plot(opts, weight_dir, list_yTest_allClasses, list_predictions_test_allNodes_allClasses, list_truth_Test_allClasses, list_PhysicalWeightsTest_allClasses, list_xTest_allClasses):
    """
    Plot difference between prediction and truth.
    """

    if opts["strategy"] not in ["regressor", "ROLR", "RASCAL"]: return #Only useful for regressors

# //--------------------------------------------
#DATA

    # print(list_yTest_allClasses[0].shape); print(list_yTest_allClasses[1].shape)
    if opts["nofOutputNodes"] == 1:
        truth_data = np.squeeze(np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1])))
    else:
        truth_data = np.squeeze(np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1]))[:,0])

    if len(list_predictions_test_allNodes_allClasses[0]) > 1: #Concatenate 'sig' and 'bkg'
        pred = np.squeeze(np.concatenate((list_predictions_test_allNodes_allClasses[0][0],list_predictions_test_allNodes_allClasses[0][1])))
        class_data = np.squeeze(np.concatenate((list_truth_Test_allClasses[0],list_truth_Test_allClasses[1])))
    else: #For simple regressor, there may be only 1 process class to consider
        pred = np.squeeze(list_predictions_test_allNodes_allClasses[0][0])
        class_data = np.squeeze(list_truth_Test_allClasses[0])

    if opts["comparVarIdx"] >= 0: #Also plot kinReco results
        truth_data = np.squeeze(np.concatenate((np.squeeze(list_yTest_allClasses[0]),np.squeeze(list_yTest_allClasses[0]))))
        pred = np.squeeze(np.concatenate((list_predictions_test_allNodes_allClasses[0][0],list_xTest_allClasses[0][:,opts["comparVarIdx"]])))
        class_data = np.squeeze(np.concatenate((list_truth_Test_allClasses[0],list_truth_Test_allClasses[0]-1))) #Set different (arbitrary) class values for predictions / comparison variable

# //--------------------------------------------

    #Transform classifier -> LR
    # if opts["strategy"] is "classifier" or opts["strategy"] is "CARL":
    #     truth = r_from_s(truth)
    #     pred = r_from_s(pred)

    nbins = 30
    xmin = -2; xmax= 2
    norm = True #True <-> normalize histos to unity

    hpull = TH1F('Pull', 'Pull', nbins, xmin, xmax); hpull.Sumw2(); hpull.SetDirectory(0)
    for idx in range(0, len(list_yTest_allClasses[0])):
        # print('1 idx', idx)
        tmp = (pred[idx] - truth_data[idx]) / truth_data[idx]
        hpull.Fill(tmp, 1.)
        # print('pred ', pred[idx], 'truth ', truth_data[idx], ' => ', tmp)

    hpull.SetFillColor(18)
    hpull.SetLineColor(1)
    if norm: hpull.Scale(1./hpull.Integral())

    c = ROOT.TCanvas()
    hpull.Draw("hist")

    if opts["strategy"] is "regressor" and opts["comparVarIdx"] >= 0: #Superimpose pull plot for comparison variable
        hpull_comp = TH1F('Pull', 'Pull', nbins, xmin, xmax); hpull_comp.Sumw2(); hpull_comp.SetDirectory(0)
        for idx in range(len(list_yTest_allClasses[0]), len(pred)):
            # print('2 idx', idx)
            tmp = (pred[idx] - truth_data[idx]) / truth_data[idx]
            hpull_comp.Fill(tmp, 1.)
            # print('comp ', pred[idx], 'truth ', truth_data[idx], ' => ', tmp)

        # hpull_comp.SetFillColor(18)
        hpull_comp.SetLineColor(2)
        if norm: hpull_comp.Scale(1./hpull_comp.Integral())
        hpull_comp.Draw("hist sames") #'sames': Same as "SAME" and draw the statistics box

    c.Update()
    plotname = weight_dir + "Pull_plot.png"
    c.SaveAs(plotname)
    print(colors.fg.lightgrey, "\nSaved pull plot as :", colors.reset, plotname)

    return


  ####  ##### #    # ###### #####   ####
 #    #   #   #    # #      #    # #
 #    #   #   ###### #####  #    #  ####
 #    #   #   #    # #      #####       #
 #    #   #   #    # #      #   #  #    #
  ####    #   #    # ###### #    #  ####

# Adapted from code by Sebastian
#TBD -- adapt ranges, etc.
def doEvaluationPlots(yTest, yPredicted, weightTest, weight_dir):

    if not os.path.exists(weight_dir):
        os.makedirs(weight_dir)

    #-- Truth/pred distributions
    # f = plt.figure()
    # values, bins, patches = plt.hist(yTest, bins=30, range=(0,1), label="true", alpha=0.5)
    # values2, bins2, patches2 = plt.hist(yPredicted, bins=30, range=(0,1), label="reco", alpha=0.5)
    # plt.legend(loc = "best")
    # plt.ylabel('Events')
    # plt.xlabel('Variable')
    # f.savefig(weight_dir+"distributions_truePred.png")

    #-- Pred VS Truth and (Truth-Pred) VS Truth scatter plots
    histo = ROOT.TH2F("", "", 500, 0, 1, 1000, -1, 1)
    histo.SetDirectory(0)
    histoRecoGen = ROOT.TH2F("", "", 500, 0, 1, 500, 0, 1)
    histoRecoGen.SetDirectory(0)
    histoRecoGen2 = ROOT.TH2F("", "", 10, 0, 1, 10, 0, 1)
    histoRecoGen2.SetDirectory(0)
    # print (yTest.shape, yPredicted.shape)
    for true, pred, weight in zip(yTest, yPredicted, weightTest):
        diff = true-pred
        histo.Fill(diff, true, weight)
        histoRecoGen.Fill(pred, true, weight)
        histoRecoGen2.Fill(pred, true, weight)

    c = ROOT.TCanvas("c1","c1",800,800)

    histo.GetXaxis().SetRangeUser(350,1500)
    histo.GetYaxis().SetRangeUser(-1200,1200)
    histo.SetStats(0)
    histo.Draw("colz")
    c.SaveAs(weight_dir+"TruthRecoDiff2D.png")
    c.Clear()

    histoRecoGen.Draw("colz")
    corrLatex = ROOT.TLatex()
    corrLatex.SetTextSize(0.65 * corrLatex.GetTextSize())
    corrLatex.DrawLatexNDC(0.65, 0.85, str(np.round(histoRecoGen.GetCorrelationFactor(),3)))
    c.SaveAs(weight_dir+"TruthReco2D.png")
    c.Clear()

    #-- Proportion of events per 2D truth/pred bin
    histoRecoGen2.Scale(1./histoRecoGen2.Integral())
    histoRecoGen2.Scale(100.)
    ROOT.gStyle.SetPaintTextFormat("1.2f");
    histoRecoGen2.SetStats(0)
    histoRecoGen2.SetMarkerSize(1)
    histoRecoGen2.Draw("colz text")
    c.SaveAs(weight_dir+"Response.png")
    c.Clear()

    # style.style1d()
    # s = style.style1d()
    c=ROOT.TCanvas()
    Xnb=20
    Xr1=0.
    Xr2=1.
    dXbin=(Xr2-Xr1)/((Xnb));
    titleRMSVsGen_ptTop_full ="; Truth;RMS"
    titleMeanVsGen_ptTop_full ="; Truth;Mean"
    h_RMSVsGen_=ROOT.TH1F()
    h_RMSVsGen_.SetDirectory(0)
    h_RMSVsGen_.SetBins(Xnb,Xr1,Xr2)
    h_RMSVsGen_.SetTitleOffset(2.0)
    h_RMSVsGen_.GetXaxis().SetTitleOffset(1.20)
    h_RMSVsGen_.GetYaxis().SetTitleOffset(1.30)
    h_RMSVsGen_.SetTitle(titleRMSVsGen_ptTop_full);
    h_RMSVsGen_.SetStats(0)
    h_meanVsGen_=ROOT.TH1F()
    h_meanVsGen_.SetDirectory(0)
    h_meanVsGen_.SetBins(Xnb,Xr1,Xr2)
    h_meanVsGen_.SetTitleOffset(2.0)
    h_meanVsGen_.GetXaxis().SetTitleOffset(1.20)
    h_meanVsGen_.GetYaxis().SetTitleOffset(1.30)
    h_meanVsGen_.SetTitle(titleMeanVsGen_ptTop_full)
    h_meanVsGen_.SetStats(0)
    for i in range(Xnb):
        h_RMSVsGen_.SetBinContent(i+1,(histo.ProjectionY("_py",histo.GetXaxis().FindFixBin(Xr1+i*dXbin) ,histo.GetXaxis().FindFixBin(Xr1+(i+1)*dXbin),"")).GetRMS());
        h_RMSVsGen_.SetBinError(i+1,(histo.ProjectionY("_py",histo.GetXaxis().FindFixBin(Xr1+i*dXbin) ,histo.GetXaxis().FindFixBin(Xr1+(i+1)*dXbin),"")).GetRMSError());
        h_meanVsGen_.SetBinContent(i+1,(histo.ProjectionY("_py",histo.GetXaxis().FindFixBin(Xr1+i*dXbin) ,histo.GetXaxis().FindFixBin(Xr1+(i+1)*dXbin),"")).GetMean());
        h_meanVsGen_.SetBinError(i+1,(histo.ProjectionY("_py",histo.GetXaxis().FindFixBin(Xr1+i*dXbin) ,histo.GetXaxis().FindFixBin(Xr1+(i+1)*dXbin),"")).GetMeanError());
    h_RMSVsGen_.SetStats(0)
    h_RMSVsGen_.Draw()
    c.SaveAs(weight_dir+"rms.png")
    c.Clear()
    h_meanVsGen_.SetStats(0)
    h_meanVsGen_.Draw()
    c.SaveAs(weight_dir+"mean.png")
