#TBD : standalone code to load DNN and test it as desired (independant from main code)

import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from Utils.ColoredPrintout import colors
from Utils.Helper import close_event

# //--------------------------------------------
# //--------------------------------------------

#Other ideas:
#- https://seaborn.pydata.org/tutorial/distributions.html
#- https://stackoverflow.com/questions/20105364/how-can-i-make-a-scatter-plot-colored-by-density-in-matplotlib
def Plot_LR_Pred_vs_Truth(opts, list_features, list_labels, list_yTrain_allClasses, list_yTest_allClasses, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_truth_Test_allClasses, weight_dir):

    if opts["strategy"] not in ["ROLR", "RASCAL"]: return

    nodename='r'
    if opts['regress_onLogr'] == True: nodename='log(r)'

    xmin=-3, xmax=+3
    ymin=-3, ymax=+3

    mycol = 'g'

# //--------------------------------------------

    fig = plt.figure('splot', figsize=(10, 10))
    timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
    timer.add_callback(close_event)
    plt.title('Predicted VS True '+nodename)
    plt.ylabel(r'Learned '+nodename+'(x|$\theta_0,\theta_1$)', fontsize=15) # add 'r' in front <-> interpreted as raw string #color='darkorange'
    plt.xlabel(r'True '+nodename+'(x|$\theta_0,\theta_1$)', fontsize=15) # add 'r' in front <-> interpreted as raw string

    if opts["nofOutputNodes"] == 1: splot = sns.scatterplot(x=np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1])), y=np.concatenate((np.squeeze(list_predictions_test_allNodes_allClasses[0][0]),np.squeeze(list_predictions_test_allNodes_allClasses[0][1]))), hue=np.concatenate((list_truth_Test_allClasses[0],list_truth_Test_allClasses[1])), style=np.concatenate((list_truth_Test_allClasses[0],list_truth_Test_allClasses[1])) )
    else: splot = sns.scatterplot(x=np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1]))[:,0], y=np.concatenate((list_predictions_test_allNodes_allClasses[0][0],list_predictions_test_allNodes_allClasses[0][1])), hue=np.concatenate((list_truth_Test_allClasses[0],list_truth_Test_allClasses[1])), style=np.concatenate((list_truth_Test_allClasses[0],list_truth_Test_allClasses[1])) )

    leg_handles = splot.get_legend_handles_labels()[0]
    splot.legend(leg_handles, ['SM', 'EFT']) #title='New legend'

    ax = fig.gca()
    ax.set(xlim=(xmin, xmax))
    ax.set(ylim=(ymin, ymax))

    # Plot your initial diagonal line based on the starting
    # xlims and ylims.
    diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")

    timer.start()
    plt.show()

    plotname = weight_dir + 'ScatterPlotLR_PredvsTruth.png'
    fig.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR scatter plot as :", colors.reset, plotname)
    fig.clear(); plt.close('splot')

# //--------------------------------------------

    fig = plt.figure('regplot', figsize=(10, 10))
    ax = fig.gca(); ax.set(xlim=(xmin, xmax)); ax.set(ylim=(xmin, xmax))
    if opts["nofOutputNodes"] == 1: regplot = sns.regplot(x=np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1])), y=np.concatenate((np.squeeze(list_predictions_test_allNodes_allClasses[0][0]),np.squeeze(list_predictions_test_allNodes_allClasses[0][1]))), color=mycol)
    else: regplot = sns.regplot(x=np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1]))[:,0], y=np.concatenate((list_predictions_test_allNodes_allClasses[0][0],list_predictions_test_allNodes_allClasses[0][1])), color=mycol)
    plotname = weight_dir + 'RegPlot_PredvsTruth.png'
    fig.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR reg plot as :", colors.reset, plotname)
    fig.clear(); plt.close('regplot')

# //--------------------------------------------

    #NB: 'jointplot' func creates its own figure (as do: FacetGrid, factorplot, lmplot, PairGrid, pairplot, JointGrid, jointplot.)
    if opts["nofOutputNodes"] == 1: jplot = sns.jointplot(x=np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1])), y=np.concatenate((np.squeeze(list_predictions_test_allNodes_allClasses[0][0]),np.squeeze(list_predictions_test_allNodes_allClasses[0][1]))), color=mycol, alpha=0.5)
    else: jplot = sns.jointplot(x=np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1]))[:,0], y=np.concatenate((list_predictions_test_allNodes_allClasses[0][0],list_predictions_test_allNodes_allClasses[0][1])), color=mycol, alpha=0.5)
    # g = sns.JointGrid(x, y)
    # g.ax_marg_x.hist(x, bins=np.arange(0, 60))
    # g.ax_marg_y.hist(y, bins=np.arange(0, 1000, 10), orientation="horizontal")
    jplot.ax_marg_x.set_xlim(xmin, xmax); jplot.ax_marg_y.set_ylim(ymin, ymax)
    plotname = weight_dir + 'JointPlot_PredvsTruth.png'
    jplot.savefig(plotname)
    jplot
    print(colors.fg.lightgrey, "\nSaved LR joint plot as :", colors.reset, plotname)

# //--------------------------------------------

    if opts["nofOutputNodes"] == 1: hexplot = sns.jointplot(x=np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1])), y=np.concatenate((np.squeeze(list_predictions_test_allNodes_allClasses[0][0]),np.squeeze(list_predictions_test_allNodes_allClasses[0][1]))), color=mycol, alpha=0.5, kind="hex")
    else: hexplot = sns.jointplot(x=np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1]))[:,0], y=np.concatenate((list_predictions_test_allNodes_allClasses[0][0],list_predictions_test_allNodes_allClasses[0][1])), color=mycol, alpha=0.5, kind="hex")
    hexplot.ax_marg_x.set_xlim(xmin, xmax); jplot.ax_marg_y.set_ylim(ymin, ymax)
    plotname = weight_dir + 'JointPlotHEX_PredvsTruth.png'
    hexplot.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR joint HEX plot as :", colors.reset, plotname)

# //--------------------------------------------

    if opts["nofOutputNodes"] == 1: kdeplot = sns.jointplot(x=np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1])), y=np.concatenate((np.squeeze(list_predictions_test_allNodes_allClasses[0][0]),np.squeeze(list_predictions_test_allNodes_allClasses[0][1]))), color=mycol, alpha=0.5, kind="kde")
    else: kdeplot = sns.jointplot(x=np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1]))[:,0], y=np.concatenate((list_predictions_test_allNodes_allClasses[0][0],list_predictions_test_allNodes_allClasses[0][1])), color=mycol, alpha=0.5, kind="kde")
    kdeplot.ax_marg_x.set_xlim(xmin, xmax); jplot.ax_marg_y.set_ylim(ymin, ymax)
    plotname = weight_dir + 'JointPlotKDE_PredvsTruth.png'
    kdeplot.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR joint KDE plot as :", colors.reset, plotname)

# //--------------------------------------------

    if opts["nofOutputNodes"] == 1:
        kdesplot = sns.jointplot(x=np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1])), y=np.concatenate((np.squeeze(list_predictions_test_allNodes_allClasses[0][0]),np.squeeze(list_predictions_test_allNodes_allClasses[0][1]))), color=mycol, alpha=0.5, kind="kde")
        kdesplot.plot_joint(plt.scatter, c="b", s=30, linewidth=1, marker="+")
    else:
        kdesplot = sns.jointplot(x=np.concatenate((list_yTest_allClasses[0],list_yTest_allClasses[1]))[:,0], y=np.concatenate((list_predictions_test_allNodes_allClasses[0][0],list_predictions_test_allNodes_allClasses[0][1])), color=mycol, alpha=0.5, kind="kde")
        kdesplot.plot_joint(plt.scatter, c="b", s=30, linewidth=1, marker="+")
    kdesplot.ax_joint.collections[0].set_alpha(0)
    kdesplot.ax_marg_x.set_xlim(xmin, xmax); jplot.ax_marg_y.set_ylim(ymin, ymax)
    # kdesplot.set_axis_labels("$X$", "$Y$");
    plotname = weight_dir + 'JointPlotKDESPLOT_PredvsTruth.png'
    kdesplot.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR joint KDE+SPLOT plot as :", colors.reset, plotname)

# //--------------------------------------------

    return
