#Specific validation functions for regressor NN

import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from Utils.ColoredPrintout import colors
from Utils.Helper import close_event

# //--------------------------------------------
# //--------------------------------------------

def Plot_LR_Pred_vs_Truth(opts, list_features, list_labels, list_yTrain_allClasses, list_yTest_allClasses, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, list_truth_Test_allClasses, list_xTrain_allClasses, list_xTest_allClasses, weight_dir):
    """
    Make validation plots for regressor NN using test data. Compare predictions to true target values.
    """

    if opts["strategy"] not in ["ROLR", "RASCAL"]: return

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
    pred_data = np.squeeze(np.concatenate((list_predictions_test_allNodes_allClasses[0][0],list_predictions_test_allNodes_allClasses[0][1])))
    class_data = np.squeeze(np.concatenate((list_truth_Test_allClasses[0],list_truth_Test_allClasses[1])))
    # print(truth_data.shape); print(pred_data.shape); print(class_data.shape)

    quantile = np.quantile(pred_data, 0.995)
    xmax=abs(quantile)
    ymax=abs(quantile)

# //--------------------------------------------
#Scatterplot
#See https://seaborn.ppred_data.org/tutorial/distributions.html

    fig = plt.figure('splot', figsize=(10, 10))
    # timer = fig.canvas.new_timer(interval = 1000) #creating a timer object and setting an interval of N milliseconds
    # timer.add_callback(close_event)
    plt.title('Predicted VS True '+nodename)
    plt.xlabel(r'True '+nodename+r'(x|$\theta_0,\theta_1$)', fontsize=15) # add 'r' in front <-> interpreted as raw string
    plt.ylabel(r'Learned '+nodename+r'(x|$\theta_0,\theta_1$)', fontsize=15) # add 'r' in front <-> interpreted as raw string #color='darkorange'

    splot = sns.scatterplot(x=truth_data, y=pred_data, hue=class_data, style=class_data)
    leg_handles = splot.get_legend_handles_labels()[0]
    splot.legend(leg_handles, ['EFT', 'SM']) #NB: EFT first because class label is 0
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
    jplot
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

# //--------------------------------------------

    return
