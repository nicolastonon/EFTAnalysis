#TBD : standalone code to load DNN and test it as desired (independant from main code)

import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from Utils.ColoredPrintout import colors

# //--------------------------------------------
# //--------------------------------------------

def Plot_LR_Pred_vs_Truth(opts, list_features, list_labels, list_predictions_train_allNodes_allClasses, list_predictions_test_allNodes_allClasses, list_PhysicalWeightsTrain_allClasses, PhysicalWeights_allClasses, list_PhysicalWeightsTest_allClasses, y_train, y_test, y_process, y_process_train, y_process_test, predictions_train, predictions_test, weight_dir):

    if opts["strategy"] not in ["ROLR", "RASCAL"]: return

    # Plotting the loss with the number of iterations
    fig = plt.figure('LR', figsize=(10, 10))
    ax = fig.gca()

    plt.title('')
    plt.ylabel('Learned log r(x|$\theta_0,\theta_1$)', color='darkorange')
    plt.xlabel('True log r(x|$\theta_0,\theta_1$)')

    sns.scatterplot(x=y_test, y=np.squeeze(predictions_test))
    # sns.regplot(x=y_test, y=predictions_test)
    # ax.set(xlim=(-10, 10))
    # ax.set(ylim=(-10, 10))
    ax.set(xlim=(-0.5, 3))
    ax.set(ylim=(-0.5, 3))

    # Plot your initial diagonal line based on the starting
    # xlims and ylims.
    diag_line, = ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c=".3")

    plotname = weight_dir + 'LR_PredvsTruth.png'
    fig.savefig(plotname)
    print(colors.fg.lightgrey, "\nSaved LR correlation plot as :", colors.reset, plotname)
    fig.clear()
    plt.close('LR')

    return
