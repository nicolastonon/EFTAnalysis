<!--
```
CODE EXAMPLE
```

=== Emoji list (see https://gist.github.com/rxaviers/7360908)
:arrow_right:
:heavy_exclamation_mark:
:heavy_check_mark:
:link:
:white_check_mark:
:heavy_multiplication_x:
:x:
:negative_squared_cross_mark:
:bangbang:
:white_check_mark:
:copyright:
:clock430:
:no_entry:
:ok:
:arrow_right_hook:
:paperclip:
:open_file_folder:
:chart_with_upwards_trend:
:lock:
:hourglass:
:warning:
:construction:
:fr:
:one: :two: :hash:
:underage:
:put_litter_in_its_place:
:new:


#HOW TO HIDE CONTENTS (which can be viewed by cliking icon) :
<details>
<summary>[NameOfHiddenContent]:</summary>
[theHiddenContent]
</details>
-------------------------------------------->
:construction: **README UNDER CONSTRUCTION**


> Codes and instructions to :
> * train BDTs (TMVA) and DNNs (Keras)
> * create MVA templates
> * create control plots
> * perform analysis-specific studies
> * generate and run COMBINE datacards


#### Table Of Contents

* [Setup](https://github.com/nicolastonon/EFTAnalysis#Setup)

* [Analysis](https://github.com/nicolastonon/EFTAnalysis#Analysis)
  * [Configuration](https://github.com/nicolastonon/EFTAnalysis#Configuration)

* [Event yields](https://github.com/nicolastonon/EFTAnalysis#Event-yields)

* [ROCS](https://github.com/nicolastonon/EFTAnalysis#ROCS)


* [Combine](https://github.com/nicolastonon/EFTAnalysis#Combine)


<!-- ![Mining gold](https://images.deepai.org/converted-papers/1805.00013/x1.png) -->

*[Credit : 10.1103/PhysRevLett.121.111801]*

# Setup

This standalone Ntuple analysis code does not require CMSSW, and can be run locally offline.

However, running `Combine` requires a CMSSW environment *(e.g. CMSSW_10_2_20)*.

Instructions for setting up `Combine` properly are given in the dedicated [directory](https://github.com/nicolastonon/EFTAnalysis/COMBINE).

## Interfacing with CMSSW

Please follow the following instructions to install and interface the code with CMSSW. This has been tested on NAF servers under SL7.

```
cd PATH
cmsrel CMSSW_11_1_2 [or another release; but the code is meant to work with TF2, apparently only available for >=11_1_0]
cd CMSSW_11_1_2/src; cmsenv
#git-cms-addpkg PhysicsTools/TensorFlow #NEEDED?
#git cms-addpkg FWCore/Framework #NEEDED?
#git cms-addpkg FWCore/Utilities #NEEDED?
scram b -j 4 [very slow]
git clone https://github.com/nicolastonon/EFTAnalysis.git
cd EFTAnalysis
[Adapt Makefile.cmssw if needed; comment/uncomment relevant headers in Utils/TFModel.h]
cp Makefile.cmssw Makefile
[Add 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:DIRPATH' into your .bashrc/profile settings]
make
```

## Input ntuples

The directory [input_ntuples](https://github.com/nicolastonon/EFTAnalysis/input_ntuples) is intended to contain the ntuples to be analyzed, produced via the [PoTATo](https://gitlab.cern.ch/cms-desy-topv/potato-common/-/tree/master) framework.
Alternatively, you may hardcode the relevant path to the ntuple directory in 'Utils/Helper.cxx'.
Ntuples must be stored into sub-directories corresponding to their years of data-taking ('2016'/'2017'/'2018').

You may run the code `input_ntuples/Split_FullSamples.cxx` in order to further split the Ntuples by sub-categories (for faster access), to compute and store the per-event EFT parameterizations in private SMEFT samples, to produce the data-driven Fakes sample, etc.

# Analysis

The basic workflow is the following :

* The [Makefile](https://github.com/nicolastonon/EFTAnalysis/tree/master/Makefile) compiles all the codes. Use [Makefile.cmssw](https://github.com/nicolastonon/EFTAnalysis/tree/master/Makefile.cmssw) instead if you're working under CMSSW.

* The [TopEFT_analysis](https://github.com/nicolastonon/EFTAnalysis/tree/master/TopEFT_analysis.cxx) class contains all the main analysis functions (make templates, plots, studies, etc.)

* The [Helper](https://github.com/nicolastonon/EFTAnalysis/tree/master/Helper.cxx) code contains additional, helper functions.

* The [analysis_main.cxx](https://github.com/nicolastonon/EFTAnalysis/tree/master/analysis_main.cxx) code contains the *main()* which runs the function calls. There are also defined all the customisable user-options.

**The user should only modify the [analysis_main.cxx](https://github.com/nicolastonon/EFTAnalysis/tree/master/analysis_main.cxx) code, compile, and run it.**

## Configuration

Modify the [analysis_main.cxx](https://github.com/nicolastonon/EFTAnalysis/tree/master/analysis_main.cxx) code to configure the analysis.

The configuration interface relies heavily on the use of vectors, and is intended to be self-explanatory.

In particular, you can easily configure the :
* general options (your signal, etc.)
* list of data-taking years to consider
* skimming cuts to apply to the ntuples
* list of samples
* list of BDT input variables (if used)
* list of additional variables, not used in the MVA but only e.g. for control plots
* list of shape systematics to compute and store

Select the functions to call by setting the corresponding booleans accordingly.

## Running the code

Copy either 'Makefile.local' or 'Makefile.cmssw' to 'Makefile', depending on your case.

```
make
./analysis_main.exe

#Optional -- specify a year and region, e.g. :
./analysis_main.exe 2016 signal
```

:arrow_right: Output root files (containing templates, control histograms, TMVA control file, etc.) are stored in the [outputs](https://github.com/nicolastonon/EFTAnalysis/tree/master/outputs) directory.

:arrow_right: Output plots are stored in the [plots](https://github.com/nicolastonon/EFTAnalysis/tree/master/plots) directory.

# Event yields

The code [Utils/Yield_Table.cxx](https://github.com/nicolastonon/EFTAnalysis/tree/master/Utils/Yield_Table.cxx) opens the input ntuples, and automatically computes the event yields for all selected processes.
Modify the *main()* options to select the data-taking year, processes, etc.

```
make
./Yield_Table.cxx.exe
```

# ROCS

Move to the [ROCS](https://github.com/nicolastonon/EFTAnalysis/tree/master/ROCS) directory.

The code [Compare_ROC_curves.cxx](https://github.com/nicolastonon/EFTAnalysis/tree/master/ROCS/Compare_ROC_curves.cxx) will read histogram files to plot the corresponding ROC curves.
It makes it easy to superimpose and compare several ROCS.
*NB1: plotting functions are implemented in [ROC_Plotter.cxx](https://github.com/nicolastonon/EFTAnalysis/tree/master/ROCS/ROC_Plotter.cxx)*
*NB2: specific naming conventions must be enforced for this code to work.*

Modify the *main()* options to define the paths of the input rootfiles, etc.

```
./make_plot.sh
```

# Combine

See the dedicated [README](https://github.com/nicolastonon/EFTAnalysis/COMBINE).
