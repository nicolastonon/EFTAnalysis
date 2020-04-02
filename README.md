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


![Mining gold](https://images.deepai.org/converted-papers/1805.00013/x1.png)
*[Credit : 10.1103/PhysRevLett.121.111801]*

# Setup

The analysis of Ntuples does not require CMSSW, and can be run locally offline.

However, running `Combine` requires a CMSSW environment *(e.g. CMSSW_10_2_20)*.

Instructions for setting up `Combine` properly are given in the dedicated [directory](https://github.com/nicolastonon/EFTAnalysis/COMBINE).

## Input ntuples

The directory [input_ntuples](https://github.com/nicolastonon/EFTAnalysis/input_ntuples) must contain the ntuples to analyze, produced via the [PoTATo](https://gitlab.cern.ch/joknolle/potato) framework.
They must be organized into sub-directories corresponding to their years of data-taking (2016/2017/2018).

# Analysis

The basic workflow is the following :

* The [Makefile](https://github.com/nicolastonon/EFTAnalysis/tree/master/Makefile) compiles all the codes.

* The [TopEFT_analysis](https://github.com/nicolastonon/EFTAnalysis/tree/master/TopEFT_analysis.cxx) class contains all the main analysis functions.

* The [Helper](https://github.com/nicolastonon/EFTAnalysis/tree/master/Helper.cxx) code contains additional, helper functions.

* The [analysis_main.cxx](https://github.com/nicolastonon/EFTAnalysis/tree/master/analysis_main.cxx) code contains the *main()* which runs the function calls.

**The user should only modify the [analysis_main.cxx](https://github.com/nicolastonon/EFTAnalysis/tree/master/analysis_main.cxx) code, compile, and run it.**

## Configuration

Modify the [analysis_main.cxx](https://github.com/nicolastonon/EFTAnalysis/tree/master/analysis_main.cxx) code to configure the analysis.

The configuration interface relies heavily on the use of vectors, and is intended to be self-explanatory.

In particuler, you can easily configure the :
* general options (your signal, etc.)
* list of data-taking years to consider
* skimming cuts to apply to the ntuples
* list of leptonic sub-channels
* list of samples
* list of BDT input variables
* list of additional variables, not used in the MVA but only e.g. for control plots
* list of shape systematics to compute and store

Select the steps to perform by setting the corresponding booleans accordingly.

## Running the code

```
make
./analysis_main.exe

#Optional -- specify a year, e.g. :
./analysis_main.exe 2016
```

:arrow_right: Output root files (containing templates, control histograms, TMVA control file, ...) are stored in the [outputs](https://github.com/nicolastonon/EFTAnalysis/tree/master/outputs) directory.

:arrow_right: Output plots are stored in the [plots](https://github.com/nicolastonon/EFTAnalysis/tree/master/plots) directory.

# Event yields

The code [Produce_Cutflow.cxx](https://github.com/nicolastonon/EFTAnalysis/tree/master/Produce_Cutflow.cxx) will read the input ntuples, and automatically compute the event yields for all processes.
Modify the *main()* options to select the data-taking year, processes, etc.

```
make
./Produce_Cutflow.exe
```
# ROCS

Move to the [ROCS](https://github.com/nicolastonon/EFTAnalysis/tree/master/ROCS) directory.

The code [Compare_ROC_curves.cxx](https://github.com/nicolastonon/EFTAnalysis/tree/master/ROCS/Compare_ROC_curves.cxx) will read TMVA control files to plot the corresponding ROC curves.
It makes it easy to superimpose and compare several ROCS.

Modify the *main()* options to define the paths of the input rootfiles, etc.

```
./make_plot.sh
```

# Combine

See the dedicated [README](https://github.com/nicolastonon/EFTAnalysis/COMBINE).
