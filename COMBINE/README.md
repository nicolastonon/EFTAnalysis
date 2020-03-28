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
-------------------------------------------->

#### Useful links

:link: `Combine` page (browse the *Running combine* tab, or use the *Search* function): http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/

:link: CMS Combine Hypernews : https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination.html

:link: CombineHarvester page (external library) : http://cms-analysis.github.io/CombineHarvester/

_____________________________________________________________________________
#### Table Of Contents

* [Setup](https://github.com/nicolastonon/EFT-Simu-Pheno/tree/master/myAnalysis/COMBINE#Setup)

* [Generate datacards](https://github.com/nicolastonon/EFT-Simu-Pheno/tree/master/myAnalysis/COMBINE#setup)
  * [Details](https://github.com/nicolastonon/EFT-Simu-Pheno/tree/master/myAnalysis/COMBINE#details)

* [Combine commands](https://github.com/nicolastonon/EFT-Simu-Pheno/tree/master/myAnalysis/COMBINE#combine-commands)
  * [Expected significance](https://github.com/nicolastonon/EFT-Simu-Pheno/tree/master/myAnalysis/COMBINE#expected-significance-profile-likelihood)
  * [Observed significance](https://github.com/nicolastonon/EFT-Simu-Pheno/tree/master/myAnalysis/COMBINE#observed-significance-profile-likelihood)
  * [Asymptotic frequentist limits](https://github.com/nicolastonon/EFT-Simu-Pheno/tree/master/myAnalysis/COMBINE#asymptotic-frequentist-limits)
  * [Postfit templates, nuisances, etc.](https://github.com/nicolastonon/EFT-Simu-Pheno/tree/master/myAnalysis/COMBINE#postfit-templates-nuisances-etc-mlf)
  * [Nuisance parameters](https://github.com/nicolastonon/EFT-Simu-Pheno/tree/master/myAnalysis/COMBINE#nuisance-parameters)


_____________________________________________________________________________


# Setup

## Install Combine

`Combine` is a top-level CMSSW Package (i.e. it should be located at `$CMSSW_BASE/src/HiggsAnalysis/CombinedLimit`)
Follow the official [installation instructions](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#cc7-release-cmssw_10_2_x-recommended-version).

## Install CombineHarvester

`CombineHarvester` is a top-level CMSSW Package (i.e. it should be located at `$CMSSW_BASE/src/CombineHarvester`) which contains some additional functions and scripts (e.g. to plot systematics impacts or access postfit distributions).
The CMSSW version that should be used with `CombineHarvester` is driven by the recommendation for the HiggsAnalysis/CombinedLimit package, which is also required.

To install it, follow the instructions from the [CombineHarvester main page](http://cms-analysis.github.io/CombineHarvester/) or the [CombineHarvester GitHub page](https://github.com/cms-analysis/CombineHarvester) :

```
cd CMSSW_$RELEASE/src

#git clone via ssh:

bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/master/CombineTools/scripts/sparse-checkout-ssh.sh)

#git clone via https:

bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/master/CombineTools/scripts/sparse-checkout-https.sh)
```

Make sure to recompile the code :
```
cd $CMSSW_BASE/src/
scramv1 b clean
scramv1 b
```


# Generate datacards

Move to the [datacards](https://github.com/nicolastonon/EFT-Simu-Pheno/tree/master/myAnalysis/COMBINE/datacards) directory to generate the final datacards.

*TL;DR* :

```
make
./Generate_Datacard_Template.exe
./ScriptProducer_GenerateAllDatacards.exe
```

## Details

- [Generate_Datacard_Template.cxx](Generate_Datacard_Template.cxx) : creates the template to be used to produce all the final datacards.
Modify directly the *main()* to adapt it to your in needs.
In particular, you can set there the list of processes, of systematic uncertainties (process normalization, normalization-only systematic, shape systematic).

  Running this code will generate the template datacard 'Template_Datacard.txt' :

  ```
  make
  ./Generate_Datacard_Template.exe
  ```

- [Parser_Datacard_Template.py](Parser_Datacard_Template.py)* : parses the template datacard (see above), and replaces specific key-words with the desired values (channel name, systematic value, etc.)
This code is called when running *ScriptProducer_GenerateAllDatacards.exe*.

- [ScriptProducer_GenerateAllDatacards.cxx](ScriptProducer_GenerateAllDatacards.cxx) : asks the user to set options via command-line, and takes care of generating the final combined datacard.

  ```
  make
  ./ScriptProducer_GenerateAllDatacards.exe
  ```


# Combine commands

:information_source: *Make sure you're using the right datacard name in the commands ; `COMBINED_datacard.txt` is only a dummy name*

:information_source: *Make sure you have placed the relevant template file in the [templates](https://github.com/nicolastonon/EFT-Simu-Pheno/tree/master/myAnalysis/COMBINE/templates) directory, following the correct naming convention.*

:information_source: *If you will use `Combine` a lot, you'll find it hepful to set-up a few aliases for these commands...*

- To derive limits and significance, we use the **Asymptotic likelihood** methods :
  - `AsymptoticLimits`: limits calculated according to the asymptotic formulas in [arxiv:1007.1727](https://arxiv.org/abs/1007.1727)
  - `Significance`: simple profile likelihood approximation, for calculating significances.

- To access postfit informations (shapes, uncertainties, etc.), we also use the **Fitting** methods :
  - `FitDiagnostics` [[doc](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#fitting-diagnostics)] : performs maximum likelihood fits to extract the signal yield and provide diagnostic tools such as pre and post-fit models and correlations
  - `MultiDimFit` [[doc](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/commonstatsmethods/#likelihood-fits-and-scans)] :  perform multi-dimensional fits and likelihood based scans/contours using models with several parameters of interest.



## Expected significance (Profile Likelihood)

The expected significance can be computed from an Asimov dataset of signal+background. There are two options for this :

- *a-posteriori* expected: will depend on the observed dataset.
- *a-priori* expected (the default behavior): does not depend on the observed dataset, and so is a good metric for optimizing an analysis when still blinded.

The a-priori expected significance from the Asimov dataset is calculated as :
```
combine -M Significance -t -1 --expectSignal=1 COMBINED_datacard.txt
```

:information_source: `-t -1` <-> Use Asimov Dataset. To use toys instead, use `-t N`, with N the number of toys

In order to produced the a-posteriori expected significance, just generate a post-fit Asimov (i.e add the option `--toysFreq` in the command above).
The output format is the same as for observed signifiances: the variable limit in the tree will be filled with the significance (or with the p-value if you put also the option `--pvalue`)


## Observed significance (Profile Likelihood)
The observed significance is calculated using the `Significance` method, as :

```
combine -M Significance COMBINED_datacard.txt
```

In case the fit has difficulties to converge, consider adding the following options :

```
combine -M Significance COMBINED_datacard.txt --rMin=0 --rMax=20 --X-rtd ADDNLL_RECURSIVE=0 --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerTolerance 0.01 --cminPreScan
```

:information_source: *Use option `-v X` (with X>=2) for verbose mode. With X=4, can access postfit POI values & get warning/error messages*

:information_source: *The output root file will contain the significance value in the branch limit. To store the p-value instead, include the option --pval*


## Asymptotic Frequentist Limits

### Un-blinded limits

:cl: The `AsymptoticLimits` method allows to compute quickly an estimate of the observed and expected limits, which is fairly accurate when the event yields are not too small and the systematic uncertainties don't play a major role in the result. The limit calculation relies on an asymptotic approximation of the distributions of the LHC test-statistic, which is based on a profile likelihood ratio, under signal and background hypotheses to compute CLs+b, CLb and therefore CLs=CLs+b/CLb - i.e it is the asymptotic approximation of computing limits with frequentist toys. This method is so commonly used that it is the default method.

The method can be run using :
```
combine -M AsymptoticLimits COMBINED_datacard.txt
```
The program will print out the limit on the signal strength r (number of signal events / number of expected signal events) e .g. `Observed Limit: r < 1.6297 @ 95% CL` , the median expected limit `Expected 50.0%: r < 2.3111` and edges of the 68% and 95% ranges for the expected limits.

:information_source: *By default, the limits are calculated using the CLs prescription, as noted in the output, which takes the ratio of p-values under the signal plus background and background only hypothesis. This can be altered to using the strict p-value by using the option `--rule CLsplusb` You can also change the confidence level (default is 95%) to 90% using the option `--cl 0.9` or any other confidence level. You can find the full list of options for `AsymptoticLimits` using `--help -M AsymptoticLimits`*



### Blinded limits

The `AsymptoticLimits` calculation follows the frequentist paradigm for calculating expected limits. This means that the routine will first fit the observed data, conditionally for a fixed value of r and set the nuisance parameters (NPs) to the values obtained in the fit for generating the Asimov data, i.e it calculates the post-fit or a-posteriori expected limit. In order to use the pre-fit nuisance parameters (to calculate an a-priori limit), you must add the option --noFitAsimov or `--bypassFrequentistFit`.

For blinding the results completely (i.e not using the data) you must include the option `--run blind` :
```
combine -M AsymptoticLimits --noFitAsimov --run blind COMBINED_datacard.txt
```

:no_entry_sign: From the instructions :
> **You should never use -t -1 to get blind limits!**

*(Otherwise, the "nominal" nuisance parameter values are taken from fits to the data and are therefore not "blind" to the observed data by default (following the fully frequentist paradigm) ).*



## Postfit templates, nuisances, etc. (MLF)

:heavy_exclamation_mark: You need to create output dir. first, to avoid segfault.

For details on recipes and outputs, it is advised to refer to the documentation page ([FitDiagnostics](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#plotting)).

Run the following command :
```
combine -M FitDiagnostics -t -1 --expectSignal 1 --saveShapes --saveWithUncertainties --plots --out outputs COMBINED_datacard.txt #expected
```
or if you want to unblind the data :
```
combine -M FitDiagnostics --saveShapes --saveWithUncertainties --plots --out outputs COMBINED_datacard.txt #observed
```
This will produce plots and a `fitDiagnostics.root` file in the specified output directory. This file contains many informations which can be later used to produce postfit plots, retrieve postfit uncertainties, postfit NPs, etc.

- `--plots` : will produce pre/post-fit plots in the output directory (however these plots are not really satisfying : backgrounds are not split, etc., so you probably want to produce them by yourself).

- `--saveShapes` : will produce additional folders in `fitDiagnostics.root` containing the pre/post-fit distributions for each category/process. The folders are named `shapes_prefit`/`shapes_fit_sb`/`shapes_fit_b` (corresponding to prefit/fit with null signal/fit with floating signal). The distributions are :
  - `data` : TGraphAsymmErrors containing the observed data (or toy data if using `-t`). The vertical error bars correspond to the 68% interval for a Poisson distribution centered on the observed count.
  - `$PROCESS` : TH1F for each process in channel, named as in the datacard.
  - `total_signal`, `total_background`, `total`.

  The above distributions are provided for each channel included in the datacard, in separate sub-folders, named as in the datacard: there will be one sub-folder per channel.

  :information_source: *The pre-fit signal is by default for r=1 but this can be modified using the option --prefitValue.*

  :information_source: *The histograms/data points will have the bin number as the x-axis and the content of each bin will be a number of events.*

- `--saveWithUncertainties` : add uncertainties on the shapes. These uncertainties are generated by re-sampling of the fit covariance matrix, thereby accounting for the full correlation between the parameters of the fit.

  :information_source: *It may be tempting to sum up the uncertainties in each bin (in quadrature) to get the total uncertainty on a process however, this is (usually) incorrect as doing so would not account for correlations between the bins. Instead you can refer to the uncertainties which will be added to the post-fit normalizations described above.*

  :information_source: *This [HN answer](https://hypernews.cern.ch/HyperNews/CMS/get/higgs-combination/991/1.html) further explains how the prefit error bars are obtained : "random sets of nuisance parameters are generated from their respective pdfs and the shapes are reevaluated each time. The RMS of the deviation with respect to the nominal bin content is used to define the bin error."*


## Nuisance parameters

### Pre/postfit NPs and pulls [doc](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#pre-and-post-fit-nuisance-parameters-and-pulls)
*It is possible to compare pre-fit and post-fit nuisance parameters with the script `diffNuisances.py` (copied locally in the Plotting/ dir. for conveninence).
Taking as input a `fitDiagnostics.root` file, the script will by default print out the parameters which have changed significantly w.r.t. their initial estimate. For each of those parameters, it will print out the shift in value and the post-fit uncertainty, both normalized to the input values, and the linear correlation between the parameter and the signal strength.
The output by default will be the changes in the nuisance parameter values and uncertainties, relative to their initial (pre-fit) values (usually relative to initial values of 0 and 1 for most nuisance types).*

- Run the `Combine` script from the `Plotting/` dir. :
  ```
  cd Plotting
  python diffNuisances.py -g nuisances.root ../datacards/outputs/fitDiagnostics.root  --all --abs
  ```
  :arrow_right: This will create a `nuisances.root` file containing prefit/postfit infos on the NPs.

- You can then produce a plot out of it with :
  ```
  root -l Draw_Canvas.C
  ```

  :arrow_right: This creates 2 png plots.


### Impact of NPs [[doc](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#nuisance-parameter-impacts)]

*It is possible to estimate the impact of the nuisance parameters.
This is effectively a measure of the correlation between the NP and the POI, and is useful for determining which NPs have the largest effect on the POI uncertainty.*

:information_source: This can be quite long and un-readable if all the NPs are included.
You might want to remove the `autoMCStats` lines from the datacard in order to ignore the statistical uncert. for all bins...

:information_source: You can rename parameters by editing rename.json, for example :
```
{
  "r" : "#mu",
  "lumi" : "luminosity"
}
```


- Run the bash script :
  ```
  cd Plotting
  ./Plot_Syst_Impacts.sh [DATACARD_NAME_NOEXT] [EXP_OBS]
  ```
  :information_source: The 1rst argument is the name of the datacard without its file extension. The 2nd argument ie either "exp" or "obs", depending if you want to blind the data or not.
  :arrow_right: This creates a pdf plot representing the pulls of the NPs.


<!--
# OBSOLETE (needs update)

:heavy_exclamation_mark: **Here are outdated commands, for the record only**

### Breakdown of systematics

- To do the break-down of systematics : [Instructions](https://indico.cern.ch/event/577649/contributions/2388797/attachments/1380376/2098158/HComb-Tutorial-Nov16-Impacts.pdf)

- To estimate stat. uncertainty contribution : [Instructions](https://indico.cern.ch/event/577649/contributions/2388797/attachments/1380376/2098158/HComb-Tutorial-Nov16-Impacts.pdf)

- Used to be these commands (**needs update**) :
  ```
  text2workspace.py datacard.txt
  combine -M MultiDimFit --algo grid --points 50 --rMin -1 --rMax 4  datacard.root -m 125 -n nominal --expectSignal=1 -t -1
  combine -M MultiDimFit --algo none --rMin -1 --rMax 4  datacard.root -m 125 -n bestfit  --saveWorkspace --expectSignal=1 -t -1
  combine -M MultiDimFit --algo grid --points 50 --rMin -1 --rMax 4 -m 125 -n stat higgsCombinebestfit.MultiDimFit.mH125.root --snapshotName MultiDimFit --freezeParameters all --expectSignal=1 -t -1
  plot1DScan.py --others 'higgsCombinestat.MultiDimFit.mH125.root:Freeze all:2' --breakdown syst,stat higgsCombinenominal.MultiDimFit.mH125.root
  ```
-------------------------------------------->
