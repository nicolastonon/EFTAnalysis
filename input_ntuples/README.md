After getting new NTuples, always run the [Split_FullSamples.cxx] code in order to:
1) Split full samples into sub-samples per sub-regions (so that then you only process the relevant subsets of events)
2) Create the [NPL.root] sample from the [DATA.root] sample (events satisfying isFake==1)
3) Create the [NPL_MC.root] sample from all MC samples (events satisfying (isPromptMC==1 && isFake==1))
