After getting new NTuples, run the [Split_FullSamples.cxx] code in order to:
1) Create 'NPL_DATA' (from data) and 'NPL_MC' (from MC) ntuples
2) Store the EFT parameterizations in the private SMEFT samples if desired (faster reading)
3) ...

*Obsolete* Ideally, would also split full ntuples by sub-categories to then speed up the reading. But this consumes too much disk space, can't be done locally.
