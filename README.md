CMSSW_8_0_8
Description on Wiki for AlignmentTool
https://github.com/chiyi/ESAlignmentWorks/wiki/AlignmentTool

Provide DB files

1. Checkout package:
	```
	cmsrel CMSSW_8_0_8
	cd CMSSW_8_0_8/src
	cmsenv
	git cms-addpkg Geometry/EcalAlgo
	git cms-addpkg Geometry/CaloEventSetup
	git clone https://github.com/alphatsai/AlignmentTool.git
	scram b -j16
	```

2. Produce DB files which contain ES coordinate:
	```
	cd Geometry/CaloEventSetup/test
	vi TestWriteESAlignments.cc
	```

3. Add the aligned coordinate which are with respect to ideal coordinate:

	Example in Run1: https://hypernews.cern.ch/HyperNews/CMS/get/ecal-calibration/560.html
	Example under ESAlignTool/genDBExample/TestWriteESAlignments.cc
	
