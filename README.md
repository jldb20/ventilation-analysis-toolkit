# Ventilation Analysis Toolkit

Data and MATLAB code for analysing ventilator circuits.  
These contents have grown out of the need for tools to analyse experiments on Ventilator Circuts conducted at the Royal United Hospitals in Bath with the University of Bath, UK, in 2020. In this repository are the raw datasets and the tools used to analyse them. These tools have evolved rapidly and are not always the tidiest, most robust, or most efficient code - help improving them is welcome!

## Quick start examples

If you want to see the identification results for some or all of the 'v'-prefixed restrictor tests, then:
1. Run one of:
    restrictorSweep('v',{112:114},true,[],0.2); % select results
    restrictorSweep('v',{1:3,4:6,7:9,10:12,13:15,16:18,19:21,22:24,25:27,28:30,31:33,34:36,37:39,40:42,109:111,112:114,115:117},true,[],0.2); % all results
2. Paste tab-delimited clipboard data directly into Excel or a text editor. (On *nix you can copy output from the concole to paste; MATLAB clipboard functions may not work?)
3. Look at figures in ./reference_material/restrictor_graphs to see results.

For more detail on the processing of a specific set of tests, try:
    identifyRestriction2({'v112','v113','v114'},0.2,[],'standard')

And to understand the preprocessing of one of these sets try:
    sig=importWrapper('v112'); splitCycle(sig,[],0.2)
or for a bidirectional flow check out
    sig=importWrapper('v43'); splitCycle(sig,[],0.2)

If you just want raw data, have a look at importWrapper, or importVT e.g.
    [sig,par] = importWrapper('v43');
          sig = importVT('FlowMeter1/v43.sig']);


## Contents

### main functions
file                   | description  
---                    | ---  
importVT.m             | imports raw data from Fluke VT+ .sig and .par files  
importWrapper.m        | calls importVT and does some preprocessing on the data so it's easier to use  
splitCycle.m           | segments signal data into complete cycles, and splits inspiration and expiration phases. Also produces nice plots of the data, and some summary stats for the signals.  
identifyRestriction.m  | tries to identify a restriction between two flow meters, either for one test case or for several datasets amalgamated. This uses metohd 1, trying to align the two signals intheir entirety. The time signal seems to drift between the two though, so a second method is being investigated.  
identifyRestriction2.m | This is now the preferred method for identifying restrictions. It uses the pressure difference from phase-averaged signals from two meters.
doCorrelate.m          | correlates two signals and returns aligned versions of the signals. Used for method 1 in the restriction identification.  
restrictorSweep.m      | sweeps through multiple test cases (or groups of tests) and compiles outputs in tab-delimited form. Can also save graphs showing identification.

### supporting functions
file                  | description  
---                   | ---  
matchXaxes.m          | keeps x-axes aligned on plots in the same figure window when zooming, panning etc.  
linkedPoints.m        | allows clicking on points in one graph to highlight respective points in another  
fitFunc.m             | does the curve fitting for identifyRestriction().  
plotRestriction.m     | does the plotting for identifyRestriction().  

### datasets
file                  | description  
---                   | ---  
./FlowMeter1/         | directory with data files from ventilator component tests of 2020-03-29 (prefixes 'M' and 'S') and 2020-04-06 (prefix 'v').  
./FlowMeter2/         |directory with data files from ventilator component tests of 2020-03-29 (prefixes 'M' and 'S') and 2020-04-06 (prefix 'v').  

### reference material
(under ./reference_material/ directory)  
file                           | description  
---                            | ---  
2020-04-09_restrictor_diagrams | initial plots of restrictor coeffiecient estiamtes. Includes some useful test case numbers and descriptions.
restrictor_graphs/             | folder where output of restrictorSweep() is saved - shows identification curves for restrictors, and identified coefficients.

## License
Copyright Jonathan du Bois, University of Bath, 2020

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
