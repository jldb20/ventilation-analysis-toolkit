# Ventilation Analysis Toolkit
Copyright Jonathan du Bois, University of Bath, 2020

## Contents

### main functions
file                  | description  
---                   | ---  
importVT.m            | imports raw data from Fluke VT+ .sig and .par files  
importWrapper.m       | calls importVT and does some preprocessing on the data so it's easier to use  
splitCycle.m          | segments signal data into complete cycles, and splits inspiration and expiration phases. Also produces nice plots of the data, and some summary stats for the signals.  
identifyRestriction.m | tries to identify a restriction between two flow meters, either for one test case or for several datasets amalgamated. This uses metohd 1, trying to align the two signals intheir entirety. The time signal seems to drift between the two though, so a second method is being investigated.  
doCorrelate.m         | correlates two signals and returns aligned versions of the signals. Used for method 1 in the restriction identification.  

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
./FlowMeter1/           | directory with data files from ventilator component tests of 2020-03-29 (prefixes 'M' and 'S') and 2020-04-06 (prefix 'v').  
./FlowMeter2/           |directory with data files from ventilator component tests of 2020-03-29 (prefixes 'M' and 'S') and 2020-04-06 (prefix 'v').  

### reference material
(under ./reference_material/ directory)  
file                  | description  
---                   | ---  
2020-04-09_restrictor_diagrams | initial plots of restrictor coeffiecient estiamtes. Includes some useful test case numbers and descriptions.

## License
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
    
