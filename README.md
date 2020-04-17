# Ventilation Analysis Toolkit
Copyright Jon du Bois, 2020

## Contents

### main functions
importVT.m imports raw data from Fluke VT+ .sig and .par files
importWrapper.m calls importVT and does some preprocessing on the data so it's easier to use
splitCycle.m segments signal data into complete cycles, and splits inspiration and expiration phases. Also produces nice plots of the data.

### supporting functions
matchXaxes.m keeps x-axes aligned on plots in the same figure window when zooming, panning etc.

### datasets
FlowMeter1/ directory with data files from ventilator component tests of 2020-03-29 (prefixes 'M' and 'S') and 2020-04-06 (prefix 'v').
FlowMeter2/ directory with data files from ventilator component tests of 2020-03-29 (prefixes 'M' and 'S') and 2020-04-06 (prefix 'v').


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
    
