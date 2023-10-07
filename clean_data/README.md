# README 

Data used for the Bimler 2023 paper 'Plant interaction networks reveal the limits of our understanding of diversity maintenance', see methods in the main text and Supplementary Methods S1.1 for details on the data collection and experimental design.

Data primarily consists of individual observations of plant fecundity and the identity and abundance of surrounding herbaceous plant neighbours, collected from a natural, annual wildflower community in Western Australia. These data were binned into three environmental categories based on the percentage of overhead canopy cover: 'open' (0 to 7.9% canopy cover), 'intermediate' (8 to 17.9% canopy cover), and 'shady' (18 to 40% canopy cover), each with their own separate data files. These categories are referred to in the names of the data files below as 1, 2, and 3, respectively.


## File structure

### Community data split into environmental categories:

The files below exist for each environmental category 1, 2, and 3, as given at the end of the file names. We describe each type of file only once to avoid repetition. 

**fecundities*.csv:** rows are observations of individual seed production and neighbour abundances. The first four columns refer to: 
plot = plot number (identifier) from 1 to 100
seeds = observation of indivdual seed production
focal = four-letter focal species code (identifier)
focalID = individual plant identifier 
Each column following those four is named after a species code and contains the abundance of each neighbouring species for that individual observation. Neighbour columns are organised alphabetically, with focal species first then non-focal neighbour species. This should match the size and order of the associated key_neighbourID*.csv file. 

**key_speciesID*.csv:** single-column list of four-letter species codes identifying focal species in that environmental category, in alphabetical order.

**key_neighbourID*.csv:** single-column list of four-letter species codes identifying all neighbouring species in that environmental category. Focal species comes first in alphabetical order, followed by non-focal neighbours also in alphabetical order. 

### Plot and demographic data:

**plot_data.csv:** Environmental, diversity and density data associated with each plot. The column names refer to: 
Plot = plot number (identifier) from 0 to 100
Density = whether the plot was left at full density (100) or thinned to 60% (60) or 30% (30) density without targeting any particular species.
Per.Litt = percentage litter cover
Per.Bare = percentage bare ground
Per.Woody = percentage woody debris (logs, fallen branches) 
Soil.P = soil phosphorus as measured by Colwell extraction
Soil.Water = soil moisture content
Canopy.Cover = percentage canopy cover 
Class.Canopy = environmental category whereby 1='open', 2='intermediate' and 3='shady'
evenness = plot-level Shannon evenness index
diversity = total plot-level species richness 
total.density = total plot-level abundance of all species

**plot_species_abundances.csv:** rows are abundance of each species present in each plot. The columns refer to:
plot = plot number (identifier) from 0 to 100
species = four-letter species code
Num_indivs = number of individuals of that species in that plot
Num_indiv_seed = number of individuals for which seeds were collected
Class.Canopy = environmental category whereby 1='open', 2='intermediate' and 3='shady'

**seed_rates.csv:** estimates of focal species seed germination and seed survival rates. The community mean is substituted for those species which did not have reliable estimates (STPA, PEDU, WAAC, EROD, GITE survival rate).
The first column is intended to be used as row names when read into R. 
name = focal species name
code = four-letter species code
germ = seed germination rate from 0 to 1
surv = seed survival rate from 0 to 1

**species_list.csv:** list of species present in the dataset. Columns give, in order: four-letter species code, species name, family, whether they are considered focal or not, and size of the interaction neighbourhood

