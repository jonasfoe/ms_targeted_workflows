
# MS Targeted Workflows

This is a repository of scripts enabling targeted LC-MS peptide detection.
Workflows include peptide stable isotope labelling, acquisition parameter optimisation (NCE, FAIMS CV, RF Lens), scheduling for PRM acquisition, spectral library (chromatogram library) generation, and peptide detection diagnostics (XIC plots, peptide sequence coverage, spectral angle overview).
Skyline MS Software is used for peptide peak picking and Quarto notebooks (.qmd) are used for interactive evaluation of accessory scripts.

## Requirements

The following software is expected to be installed:
* R ≥ v4.2.0 (https://cran.r-project.org/)
* julia ≥ v1.9.0 (https://julialang.org)
* Quarto CLI ≥ v1.3.427 (https://quarto.org/)
* Skyline MS ≥ v22.2 (https://skyline.ms/)

R and julia projects need to be initilized from the root folder (can be done with 'debug_julia_and_R.cmd').
RStudio is probably the best option for running the .qmd notebooks.

## Modules

### !Peptide Lists

Keeps lists of target peptides for mixing and matching in the assays.
Optimization data is collected for peptides, such as optimal collision energy and relative LC retention time.

### !OE Assays

Helps set up inclusion lists / schedules for targeting peptides on orbitrap instruments.
Creates scan-by-scan tables for direct infusion MS, inclusion lists for targeted DDA experiments and scheduling for PRM exeriments

### !DI-MS Analyzers

Generate optimisation results (.pdf plots and .csv with optimisations) from direct infusion experiments.

### !Reports

A collection of report files that can be used to generate detection diagnostics from skyline files.
To use to .cmd files, use one folder per skyline file and call it id_all.sky

### !Skyline Templates

Template files for skyline. .sky template files can be used for evaluation of .raw files. .skyr are templates for the reports that are used to export skyline results.

## Workflows

### Prepare stable isotope labelling

Edit '!Peptide Lists\\!prepare_peptide_order_peptides.txt' to add peptides and mark heavy labels by using lowercase letters.
Then iteratively use '!Peptide Lists\\!prepare_peptide_order.qmd' to get feedback.
If the peptides end up in split aliquots, the splits must be submitted to new \*.txt peptide lists to set the split in stone.

### Add new peptides / new peptide lists

Create a new .txt file in '!Peptide Lists\\' of any name, containing the (new) peptides.
See '!Peptide Lists\\prtc_jptrt_26.txt' for an example.
Run '!Peptide Lists\\!make_precursor_lists.qmd' and '!Peptide Lists\\!make_optimizations_oe.qmd' successively.
This places the complete peptide and precursors lists into the subfolders 'peptides\\', 'precursors\\' and 'optimized_oe\\'.

### Perform direct infusion parameter optimisation

To get scan-by-scan targeting tables for the orbitrap ms run one of the files '!OE Assays\\!make_static_*_scan.qmd'.
To do so set 'peptide_set:' in the 'params:' section of the preamble and optionally adjust the settings in the 'settings' chunk.
The output tables will be placed in '!OE Assays\\static_\*_scan\\'.

Perform direct infusion experiments that cycle through all targets.

Evaluate the resulting .raw file using one of '!DI-MS Analyzers\\!README_\*_analysis.qmd'.
To do so set 'rawfile:' in the 'params:' section of the preamble and optionally adjust the settings in the 'settings' chunk.
The analyzer placed .pdfs and a .csv table under the name of the rawfile in the '!DI-MS Analyzers\\' folder.
Optionally, manually confirm the .csv contents by checking the quality of the results in the .pdf files.
Submit optimized values by copying the .csv to '!Peptide Lists\\optimized_oe\\readin' and running '!Peptide Lists\\!make_optimizations_oe.qmd'.

### Set up a inclusion list DDA experiment to find peptides in LC-MS

Run '!OE Assays\\!make_inclusion.qmd' and use the generated inclusion list table with a DDA experiment.
To do so set the 'peptide_set' variable and optionally further settings in the 'settings' and 'filters' chunks.
This places the inclusion list into the subfolder 'inclusion_lists\\'.

### Submit rt information for peptides

Set up a id_all.sky file with iRT peptides and target peptides in a dedicated folder.
Use '!export_precursor_data.cmd' in the folder with the .sky file.

Import the information by evaluating '!Peptide Lists\\!submit_rt_oe.qmd'.
To do so set 'detection_folder:' in the 'params:' section of the preamble.

### Schedule a PRM experiment

Run '!OE Assays\\!make_schedule.qmd' use the generated scheduling table with a PRM experiment.
Typically, a recent iRT reference file is used to produce accurate scheduling.
Confirm this by setting the 'irt_latest_folder' variable to an appropriate folder.'
Targets can be mixed and matched from the peptide lists in the '!Peptide Lists' folder.
A heavy and a light schedule are generated.
The heavy schedule is typically used to acquire reference libraries using heavy synthetic peptides.
The light schedule is typically used for target peptide detection runs.

### Evaluate a reference libary PRM run

Create a new subfolder somewhere in the main project directory and create your id_all.sky file there.
Import results using a skyline file with neutral losses enabled by default, ions from ion 1 to last ion and complete ion types (here y, b, a, p).
See '!Skyline Templates\spectral_library.sky' for example.
Ensure a set of retention time peptides is detected, too.
Copy evaluation helpers from '!Reports\': '!export_speclib_data.cmd', '!gen_sequence_coverage_plots.cmd'.
Then export a spectral library and run '!export_speclib_data.cmd'.

### Evaluate a target peptide PRM run

Create a new subfolder somewhere in the main project directory and create your id_all.sky file there.
Import results using a skyline file with spectral library active for all peptides.
See '!Skyline Templates\prm_target_detection.sky' for example.
On initial assessment for detections, only leave plausible peaks detected.
Then change settings to not 'Auto-select all matching transitions' and expand transitions to top 25 from library, including y, b, a, p from ion 1 to last ion.
Now double check plausible peaks by activating all remaining transitions.
Ensure a set of retention time peptides is detected, too.
Copy evaluation helpers from '!Reports\': '!PRM_Report_TAG.yaml', '!make_reports.cmd', '!gen_sequence_coverage_plots.cmd'.
Fill out !PRM_Report_TAG.yaml' and run '!make_reports.cmd'.
