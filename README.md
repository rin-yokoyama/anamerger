# anamerger

This project was originally developped for BRIKEN experiment.

### Required dependencies

* CERN ROOT (root-project/root) v6.14.00+
* yaml-cpp (jbeder/yaml-cpp) v0.6+
* BrikenTools
* cmake3
* c++11

### How to build


```
mkdir build
cd build
cmake [-DBRIKEN_TOOLS_DIR='path_to_your_BrikenTools_dir'] [-DMERGER_WORK_DIR='path_to_your_work_dir']
make install
```

Make sure to provide your working directory for the merger. It will generate script files to run the programs in your working directory. 
The default install path is ${CMAKE_SOURCE_DIR}/install
ROOT dictionary will also be installed in install/lib

### How to use script files

Following files will be generated under the install directory.
* config/
* scripts/
* mergerLists/
* reference_cuts/

Copy these directories into your working directory.
Then run following scripts to generate config files after changing tmplate yaml files if necessary.
* config/configallnuclides.sh
```
cd config
chmod a+x config_gen.sh
sh configallnuclides.sh
```
anamerger_config.yaml is the template file. When you want to change any common parameters, edit this file before you run configallnuclides.sh

Create a directory for log outputs named "logs" in your working directory.

Generate a mergerList file.

mergerList files are the list of input file names with the absolute paths. You can change the file names or the file paths to match your work environment.

When you want to submit jobs on a computing facility,
```
qsub scripts/run_anamerger.pbs
```
You can check the IDs of running jobs by qstat command.

### How to run the excutable

There will be a module file for Environment Modules installed in install/share/modulefiles/
You can use
```
module use -a "path-to-install/share/modulefiles"
module load anamerger
```
to setup environment.
Or you can simply source the setup.sh
```
source install/share/setup.sh
```

anamerger usage:
```
anamerger_main -o [output_rootfile] -c [config yaml file]
```

### Parameters in the config file 

* TreeName: Name of the input tree
* MergerListName: Name of the mergerList.txt file (relative path)
* ReferenceCutsName: Name of the reference cut file (absolute path when use proof)
* UseProof: If True, run on a proof lite server. Otherwise, run in a single thread.
* NumWorkers: Number of parallel workers. Available if UseProof is True.
* PrintFrequency: Defines how often the program to print the progress.
* HistogramGroups: A list of histogram group names.
* - PID: Particle identification plots gated by implant layers
* - ISOMER: Isomer gamma spectra
* - BETA: histograms related to beta events
* - DECAYCURVE: decay curves (requires BETA)
* - GAMMA_ET: Egamma vs Tbeta-Timp histograms (requires BETA)
* - GAMMA_GAMMA: gamma-gamma (requires BETA)
