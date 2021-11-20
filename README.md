# Riboregulator Designer

Riboregulator 2.0 is a new design of riboregulators, that aim to overcome the limited specificty of toehold switches by placing the trigger binding site in the loop region of the hairping structure. The target RNA must compete with the stem of the hairpin for binding, this competition strongly disfavours mismatches, analogous to molecular beacons and may be able to discriminate between RNAs with a single nucleotide difference. Note, these riboregulators are currently untested in vitro.  

This tool is for designing and analysing riboregulator 2.0 structures for a specified RNA sequence. This tool runs in python and employs the NUPACK software to analyse RNA structures. 

The design of the riboregulators can be changed to design any desired riboregulator structure if so desired. 

## Dependencies

See "Using Riboregulator Designer" for more information. 

### NUPACK 3.0 or later (3.1, 3.2)

NUPACK  can be downloaded from http://nupack.org/.

As NUPACK only runs reliably on Linux, it is recommeneded that Windows users run via Windows Subsystem for Linux (WSL) 2, which we have verified to work without any problems. The documentation for WSL is here: https://docs.microsoft.com/en-us/windows/wsl/

### Python 3/Jupyter

The methods are demonstrated in a Jupyter notebook. Once you have python 3 installed, you can install the requirements using 

`pip install -r requirements.txt` 


## References

Large parts of the code are based on the code written in https://github.com/elanibal/NupackSensors.

The python wrappers are modified from https://github.com/DNA-and-Natural-Algorithms-Group/multistrand to run in python 3 (see nupack folder).

The scoring metrics for the designs were taken from  Ma, D. et al (Low-cost detection of norovirus using paper-based cell-free systems and synbody-based viral enrichment. Synth.  Biol.3, ysy018 (2018)), however, the scoring equation is for toehold switches and so a scoring equation for the the riboregulator 2.0 design is needed.  
