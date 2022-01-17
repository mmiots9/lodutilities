<h1 style='text-align:center'> lodutilities package 🔬📊</h1>
---

Hi everyone 👋🏻!<br>
In this package, you can find some useful functions and files to perform your analysis. 

<h2>Installation</h2>
To install the package run <code>devtools::install_github("lab-lodato/lodutilities")</code>.  
It will install the package as well as its dependencies.

<h2>Functions</h2>
Here, you can find the explanation of the various functions.

<h3>PCRAnalysis</h3>
This function allows to analyze qPCR data. In particular, it takes as input the output file (xlsx) of the analysis and returns a sheet (on the same file or on another) containing a table for each sample with these columns:  
target | Rep1...Repn | Mean | sd | deltaCT | foldchange

<h5>Input</h5>
The inputs of the function are:
<ul>
<li><b>inputfile:</b> String containing the name of the input file. If left NA (default), a prompt will ask you to select the file</li>
<li><b>outputfile:</b> String containing the name of the output file. If left NA (default), the outputfile will be the same as the inputfile (NOT overwriting it)</li>
<li><b>max_rep_diff:</b> Number indicating the max range between min and max CT replicate (default 0.8)</li>
<li><b>housekeeping:</b> Name of the housekeeping gene (default GADPH)</li>
<li><b>sheetname:</b> String containing the name of the new sheet (default Analysis). It should NOT exist in outputfile yet</li>
</ul>

<h5>Output</h5>
The output is an xlsx sheet, containing one table for each sample with these columns: <br>
target | Rep1...Repn | Mean | sd | deltaCT | foldchange

deltaCT and foldchange are calculated from the housekeeping gene inserted as input. Replicates which exceed max_rep_diff are deleted and not used for calculations.
