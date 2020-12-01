# LBPA
Length-Based Pseudo-cohort Analysis (LBPA)
(https://doi.org/10.1016/j.fishres.2020.105810)

Stock status for many medium- and small-scale fisheries is unknown due, for example, to a lack of catch data and
the absence of scientific observer programs. However, length-frequency data are often available for such fisheries
because they are the cheapest and easiest data to obtain. Various stock assessment methods have been developed
that use length-frequency data and make equilibrium assumptions regarding both recruitment and fishing
mortality. These assumptions raise questions regarding the reliability of the results, particularly when the
method is applied to a single sample of length-frequency. We developed a Length-Based Pseudo-cohort Analysis
(LBPA) model whose parameters can be estimated using multiple length frequencies and penalized maximum
likelihood, under the assumption that using more than one length-frequency sample reduces the effects of the
equilibrium conditions assumed in the model. This work provides guidelines that should
be considered when using length-based pseudo-cohort models for data-poor fisheries.

When run the model, be keep to have into same directory/folder next files

- LBPA.exe
- LBPA.dat
- LBPA.tpl

*.tpl file must be compiled in ADMB for producing the executable file (lbpa.exe). Alternatively you will find the executable file for windows 64 bits with other extention (.mexe). Be keep to rename the extension file .mexe by .exe. Depending on your OS, you can compile the source code (.tpl) as be necessary

After running the model (double click on lbpa.exe) you will find several files, being the most interest: lbpa.rep (results report), lbpa.par, lbpa.std and lbpa.cor.  

Data file (lbpa.dat) is explained by itself. There are several comments related to parameters and data inside.

For any question: cristian.canales.r@pucv.cl
