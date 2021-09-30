# Directory structure

- data/ : data used as input in the openmole scripts
- models/ : netlogo models used in the scripts
- openmole/ : openmole scripts
- output/ : output of simulation experiments (not all openmole output, but only the files that are actually used to make the reports).
- reports/: scripts to generate figures and reports
- post_proc/: scripts to analyse results from OpenMOLE


# Usage

## Reports

`reports/abc.Rmd` plots the posterior distributions from the results of the openmole script. To use it, edit the line preamble section at the top of the file, in particular the datafile to set it to the most up to date posterior sample (last file that was output by the openmole script).

To generate an html file, run in a shell

```sh
Rscript render.R
```



