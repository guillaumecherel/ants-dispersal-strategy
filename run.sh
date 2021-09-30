#!/usr/bin/bash

set -euo pipefail

echo "Running openmole ABC job."
# openmole-automate openmole-job-abc.toml

echo "Creating marginals figures."
Rscript reports/result_abc_marginals.R

echo "Resampling."
Rscript reports/result_abc_resample.R

echo "Running openmole resampling job."
cp output/resampleABCs_4_params_v3.csv openmole/
# openmole-automate openmole-job-resample.toml
rm openmole/resampleABCs_4_params_v3.csv

echo "Computing errors figures."
Rscript reports/result_abc_errors.R

