# ivtFOODIE Kd fitting

Start with the Snakemake workflow and sample-sheet instructions in [docs/README.md](docs/README.md). [Input and output notes](docs/INPUT_OUTPUT.md) describe the expected ratio, concentration and motif files. Set the reference FASTA, sample ratio inputs, output and cluster paths in the selected configuration. For the running Kd helper, pass a reference FASTA index or chromosome sizes file as the first argument to the shell script.

This subproject follows the xielab working-tree workflow, including the conversion quality-control rule. Mapping ratio outputs are the inputs to this workflow.
