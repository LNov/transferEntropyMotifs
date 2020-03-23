Instructions:

1) Open FigX folder (where X is the figure number)

2) Run `pypet_master_figX.py` (where X is the figure number) to generate the raw results. A new time-stamped subfolder will be generated (we'll call it `<my_new_subfolder>` in the following) and the raw results will be stored in the `traj.hdf5` file. This can take a while, especially for Fig2 and Fig4. To save time, skip this step and unzip the raw results (`traj.zip`) in the existing subfolder.

4) Run `pypet_postprocessing_figX.py` to postprocess the raw results. This will generate a pandas DataFrame (`postprocessing.pkl`) which will be stored in the same subfolder as the raw results. If you generated the raw results from scratch and they are stored in a new subfolder, you first need to edit `pypet_postprocessing_figX.py` and set `traj_dir` as the new subfolder name (e.g. `traj_dir = '<my_new_subfolder>'`).

5) Run `plotting_figX.py` to generate the plot. The new figure will be stored in the same subfolder as the postprocessed data (as a PDF file). Once again, if you generated the raw results from scratch and they are stored in a new subfolder, you first need to edit `plotting_figX.py` and set `traj_dir` as the new subfolder name (e.g. `traj_dir = '<my_new_subfolder>'`).