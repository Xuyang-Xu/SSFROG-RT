from RT_utils import run_RT, plot_rt_summary, load_dat_colormap
import matplotlib.pyplot as plt

custom_cmap = load_dat_colormap(r"D:\FROG_RT_Xuyang\color_maps\Wh_rainbow.dat")

result, eng = run_RT(
    file_dir=r"Z:\Data\2026_02_15_PlasmaMirrors_PhaseMeasurement\kHz_FROG_ExpTable",
    search_str=r"duration_kHz_2ndorder_37800fs^2_3rdorder_-54760.0fs^3_20260215_182358.Tiff",
    matlab_code_dir=r"Z:\Code\Experimental_Analysis\FROG\Frog_Automated",
    max_iter=200,
    plot_iter=1000,
    temporal_calibration=2.0,
    spectral_calibration=0.4614,
    background=5,
    central_wavelength=386.78,
    sz=256,
    width=70,
)

print("Selected file:", result["selected_file"])
print("Final G:", result["final_G"])
print("Best G :", result["best_G"])

fig, axes = plot_rt_summary(
    result,
    trace_cmap= custom_cmap,
    diff_cmap="RdBu_r",
    temporal_xlim=(-150, 150),
    spectral_xlim=(725, 875),
)
plt.show()

eng.quit()