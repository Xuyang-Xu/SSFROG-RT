from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable, Optional
import time

import numpy as np
import tifffile
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

def save_array_as_tiff(
    array: np.ndarray,
    name: str,
    folder: str | Path,
    *,
    dtype: np.dtype | None = None,
    overwrite: bool = True,
) -> Path:
    """
    Save a 2D NumPy array as a TIFF file.

    Parameters
    ----------
    array : np.ndarray
        Input 2D image array.
    name : str
        Output file name, with or without '.tiff' / '.tif'.
    folder : str or Path
        Folder where the TIFF file will be saved.
    dtype : np.dtype or None
        Optional output dtype. If None, keep original dtype when possible.
    overwrite : bool
        Whether to overwrite an existing file.

    Returns
    -------
    Path
        Full path to the saved TIFF file.
    """
    array = np.asarray(array)

    if array.ndim != 2:
        raise ValueError(f"input_array must be 2D, got shape {array.shape}")

    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)

    name = str(name)
    if not name.lower().endswith((".tif", ".tiff")):
        name = f"{name}.tiff"

    save_path = folder / name

    if save_path.exists() and not overwrite:
        raise FileExistsError(f"File already exists: {save_path}")

    out = array
    if dtype is not None:
        out = out.astype(dtype)

    tifffile.imwrite(save_path, out)

    return save_path

def _load_matlab_engine():
    try:
        import matlab
        import matlab.engine
        return matlab, matlab.engine
    except Exception as e:
        raise ImportError(
            "MATLAB Engine for Python is required in the current Python environment."
        ) from e


def _matlab_array_to_numpy(x) -> np.ndarray:
    """
    Convert a MATLAB array returned by matlab.engine to a NumPy array
    while preserving MATLAB column-major ordering.
    """
    if hasattr(x, "_data") and hasattr(x, "size"):
        return np.array(x._data, dtype=float).reshape(x.size, order="F")
    return np.array(x, dtype=float)


def _get_real_from_workspace(eng, var_name: str) -> np.ndarray:
    return _matlab_array_to_numpy(eng.workspace[var_name])


def _get_complex_from_workspace(eng, var_name: str) -> np.ndarray:
    real_part = _matlab_array_to_numpy(eng.eval(f"real({var_name});", nargout=1))
    imag_part = _matlab_array_to_numpy(eng.eval(f"imag({var_name});", nargout=1))
    return real_part + 1j * imag_part


def _as_1d(x: np.ndarray) -> np.ndarray:
    return np.squeeze(np.asarray(x))


def run_RT(
    file_dir: str | Path | None = None,
    search_str: str | None = None,
    *,
    input_array: np.ndarray | None = None,
    input_name: str | None = None,
    input_save_dir: str | Path | None = None,
    input_dtype: np.dtype | None = None,
    overwrite_input_tiff: bool = True,
    max_iter: int = 500,
    plot_iter: int = 1000,
    temporal_calibration: float = 2.0,
    spectral_calibration: float = 0.4614,
    background: float = 5.0,
    central_wavelength: float = 386.78,
    sz: int = 256,
    width: int = 70,
    file_index: int = 0,
    matlab_code_dir: str | Path | None = None,
    extra_matlab_paths: Optional[Iterable[str | Path]] = None,
    eng=None,
    close_engine: bool = False,
    verbose: bool = True,
) -> tuple[dict[str, Any], Any]:
    """
    Run one RT/FROG retrieval by calling the MATLAB function
    Copy_of_GrenouilleAlgorithm_V2 from Python.

    Two input modes are supported:

    Mode 1: Search an existing TIFF
        file_dir + search_str

    Mode 2: Create a TIFF from a NumPy array, save it, then run retrieval
        input_array + input_name + input_save_dir

    The reported runtime excludes plotting.
    """
    matlab, matlab_engine = _load_matlab_engine()

    # Validate mutually exclusive input modes
    use_array_mode = input_array is not None

    if use_array_mode:
        if input_name is None or input_save_dir is None:
            raise ValueError(
                "When input_array is provided, input_name and input_save_dir "
                "must also be provided."
            )
    else:
        if file_dir is None or search_str is None:
            raise ValueError(
                "Either provide (file_dir and search_str), or provide "
                "(input_array, input_name, input_save_dir)."
            )

    owns_engine = eng is None
    if eng is None:
        eng = matlab_engine.start_matlab()

    t0 = time.perf_counter()

    try:
        if matlab_code_dir is not None:
            matlab_code_dir = str(Path(matlab_code_dir).expanduser().resolve())
            eng.addpath(eng.genpath(matlab_code_dir), nargout=0)
            eng.cd(matlab_code_dir, nargout=0)

        if extra_matlab_paths:
            for p in extra_matlab_paths:
                p = str(Path(p).expanduser().resolve())
                eng.addpath(eng.genpath(p), nargout=0)

        # -------------------------------------------------
        # Input mode selection
        # -------------------------------------------------
        if use_array_mode:
            saved_tiff_path = save_array_as_tiff(
                array=input_array,
                name=input_name,
                folder=input_save_dir,
                dtype=input_dtype,
                overwrite=overwrite_input_tiff,
            )

            file_dir = str(saved_tiff_path.parent)
            search_str = saved_tiff_path.name

            if verbose:
                print(f"Saved input TIFF: {saved_tiff_path}")

        trace_pattern = str((Path(file_dir).expanduser() / search_str).resolve())

        eng.workspace["trace_pattern_py"] = trace_pattern
        eng.workspace["file_index_py"] = float(file_index + 1)

        eng.workspace["MaxIter_py"] = matlab.double([float(max_iter)])
        eng.workspace["PlotIter_py"] = float(plot_iter)
        eng.workspace["TemporalCalibration_py"] = float(temporal_calibration)
        eng.workspace["SpectralCalibration_py"] = float(spectral_calibration)
        eng.workspace["BackGround_py"] = float(background)
        eng.workspace["CentralWavelength_py"] = float(central_wavelength)
        eng.workspace["SZ_py"] = float(sz)
        eng.workspace["WIDTH_py"] = float(width)

        matlab_cmd = """
        TifFiles = dir(trace_pattern_py);
        assert(~isempty(TifFiles), ['No file matched: ', trace_pattern_py]);

        [~, idx] = sort([TifFiles.datenum], 'ascend');
        TifFiles = TifFiles(idx);

        assert(numel(TifFiles) >= round(file_index_py), ...
            'file_index is out of range for the matched files.');
        TraceFileInfo = TifFiles(round(file_index_py));

        [Asig, Esig, tau, Et, lam1, lam2, w1, w2, Ew, FWHMt, FWHMw, Z, G] = ...
            Copy_of_GrenouilleAlgorithm_V2( ...
                MaxIter_py, PlotIter_py, TemporalCalibration_py, ...
                SpectralCalibration_py, BackGround_py, ...
                CentralWavelength_py, SZ_py, WIDTH_py, TraceFileInfo);
        """
        eng.eval(matlab_cmd, nargout=0)

        Asig = _get_real_from_workspace(eng, "Asig")
        Esig = _get_complex_from_workspace(eng, "Esig")
        tau = _as_1d(_get_real_from_workspace(eng, "tau"))
        Et = _as_1d(_get_complex_from_workspace(eng, "Et"))
        lam1 = _as_1d(_get_real_from_workspace(eng, "lam1"))
        lam2 = _as_1d(_get_real_from_workspace(eng, "lam2"))
        w1 = _as_1d(_get_real_from_workspace(eng, "w1"))
        w2 = _as_1d(_get_real_from_workspace(eng, "w2"))
        Ew = _as_1d(_get_complex_from_workspace(eng, "Ew"))
        FWHMt = float(eng.workspace["FWHMt"])
        FWHMw = float(eng.workspace["FWHMw"])
        Z = _as_1d(_get_real_from_workspace(eng, "Z"))
        G = _as_1d(_get_real_from_workspace(eng, "G"))
        selected_file = eng.eval("TraceFileInfo.name;", nargout=1)

        It = np.abs(Et) ** 2
        Iw = np.abs(Ew) ** 2

        jetmax = int(np.argmax(It))
        jewmax = int(np.argmax(Iw))

        etmax = float(It[jetmax])
        ewmax = float(Iw[jewmax])

        phiet = -np.unwrap(np.angle(Et))
        phiew = -np.unwrap(np.angle(Ew))

        phiet_rel = phiet - phiet[jetmax]
        phiew_rel = phiew - phiew[jewmax]

        Esig_abs = np.abs(Esig)
        difference = Esig_abs - Asig

        G_curve = np.asarray(G, dtype=float)
        Z_curve = np.asarray(Z, dtype=float)

        final_G = float(G_curve[-1]) if G_curve.size else np.nan
        best_G = float(np.min(G_curve)) if G_curve.size else np.nan
        final_Z = float(Z_curve[-1]) if Z_curve.size else np.nan

        rmse_trace = float(np.sqrt(np.mean((Esig_abs - Asig) ** 2)))
        mae_trace = float(np.mean(np.abs(Esig_abs - Asig)))

        runtime_seconds = time.perf_counter() - t0

        result = {
            "selected_file": selected_file,
            "trace_pattern": trace_pattern,
            "input_mode": "array_to_tiff" if use_array_mode else "existing_tiff",
            "generated_input_tiff": str(saved_tiff_path) if use_array_mode else None,
            "params": {
                "max_iter": max_iter,
                "plot_iter": plot_iter,
                "temporal_calibration": temporal_calibration,
                "spectral_calibration": spectral_calibration,
                "background": background,
                "central_wavelength": central_wavelength,
                "sz": sz,
                "width": width,
                "file_index": file_index,
            },
            "Asig": Asig,
            "Esig": Esig,
            "Esig_abs": Esig_abs,
            "difference": difference,
            "tau": tau,
            "Et": Et,
            "It": It,
            "It_norm": It / etmax if etmax != 0 else It,
            "phiet": phiet,
            "phiet_rel": phiet_rel,
            "lam1": lam1,
            "lam2": lam2,
            "w1": w1,
            "w2": w2,
            "Ew": Ew,
            "Iw": Iw,
            "Iw_norm": Iw / ewmax if ewmax != 0 else Iw,
            "phiew": phiew,
            "phiew_rel": phiew_rel,
            "FWHMt": FWHMt,
            "FWHMw": FWHMw,
            "Z_curve": Z_curve,
            "G_curve": G_curve,
            "final_G": final_G,
            "best_G": best_G,
            "final_Z": final_Z,
            "rmse_trace": rmse_trace,
            "mae_trace": mae_trace,
            "runtime_seconds": runtime_seconds,
            "etmax": etmax,
            "ewmax": ewmax,
            "jetmax": jetmax,
            "jewmax": jewmax,
        }

        if verbose:
            print(f"Selected file: {selected_file}")
            print(f"Runtime (no plotting): {runtime_seconds:.3f} s")
            print(f"Final G: {final_G:.6g}")
            print(f"Best G : {best_G:.6g}")

        return result, eng

    finally:
        if close_engine and owns_engine:
            eng.quit()

def load_dat_colormap(dat_file: str | Path) -> ListedColormap:
    """
    Load a matplotlib colormap from a .dat file with RGB values in [0, 1].
    """
    dat_file = Path(dat_file)
    colors = np.loadtxt(dat_file)
    return ListedColormap(colors, name=dat_file.stem)

def plot_rt_summary(
    result: dict[str, Any],
    *,
    figsize: tuple[float, float] = (18, 10),
    trace_cmap="viridis",
    diff_cmap="RdBu_r",
    trace_clim: tuple[float, float] | None = (0.0, 1.0),
    diff_clim: tuple[float, float] | None = None,
    temporal_xlim: tuple[float, float] | None = None,
    spectral_xlim: tuple[float, float] | None = None,
    save_dir: str | Path | None = None,
    filename: str | None = None,
    dpi: int = 150,
) -> tuple[plt.Figure, np.ndarray] | None:
    """
    Plot a 2x3 summary figure or save it to disk.

    If save_dir is provided, the figure is saved and not displayed.
    """
    tau = result["tau"]
    lam1 = result["lam1"]
    lam2 = result["lam2"]

    Asig = result["Asig"]
    Esig_abs = result["Esig_abs"]
    difference = result["difference"]

    It_norm = result["It_norm"]
    Iw_norm = result["Iw_norm"]
    phiet_rel = result["phiet_rel"]
    phiew_rel = result["phiew_rel"]

    FWHMt = result["FWHMt"]
    FWHMw = result["FWHMw"]
    G_curve = result["G_curve"]
    final_G = result.get("final_G", np.nan)
    runtime_seconds = result.get("runtime_seconds", np.nan)
    selected_file = result.get("selected_file", "RT_result")

    fig, axes = plt.subplots(2, 3, figsize=figsize, constrained_layout=True)

    # ===== Row 1 =====
    # Original
    ax = axes[0, 0]
    im0 = ax.imshow(
        Asig,
        extent=[tau.min(), tau.max(), lam2.min(), lam2.max()],
        origin="lower",
        aspect="auto",
        cmap=trace_cmap,
    )
    if trace_clim:
        im0.set_clim(*trace_clim)
    ax.set_title("Original Trace")
    fig.colorbar(im0, ax=ax)

    # Reconstructed
    ax = axes[0, 1]
    im1 = ax.imshow(
        Esig_abs,
        extent=[tau.min(), tau.max(), lam2.min(), lam2.max()],
        origin="lower",
        aspect="auto",
        cmap=trace_cmap,
    )
    if trace_clim:
        im1.set_clim(*trace_clim)
    ax.set_title("Reconstructed Trace")
    fig.colorbar(im1, ax=ax)

    # Difference
    ax = axes[0, 2]
    im2 = ax.imshow(
        difference,
        extent=[tau.min(), tau.max(), lam2.min(), lam2.max()],
        origin="lower",
        aspect="auto",
        cmap=diff_cmap,
    )
    vmax = np.max(np.abs(difference))
    if vmax > 0:
        im2.set_clim(-vmax, vmax)
    ax.set_title("Difference")
    fig.colorbar(im2, ax=ax)

    # ===== Row 2 =====
    # Temporal
    ax = axes[1, 0]
    ax2 = ax.twinx()
    ax.plot(tau, It_norm)
    ax2.plot(tau, phiet_rel, "--")
    ax.set_xlim(-200, 200)
    ax.set_title(f"Temporal | FWHM={FWHMt:.2f} fs")

    # Spectral
    ax = axes[1, 1]
    ax2 = ax.twinx()
    ax.plot(lam1, Iw_norm)
    ax2.plot(lam1, phiew_rel, "--")
    ax.set_xlim(700, 900)
    ax2.set_xlim(ax.get_xlim())
    ax.set_title(f"Spectral | FWHM={FWHMw:.2f} nm")

    # G vs iteration
    ax = axes[1, 2]
    if G_curve.size > 0:
        ax.plot(np.arange(1, len(G_curve)+1), G_curve)
    ax.set_title(f"G vs Iter | Final={final_G:.4g}")
    ax.grid(True)

    fig.suptitle(
        f"{selected_file} | Runtime={runtime_seconds:.3f}s",
        fontsize=13
    )

    # ===== SAVE OR SHOW =====
    if save_dir is not None:
        save_dir = Path(save_dir)
        save_dir.mkdir(parents=True, exist_ok=True)

        if filename is None:
            # auto filename
            filename = f"{Path(selected_file).stem}_RT_summary.png"

        save_path = save_dir / filename

        fig.savefig(save_path, dpi=dpi)
        plt.close(fig)   # IMPORTANT: avoid memory leak in loops
        return None

    else:
        return fig, axes