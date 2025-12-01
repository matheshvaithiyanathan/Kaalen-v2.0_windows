# Copyright (c) [2025] [Mathesh Vaithiyanathan]
# This software is licensed under the MIT License.
# See the LICENSE file for details.
import resources_rc
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import lmfit
from scipy.linalg import lstsq
from scipy.optimize import curve_fit
from scipy.special import erf
from functools import partial
import scipy
import re
from PyQt5.QtWidgets import (QApplication, QWidget, QGridLayout, QLabel, QLineEdit,
                             QPushButton, QCheckBox, QMessageBox, QVBoxLayout, QHBoxLayout, QTextEdit, QSpinBox,
                             QMainWindow, QSlider, QDialog, QSplitter, QDialogButtonBox, QInputDialog, QDoubleSpinBox,
                             QMenu, QAction, QFileDialog, QComboBox, QTabWidget, QGroupBox)
from PyQt5 import QtCore
from PyQt5.QtCore import QThread, pyqtSignal, Qt
from PyQt5.QtGui import QFont, QColor, QDoubleValidator, QIntValidator, QIcon
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import json
import os
import ctypes
import pyqtgraph as pg
import matplotlib
import pandas as pd
from scipy.interpolate import RectBivariateSpline, UnivariateSpline, interp1d, griddata
from lmfit.printfuncs import report_fit

print(f"Matplotlib Version: {matplotlib.__version__}")
print(f"PyQt5 Version: {QtCore.PYQT_VERSION_STR}")
print(f"Pandas Version: {pd.__version__}")
print(f"NumPy Version: {np.__version__}")
print(f"PyQtGraph Version: {pg.__version__}")
print(f"Scipy Version: {scipy.__version__}")
print(f"LMFIT Version: {lmfit.__version__}")


def find(in_array, target_value):
    array = in_array
    nearest_index = np.abs(array - target_value).argmin( )
    nearest_value = array [ nearest_index ]
    return nearest_index


def simple_multi_exponential_gf(time, taus):
    return np.exp(-time [ :, None ] / taus)


def convolved_exponential_analytical(time, tau, t0, delta):
    k = 1.0 / tau
    mu = t0
    delta_tilde = delta / (2 * np.sqrt(2 * np.log(2)))

    model = 0.5 * np.exp(-k * time) * np.exp(k * (mu + k * delta_tilde ** 2 / 2)) * (
            1 + erf((time - (mu + k * delta_tilde ** 2)) / (np.sqrt(2) * delta_tilde)))
    model [ time < t0 - 5 * delta_tilde ] = 0
    return model


def _format_unit_for_display(unit_string):
    unit_string = unit_string.replace("^-1", "\u207B\u00B9")
    unit_string = unit_string.replace("^-2", "\u207B\u00B2")
    unit_string = unit_string.replace("^2", "\u00B2")
    unit_string = unit_string.replace("^3", "\u00B3")
    unit_string = unit_string.replace("^-3", "\u207B\u00B3")
    unit_string = unit_string.replace("^-4", "\u207B\u2074")
    unit_string = unit_string.replace("_1", "\u2081")
    unit_string = unit_string.replace("_2", "\u2082")
    unit_string = unit_string.replace("_3", "\u2083")
    return unit_string


def _parse_label_and_unit(label_string):
    match = re.search(r'^(.*?)\s*\[(.*?)\]', label_string)
    if match:
        label = match.group(1).strip( )
        unit = match.group(2).strip( )
        return label, unit
    return label_string, ''


def _build_single_component_pfid_terms(T_new, ν, T2, ν10, ν21, r):
    c = 2.998e10

    T_sec = T_new * 1e-12
    T2_sec = T2 * 1e-12

    exp_factor = np.exp(-T_sec [ :, None ] / T2_sec)

    term1_numerator = (1 / (T2 * 1e-12)) * np.cos(2 * np.pi * c * (ν - ν10) * T_sec [ :, None ]) - \
                      2 * np.pi * c * (ν - ν10) * np.sin(2 * np.pi * c * (ν - ν10) * T_sec [ :, None ])
    term1_denominator = (2 * np.pi * c * (ν - ν10)) ** 2 + (1 / T2_sec) ** 2
    term1 = exp_factor * (term1_numerator / term1_denominator)

    term2_numerator_orig = (1 / (T2 * 1e-12)) * np.cos(2 * np.pi * c * (ν - ν10) * T_sec [ :, None ]) - \
                           2 * np.pi * c * (ν - ν21) * np.sin(2 * np.pi * c * (ν - ν10) * T_sec [ :, None ])
    term2_denominator_orig = (2 * np.pi * c * (ν - ν21)) ** 2 + (1 / T2_sec) ** 2

    term2 = -r * exp_factor * (term2_numerator_orig / term2_denominator_orig)

    return term1.ravel( ), term2.ravel( )


def build_design_matrix_pfid(T_new, ν, component_params, r_values):
    all_terms = [ ]

    for i, (T2, ν10, ν21) in enumerate(component_params):
        r = r_values [ i ]

        term1_ravel, term2_ravel = _build_single_component_pfid_terms(T_new, ν, T2, ν10, ν21, r)

        all_terms.append(term1_ravel)
        all_terms.append(term2_ravel)

    n_T = T_new.size
    n_ν = ν.size
    offset_term = np.ones((n_T, n_ν)).ravel( )
    all_terms.append(offset_term)

    return np.column_stack(all_terms)


def residual_pfid(params, T_new, ν, data, num_components):
    component_params = [ ]
    r_values = [ ]

    for i in range(num_components):
        comp_idx = i + 1
        T2 = params [ f'T2_{comp_idx}' ].value
        ν10 = params [ f'ν10_{comp_idx}' ].value
        ν21 = params [ f'ν21_{comp_idx}' ].value
        r_values.append(params [ f'r_{comp_idx}' ].value)
        component_params.append((T2, ν10, ν21))

    A = build_design_matrix_pfid(T_new, ν, component_params, r_values)
    y = data.ravel( )

    amplitudes = lstsq(A, y) [ 0 ]

    model = A @ amplitudes
    return (model - y)


class AnalysisWorker(QThread):
    finished = pyqtSignal(object)
    error = pyqtSignal(str)

    def __init__(self, params):
        super( ).__init__( )
        self.params = params

    def run(self):
        try:
            results = self.run_global_fit_analysis(**self.params)
            self.finished.emit(results)
        except Exception as e:
            self.error.emit(str(e))

    def run_global_fit_analysis(self, time_min, time_max, probe_min, probe_max,
                                manual_num_components, probes_to_plot, use_convolved_model,
                                use_svd_initial_guess, manual_tau_guesses,
                                manual_t0_guess, manual_fwhm_guess, fix_t0, fix_fwhm,
                                x_axis, y_axis, two_d_spectrum, fixed_tau_indices, fixed_long_tau):

        if x_axis is None or y_axis is None or two_d_spectrum is None:
            self.error.emit("Data not provided to the analysis worker.")
            return

        try:
            time_min_idx = find(y_axis, time_min)
            time_max_idx = find(y_axis, time_max)
            probe_min_idx = find(x_axis, probe_min)
            probe_max_idx = find(x_axis, probe_max)

            data_sliced = two_d_spectrum [ time_min_idx:time_max_idx + 1, probe_min_idx:probe_max_idx + 1 ]
            time_sliced = y_axis [ time_min_idx:time_max_idx + 1 ]
            probe_sliced = x_axis [ probe_min_idx:probe_max_idx + 1 ]
        except IndexError:
            self.error.emit(
                "Error: The provided 'Probe Min' or 'Probe Max' values are outside the range of the loaded data. "
                "Please check your input values.")
            return

        initial_tau_guesses = [ fixed_long_tau ]
        warnings = [ ]

        if use_svd_initial_guess:
            U, S, Vt = np.linalg.svd(data_sliced, full_matrices = False)
            num_components_to_fit = manual_num_components

            for i in range(num_components_to_fit):
                initial_amp = np.sign(U [ 0, i ])
                p0 = [ initial_amp, 10.0 ]

                def mono_exponential_model(t, amplitude, tau):
                    return amplitude * np.exp(-t / tau)

                try:
                    popt, _ = curve_fit(mono_exponential_model, time_sliced, U [ :, i ], p0 = p0)
                    initial_tau_guesses.append(popt [ 1 ])
                except RuntimeError:
                    warnings.append(f"Warning: Fit for component {i + 1} failed. Skipping...")

            initial_tau_guesses = [ fixed_long_tau ] + sorted([ tau for tau in initial_tau_guesses [ 1: ] if tau > 0 ])
            initial_tau_guesses_for_print = initial_tau_guesses [ : ]
            if not initial_tau_guesses:
                self.error.emit("SVD-based initial guess generation failed. Check your data and fitting range.")
                return
        else:
            if manual_tau_guesses is None or len(manual_tau_guesses) != manual_num_components:
                self.error.emit("manual_tau_guesses must be a list with length equal to manual_num_components.")
                return
            initial_tau_guesses = [ fixed_long_tau ] + manual_tau_guesses
            initial_tau_guesses_for_print = initial_tau_guesses [ : ]

        def build_design_matrix(time, taus, t0, fwhm):
            if use_convolved_model:
                A_columns = [ convolved_exponential_analytical(time, tau, t0, fwhm) for tau in taus ]
            else:
                A_columns = [ simple_multi_exponential_gf(time, np.array([ tau ])).flatten( ) for tau in taus ]
            return np.stack(A_columns, axis = 1)

        def objective_function(params, time, data):
            taus = [ params [ f'tau_{i}' ].value for i in range(len(initial_tau_guesses)) ]
            t0 = params [ 't0' ].value if use_convolved_model else None
            fwhm = params [ 'fwhm' ].value if use_convolved_model else None

            A = build_design_matrix(time, taus, t0, fwhm)
            amplitudes, _, _, _ = lstsq(A, data)
            model = A @ amplitudes
            return (model - data).ravel( )

        params = lmfit.Parameters( )
        for i, tau_val in enumerate(initial_tau_guesses):
            if i == 0:
                vary_tau = False
            else:
                vary_tau = i - 1 not in fixed_tau_indices
            params.add(f'tau_{i}', value = tau_val, min = 0.01, max = np.inf, vary = vary_tau)

        if use_convolved_model:
            params.add('t0', value = manual_t0_guess, min = -5, max = 5, vary = not fix_t0)
            params.add('fwhm', value = manual_fwhm_guess, min = 0.01, max = 4, vary = not fix_fwhm)

        minimizer = lmfit.Minimizer(objective_function, params, fcn_args = (time_sliced, data_sliced))
        result = minimizer.minimize(method = 'leastsq')

        ss_res = np.sum(result.residual ** 2)
        ss_tot = np.sum((data_sliced - np.mean(data_sliced)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)

        best_taus = [ result.params [ f'tau_{i}' ].value for i in range(len(initial_tau_guesses)) ]
        best_t0 = result.params [ 't0' ].value if use_convolved_model else None
        best_fwhm = result.params [ 'fwhm' ].value if use_convolved_model else None

        A_final = lstsq(build_design_matrix(time_sliced, best_taus, best_t0, best_fwhm), data_sliced) [ 0 ]
        best_fit = (build_design_matrix(time_sliced, best_taus, best_t0, best_fwhm) @ A_final)

        import io
        import sys
        old_stdout = sys.stdout
        sys.stdout = buffer = io.StringIO( )
        try:
            report_fit(result)
            fit_report_string = buffer.getvalue( )
        finally:
            sys.stdout = old_stdout

        return (
            best_fit, A_final, best_taus, r_squared, probe_sliced, time_sliced, data_sliced, probes_to_plot, best_t0,
            best_fwhm, None, use_convolved_model, initial_tau_guesses_for_print, fit_report_string, warnings)


class PFIDFitWorker(QThread):
    finished = pyqtSignal(object)
    error = pyqtSignal(str)

    def __init__(self, params):
        super( ).__init__( )
        self.params = params

    def run(self):
        try:
            results = self.run_pfid_fit_analysis(**self.params)
            self.finished.emit(results)
        except Exception as e:
            self.error.emit(f"PFID Fit Error: {e}")

    def run_pfid_fit_analysis(self, time_min, time_max, probe_min, probe_max,
                              num_components, T2_params, nu10_params, nu21_params, r_params,
                              interp_method, num_interp_points,
                              x_axis, y_axis, two_d_spectrum):

        if x_axis is None or y_axis is None or two_d_spectrum is None:
            raise ValueError("Data not provided to the analysis worker.")

        T2_guesses, T2_fix_flags = T2_params
        nu10_guesses, nu10_fix_flags = nu10_params
        nu21_guesses, nu21_fix_flags = nu21_params
        r_guesses, r_fix_flags = r_params

        try:
            time_min_idx = find(y_axis, time_min)
            time_max_idx = find(y_axis, time_max)
            probe_min_idx = find(x_axis, probe_min)
            probe_max_idx = find(x_axis, probe_max)

            data_raw_sliced = two_d_spectrum [ time_min_idx:time_max_idx + 1, probe_min_idx:probe_max_idx + 1 ]
            time_raw_sliced = y_axis [ time_min_idx:time_max_idx + 1 ]
            probe_raw_sliced = x_axis [ probe_min_idx:probe_max_idx + 1 ]
        except IndexError:
            raise IndexError("The provided slicing range is outside the loaded data boundaries.")

        T_orig = np.abs(time_raw_sliced)
        ω_orig = probe_raw_sliced
        data_orig = data_raw_sliced

        T_grid, ω_grid = np.meshgrid(T_orig, ω_orig, indexing = 'ij')
        points = np.column_stack((T_grid.ravel( ), ω_grid.ravel( )))
        values = data_orig.ravel( )

        T_new = np.linspace(np.min(T_orig), np.max(T_orig), num_interp_points)
        ω_new = np.linspace(np.min(ω_orig), np.max(ω_orig), num_interp_points)
        T_new_grid, ω_new_grid = np.meshgrid(T_new, ω_new, indexing = 'ij')

        interpolated_data = griddata(
            points = points,
            values = values,
            xi = (T_new_grid, ω_new_grid),
            method = interp_method,
            fill_value = 0
        )
        data_interp = interpolated_data

        params = lmfit.Parameters( )

        for i in range(num_components):
            comp_idx = i + 1

            params.add(f'T2_{comp_idx}', value = T2_guesses [ i ], min = 0.01, max = 30,
                       vary = not T2_fix_flags [ i ])

            params.add(f'ν10_{comp_idx}', value = nu10_guesses [ i ], min = np.min(ω_new), max = np.max(ω_new),
                       vary = not nu10_fix_flags [ i ])

            params.add(f'ν21_{comp_idx}', value = nu21_guesses [ i ], min = np.min(ω_new), max = np.max(ω_new),
                       vary = not nu21_fix_flags [ i ])

            params.add(f'r_{comp_idx}', value = r_guesses [ i ], min = 0.0, max = 1.0,
                       vary = not r_fix_flags [ i ])

        minimizer = lmfit.Minimizer(residual_pfid, params, fcn_args = (T_new, ω_new, data_interp, num_components))
        result = minimizer.minimize(method = 'leastsq')

        best_params_list = [ ]
        best_r_list = [ ]
        for i in range(num_components):
            comp_idx = i + 1
            best_T2 = result.params [ f'T2_{comp_idx}' ].value
            best_nu10 = result.params [ f'ν10_{comp_idx}' ].value
            best_nu21 = result.params [ f'ν21_{comp_idx}' ].value
            best_r = result.params [ f'r_{comp_idx}' ].value
            best_params_list.append((best_T2, best_nu10, best_nu21))
            best_r_list.append(best_r)

        A_final = lstsq(build_design_matrix_pfid(T_new, ω_new, best_params_list, best_r_list),
                        data_interp.ravel( )) [ 0 ]

        best_fit_flat = build_design_matrix_pfid(T_new, ω_new, best_params_list, best_r_list) @ A_final
        best_fit = best_fit_flat.reshape(data_interp.shape)

        ss_res = np.sum((best_fit_flat - data_interp.ravel( )) ** 2)
        ss_tot = np.sum((data_interp - np.mean(data_interp)) ** 2)
        r_squared = 1 - (ss_res / ss_tot)

        import io
        import sys
        old_stdout = sys.stdout
        sys.stdout = buffer = io.StringIO( )
        try:
            report_fit(result)
            fit_report_string = buffer.getvalue( )
        finally:
            sys.stdout = old_stdout

        return (best_fit, A_final, r_squared, ω_new, T_new, data_interp, fit_report_string, result, num_components)


class PFIDFitterApp(QWidget):
    def __init__(self, main_window, x_axis_label='Probe wavenumber', y_axis_label='Time', x_axis_unit='cm\u207B\u00B9',
                 y_axis_unit='ps',
                 font_size=12):
        super( ).__init__(main_window)
        self.setWindowTitle("PFID (Photon-Frequency-ID) Fitting Analysis")
        self.setObjectName("PFID Fit")

        self.main_window = main_window

        self.x_axis_label = x_axis_label
        self.y_axis_label = y_axis_label
        self.x_axis_unit = x_axis_unit
        self.y_axis_unit = y_axis_unit
        self.font_size = font_size
        self.worker_thread = None

        self.x_axis_data = None
        self.y_axis_data = None
        self.two_d_spectrum_data = None

        self.initUI( )
        self.results_data = None

    def initUI(self):
        main_layout = QGridLayout(self)

        input_panel = QGroupBox("PFID Model Parameters & Range")
        input_layout = QVBoxLayout(input_panel)

        formatted_x_unit = _format_unit_for_display(self.x_axis_unit)
        formatted_y_unit = _format_unit_for_display(self.y_axis_unit)

        grid_range = QGridLayout( )
        grid_range.addWidget(QLabel(f"{self.y_axis_label} Min (must be negative, {formatted_y_unit}):"), 0, 0)
        self.time_min_input = QLineEdit("-15.0")
        grid_range.addWidget(self.time_min_input, 0, 1)
        grid_range.addWidget(QLabel(f"{self.y_axis_label} Max (must be negative, {formatted_y_unit}):"), 0, 2)
        self.time_max_input = QLineEdit("-0.5")
        grid_range.addWidget(self.time_max_input, 0, 3)

        grid_range.addWidget(QLabel(f"{self.x_axis_label} Min ({formatted_x_unit}):"), 1, 0)
        self.probe_min_input = QLineEdit("1775")
        grid_range.addWidget(self.probe_min_input, 1, 1)
        grid_range.addWidget(QLabel(f"{self.x_axis_label} Max ({formatted_x_unit}):"), 1, 2)
        self.probe_max_input = QLineEdit("1850")
        grid_range.addWidget(self.probe_max_input, 1, 3)
        input_layout.addLayout(grid_range)

        grid_model_comp = QGridLayout( )
        grid_model_comp.addWidget(QLabel("Number of PFID Components:"), 0, 0)
        self.num_components_input = QLineEdit("1")
        self.num_components_input.setValidator(QIntValidator(1, 5))
        grid_model_comp.addWidget(self.num_components_input, 0, 1)

        fix_help_label = QLabel("Note: Append ':' to a value to fix it (e.g., 2.0:, 1810:, 0.2:)")
        fix_help_label.setStyleSheet("color: blue;")
        grid_model_comp.addWidget(fix_help_label, 0, 2, 1, 2)

        input_layout.addLayout(grid_model_comp)

        grid_guesses = QGridLayout( )
        grid_guesses.addWidget(QLabel("T₂ Guesses (ps, comma-separated):"), 0, 0)
        self.T2_guess_input = QLineEdit("2.0")
        grid_guesses.addWidget(self.T2_guess_input, 0, 1, 1, 2)

        grid_guesses.addWidget(QLabel(f"ν₁₀ Guesses ({formatted_x_unit}, comma-separated):"), 1, 0)
        self.nu10_guess_input = QLineEdit("1810")
        grid_guesses.addWidget(self.nu10_guess_input, 1, 1, 1, 2)

        grid_guesses.addWidget(QLabel(f"ν₂₁ Guesses ({formatted_x_unit}, comma-separated):"), 2, 0)
        self.nu21_guess_input = QLineEdit("1785")
        grid_guesses.addWidget(self.nu21_guess_input, 2, 1, 1, 2)

        grid_guesses.addWidget(QLabel("r (Ratio) Guesses (comma-separated):"), 3, 0)
        self.r_guess_input = QLineEdit("0.2")
        grid_guesses.addWidget(self.r_guess_input, 3, 1, 1, 2)

        input_layout.addLayout(grid_guesses)

        interp_group = QGroupBox("Interpolation Settings (for fit matrix)")
        interp_layout = QGridLayout(interp_group)

        interp_layout.addWidget(QLabel("Method:"), 0, 0)
        self.interp_method_combo = QComboBox(self)
        self.interp_method_combo.addItems([ "linear", "cubic", "nearest" ])
        self.interp_method_combo.setCurrentText("cubic")
        interp_layout.addWidget(self.interp_method_combo, 0, 1)

        interp_layout.addWidget(QLabel("Points per Dimension:"), 1, 0)
        self.interp_points_input = QSpinBox(self)
        self.interp_points_input.setRange(20, 1000)
        self.interp_points_input.setValue(100)
        interp_layout.addWidget(self.interp_points_input, 1, 1)

        input_layout.addWidget(interp_group)

        button_layout = QHBoxLayout( )
        self.run_button = QPushButton("Run PFID Fit")
        self.run_button.clicked.connect(self.run_pfid_fit)
        button_layout.addWidget(self.run_button)

        self.export_button = QPushButton("Export Results")
        self.export_button.clicked.connect(self.export_fit_results)
        self.export_button.setDisabled(True)
        button_layout.addWidget(self.export_button)

        input_layout.addLayout(button_layout)

        main_layout.addWidget(input_panel, 0, 0, 1, 1, Qt.AlignTop)

        results_panel = QGroupBox("Fit Report")
        results_layout = QVBoxLayout(results_panel)
        self.results_text_edit = QTextEdit( )
        self.results_text_edit.setReadOnly(True)
        self.results_text_edit.setMinimumHeight(300)
        results_layout.addWidget(self.results_text_edit)

        main_layout.addWidget(results_panel, 0, 1, 1, 1, Qt.AlignTop)

        self.plot_container = QWidget( )
        self.plot_layout = QHBoxLayout(self.plot_container)
        main_layout.addWidget(self.plot_container, 1, 0, 1, 2)
        main_layout.setRowStretch(1, 10)
        self.setLayout(main_layout)

    def parse_fixed_float_list_gui(self, text, name, num_components):
        vals_str = [ v.strip( ) for v in text.split(',') if v.strip( ) ]

        if len(vals_str) != num_components:
            raise ValueError(
                f"Number of {name} guesses ({len(vals_str)}) must match number of components ({num_components}).")

        values = [ ]
        fix_flags = [ ]
        for v_str in vals_str:
            is_fixed = False
            v_str_clean = v_str
            if v_str.endswith(':'):
                is_fixed = True
                v_str_clean = v_str [ :-1 ]

            try:
                values.append(float(v_str_clean))
                fix_flags.append(is_fixed)
            except ValueError:
                raise ValueError(f"Invalid number format in {name}: '{v_str}' is not a valid float.")

        return values, fix_flags

    def run_pfid_fit(self):
        if self.main_window and hasattr(self.main_window, 'data_loaded') and self.main_window.data_loaded:
            self.x_axis_data = self.main_window.current_x_values
            self.y_axis_data = self.main_window.current_y_values
            self.two_d_spectrum_data = self.main_window.current_signal_data

        if self.x_axis_data is None:
            QMessageBox.critical(self, "No Data", "No 2D data is currently loaded in the main application for fitting.")
            return

        try:
            num_components = int(self.num_components_input.text( ))
            if num_components < 1:
                raise ValueError("Number of PFID Components must be 1 or greater.")

            T2_params = self.parse_fixed_float_list_gui(self.T2_guess_input.text( ), "T2", num_components)
            nu10_params = self.parse_fixed_float_list_gui(self.nu10_guess_input.text( ), "ν10", num_components)
            nu21_params = self.parse_fixed_float_list_gui(self.nu21_guess_input.text( ), "ν21", num_components)
            r_params = self.parse_fixed_float_list_gui(self.r_guess_input.text( ), "r", num_components)

            params = {
                "time_min": float(self.time_min_input.text( )),
                "time_max": float(self.time_max_input.text( )),
                "probe_min": float(self.probe_min_input.text( )),
                "probe_max": float(self.probe_max_input.text( )),
                "num_components": num_components,
                "T2_params": T2_params,
                "nu10_params": nu10_params,
                "nu21_params": nu21_params,
                "r_params": r_params,
                "interp_method": self.interp_method_combo.currentText( ),
                "num_interp_points": self.interp_points_input.value( ),
                "x_axis": self.x_axis_data,
                "y_axis": self.y_axis_data,
                "two_d_spectrum": self.two_d_spectrum_data,
            }

            if params [ "time_min" ] >= 0 or params [ "time_max" ] >= 0:
                QMessageBox.warning(self, "Input Warning",
                                    "PFID fit typically requires a **negative** time range (e.g., -15 to -0.5 ps). Please confirm your input is correct.")

            self.results_text_edit.setText(
                "Running fit... Please wait. This may take a moment due to grid interpolation and minimization.")
            self.run_button.setDisabled(True)
            self.export_button.setDisabled(True)

            self.worker_thread = PFIDFitWorker(params)
            self.worker_thread.finished.connect(self.plot_results)
            self.worker_thread.error.connect(self.handle_error)
            self.worker_thread.start( )

        except ValueError as ve:
            QMessageBox.critical(self, "Input Error", f"Invalid number format or parameter count: {ve}")
            self.run_button.setDisabled(False)
        except Exception as e:
            QMessageBox.critical(self, "Error", f"An unexpected error occurred: {e}")
            self.run_button.setDisabled(False)

    def plot_results(self, results):
        self.run_button.setDisabled(False)
        self.export_button.setDisabled(False)

        (best_fit, A_final, r_squared, probe, time, data_interp, fit_report_string, result, num_components) = results
        self.results_data = results

        A_offset = A_final [ -1 ]

        output_text = "--- PFID Fit Results ---\n\n"
        output_text += f"Number of Components Fitted: {num_components}\n"
        output_text += f"Overall R-squared value: {r_squared:.4f}\n"

        for i in range(num_components):
            A1 = A_final [ i * 2 ]
            A2 = A_final [ i * 2 + 1 ]
            output_text += f"Component {i + 1} Amplitudes: A₁ = {A1:.2g}, A₂ = {A2:.2g}\n"

        output_text += f"Overall Offset: {A_offset:.2g}\n\n"
        output_text += "--- Full lmfit Report ---\n"
        output_text += fit_report_string
        self.results_text_edit.setText(output_text)

        for i in reversed(range(self.plot_layout.count( ))):
            self.plot_layout.itemAt(i).widget( ).setParent(None)

        font_size = self.font_size
        plt.rcParams.update({ 'font.size': font_size })

        formatted_x_unit = _format_unit_for_display(self.x_axis_unit)
        formatted_y_unit = _format_unit_for_display(self.y_axis_unit)

        plot1_container = QWidget( )
        plot1_vbox = QVBoxLayout(plot1_container)

        fig1 = Figure(figsize = (5, 5), dpi = 100)
        canvas1 = FigureCanvas(fig1)
        toolbar1 = NavigationToolbar(canvas1, plot1_container)

        ax1 = fig1.add_subplot(111)

        vmin, vmax = data_interp.min( ), data_interp.max( )
        norm = plt.Normalize(vmin = vmin, vmax = vmax)

        im1 = ax1.pcolormesh(probe, time, data_interp, cmap = 'seismic', shading = 'auto', norm = norm)
        ax1.set_title('Interpolated Original Data')
        ax1.set_xlabel(f'{self.x_axis_label} ({formatted_x_unit})')
        ax1.set_ylabel(f'|{self.y_axis_label}| ({formatted_y_unit})')
        fig1.colorbar(im1, ax = ax1, label = 'Signal')
        fig1.tight_layout( )

        plot1_vbox.addWidget(toolbar1)
        plot1_vbox.addWidget(canvas1)
        self.plot_layout.addWidget(plot1_container)

        plot2_container = QWidget( )
        plot2_vbox = QVBoxLayout(plot2_container)

        fig2 = Figure(figsize = (5, 5), dpi = 100)
        canvas2 = FigureCanvas(fig2)
        toolbar2 = NavigationToolbar(canvas2, plot2_container)

        ax2 = fig2.add_subplot(111)
        im2 = ax2.pcolormesh(probe, time, best_fit, cmap = 'seismic', shading = 'auto', norm = norm)
        ax2.set_title('Fitted Model')
        ax2.set_xlabel(f'{self.x_axis_label} ({formatted_x_unit})')
        fig2.colorbar(im2, ax = ax2, label = 'Signal')
        fig2.tight_layout( )

        plot2_vbox.addWidget(toolbar2)
        plot2_vbox.addWidget(canvas2)
        self.plot_layout.addWidget(plot2_container)

        plot3_container = QWidget( )
        plot3_vbox = QVBoxLayout(plot3_container)

        residuals = data_interp - best_fit
        fig3 = Figure(figsize = (5, 5), dpi = 100)
        canvas3 = FigureCanvas(fig3)
        toolbar3 = NavigationToolbar(canvas3, plot3_container)

        ax3 = fig3.add_subplot(111)

        res_max_abs = np.abs(residuals).max( )
        res_norm = plt.Normalize(vmin = -res_max_abs, vmax = res_max_abs)

        im3 = ax3.pcolormesh(probe, time, residuals, cmap = 'seismic', shading = 'auto', norm = res_norm)
        ax3.set_title('Residuals (Data - Fit)')
        ax3.set_xlabel(f'{self.x_axis_label} ({formatted_x_unit})')
        fig3.colorbar(im3, ax = ax3, label = 'Residual Signal')
        fig3.tight_layout( )

        plot3_vbox.addWidget(toolbar3)
        plot3_vbox.addWidget(canvas3)
        self.plot_layout.addWidget(plot3_container)

    def handle_error(self, message):
        self.run_button.setDisabled(False)
        self.export_button.setDisabled(True)
        QMessageBox.critical(self, "PFID Analysis Error", message)
        self.results_text_edit.setText(f"ERROR: {message}")

    def export_fit_results(self):
        if self.results_data is None:
            QMessageBox.warning(self, "No Data", "Please run the PFID fit analysis before exporting.")
            return

        (best_fit, A_final, r_squared, probe, time, data_interp, fit_report_string, result,
         num_components) = self.results_data

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save PFID Fit Data", "", "CSV Files (*.csv)"
        )

        if not file_path:
            return

        dir_name = os.path.dirname(file_path)
        base_name = os.path.basename(file_path).split('.') [ 0 ]

        try:
            report_path = os.path.join(dir_name, f"{base_name}_PFID_report.txt")
            with open(report_path, 'w', encoding = 'utf-8') as f:
                f.write(fit_report_string)

            header = [ "Time (abs_ps)" ] + [ f"{p:.2f} ({self.x_axis_unit})" for p in probe ]

            data_to_export = np.hstack((time [ :, np.newaxis ], data_interp))
            fit_to_export = np.hstack((time [ :, np.newaxis ], best_fit))

            data_df = pd.DataFrame(data_to_export, columns = header)
            fit_df = pd.DataFrame(fit_to_export, columns = header)

            data_df.to_csv(os.path.join(dir_name, f"{base_name}_PFID_interpolated_data.csv"), index = False,
                           encoding = 'utf-8')
            fit_df.to_csv(os.path.join(dir_name, f"{base_name}_PFID_fitted_data.csv"), index = False,
                          encoding = 'utf-8')

            QMessageBox.information(
                self, "Export Successful",
                f"PFID fit results exported to multiple files in the directory:\n{dir_name}"
            )

        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"An error occurred during data export: {e}")

    def get_fitter_state(self):
        return {
            'time_min_input': self.time_min_input.text( ),
            'time_max_input': self.time_max_input.text( ),
            'probe_min_input': self.probe_min_input.text( ),
            'probe_max_input': self.probe_max_input.text( ),
            'num_components_input': self.num_components_input.text( ),
            'T2_guess_input': self.T2_guess_input.text( ),
            'nu10_guess_input': self.nu10_guess_input.text( ),
            'nu21_guess_input': self.nu21_guess_input.text( ),
            'r_guess_input': self.r_guess_input.text( ),
            'interp_method': self.interp_method_combo.currentText( ),
            'num_interp_points': self.interp_points_input.value( ),
            'type': 'PFIDFitterApp', 'tab_title': self.objectName( )
        }

    def set_fitter_state(self, state):
        self.setObjectName(state.get('tab_title', 'PFID Fit'))
        self.setWindowTitle(state.get('tab_title', 'PFID Fit'))

        self.time_min_input.setText(state.get('time_min_input', "-15.0"))
        self.time_max_input.setText(state.get('time_max_input', "-0.5"))
        self.probe_min_input.setText(state.get('probe_min_input', "1775"))
        self.probe_max_input.setText(state.get('probe_max_input', "1850"))

        self.num_components_input.setText(state.get('num_components_input', "1"))
        self.T2_guess_input.setText(state.get('T2_guess_input', "2.0"))
        self.nu10_guess_input.setText(state.get('nu10_guess_input', "1810"))
        self.nu21_guess_input.setText(state.get('nu21_guess_input', "1785"))
        self.r_guess_input.setText(state.get('r_guess_input', "0.2"))
        self.interp_method_combo.setCurrentText(state.get('interp_method', "cubic"))
        self.interp_points_input.setValue(state.get('num_interp_points', 100))

        self.results_data = None
        self.results_text_edit.clear( )
        QMessageBox.information(self, "State Loaded",
                                "All PFID fit settings have been restored. Please click 'Run PFID Fit' to re-calculate and display the figures and report.")


class GlobalFitApp(QWidget):
    def __init__(self, x_axis_data=None, y_axis_data=None, two_d_spectrum_data=None, parent=None,
                 x_axis_label='Probe wavenumber', y_axis_label='Time', x_axis_unit='cm\u207B\u00B9', y_axis_unit='ps',
                 font_size=12):
        super( ).__init__(parent)
        self.setWindowTitle("Global Fitting Analysis")
        self.setGeometry(100, 100, 1920, 1080)
        self.x_axis_data = x_axis_data
        self.y_axis_data = y_axis_data
        self.two_d_spectrum_data = two_d_spectrum_data
        self.x_axis_label = x_axis_label
        self.y_axis_label = y_axis_label
        self.x_axis_unit = x_axis_unit
        self.y_axis_unit = y_axis_unit
        self.font_size = font_size
        self.setObjectName("Global Fit")

        self.initUI( )
        self.worker_thread = None

        self.best_fit_data = None
        self.A_final_data = None
        self.best_taus_data = None
        self.r_squared_data = None
        self.probe_data = None
        self.time_data = None
        self.probes_to_plot_data = None
        self.best_t0_data = None
        self.best_fwhm_data = None
        self.best_offset_data = None
        self.use_convolved_model_data = None
        self.initial_tau_guesses_data = None
        self.fit_report_string_data = ""
        self.data_sliced_for_export = None

        self.das_interp_method = "None"
        self.das_interp_multiplier = 1

    def initUI(self):
        main_layout = QGridLayout(self)

        formatted_x_unit = _format_unit_for_display(self.x_axis_unit)
        formatted_y_unit = _format_unit_for_display(self.y_axis_unit)

        input_panel = QWidget( )
        input_layout = QVBoxLayout(input_panel)
        grid_params = QGridLayout( )
        grid_params.addWidget(QLabel(f"{self.y_axis_label} Min ({formatted_y_unit}):"), 0, 0)
        self.time_min_input = QLineEdit("1.0")
        grid_params.addWidget(self.time_min_input, 0, 1)
        grid_params.addWidget(QLabel(f"{self.y_axis_label} Max ({formatted_y_unit}):"), 0, 2)
        self.time_max_input = QLineEdit("80.0")
        grid_params.addWidget(self.time_max_input, 0, 3)
        grid_params.addWidget(QLabel(f"{self.x_axis_label} Min ({formatted_x_unit}):"), 1, 0)
        self.probe_min_input = QLineEdit("1900")
        grid_params.addWidget(self.probe_min_input, 1, 1)
        grid_params.addWidget(QLabel(f"{self.x_axis_label} Max ({formatted_x_unit}):"), 1, 2)
        self.probe_max_input = QLineEdit("2100")
        grid_params.addWidget(self.probe_max_input, 1, 3)
        input_layout.addLayout(grid_params)

        grid_model = QGridLayout( )
        grid_model.addWidget(QLabel("Model Options:"), 0, 0)
        self.convolved_checkbox = QCheckBox("Use Gaussian convoluted model")
        self.convolved_checkbox.setChecked(False)
        grid_model.addWidget(self.convolved_checkbox, 0, 1, 1, 2)

        grid_model.addWidget(QLabel("Number of Components:"), 1, 0)
        self.num_components_input = QLineEdit("2")
        grid_model.addWidget(self.num_components_input, 1, 1, 1, 2)
        input_layout.addLayout(grid_model)

        self.t0_label = QLabel(f"t₀ ({formatted_y_unit}):")
        grid_model.addWidget(self.t0_label, 2, 0)
        self.t0_input = QLineEdit("0")
        grid_model.addWidget(self.t0_input, 2, 1)
        self.fix_t0_checkbox = QCheckBox("Fix")
        grid_model.addWidget(self.fix_t0_checkbox, 2, 2)

        self.fwhm_label = QLabel(f"FWHM ({formatted_y_unit}):")
        grid_model.addWidget(self.fwhm_label, 3, 0)
        self.fwhm_input = QLineEdit("0.15")
        grid_model.addWidget(self.fwhm_input, 3, 1)
        self.fix_fwhm_checkbox = QCheckBox("Fix")
        grid_model.addWidget(self.fix_fwhm_checkbox, 3, 2)

        self.convolved_checkbox.stateChanged.connect(self.toggle_convolved_options)
        self.toggle_convolved_options( )

        grid_guesses = QGridLayout( )
        grid_guesses.addWidget(QLabel("Initial Guesses:"), 0, 0)
        self.svd_checkbox = QCheckBox("Use SVD for inital guess")
        self.svd_checkbox.setChecked(False)
        grid_guesses.addWidget(self.svd_checkbox, 0, 1)
        grid_guesses.addWidget(QLabel("Manual τ (comma-separated):"), 1, 0)
        self.manual_guess_input = QLineEdit("5, 20")
        self.manual_guess_input.setDisabled(self.svd_checkbox.isChecked( ))
        self.svd_checkbox.stateChanged.connect(
            lambda: self.manual_guess_input.setDisabled(self.svd_checkbox.isChecked( )))
        grid_guesses.addWidget(self.manual_guess_input, 1, 1)

        grid_guesses.addWidget(QLabel("Fixed Long τ (ps):"), 2, 0)
        self.fixed_long_tau_input = QLineEdit("1000")
        grid_guesses.addWidget(self.fixed_long_tau_input, 2, 1)

        input_layout.addLayout(grid_guesses)

        das_interp_group = QGroupBox("DAS Interpolation")
        das_interp_layout = QHBoxLayout( )

        das_method_label = QLabel("Method:")
        self.das_interp_method_combo = QComboBox( )
        self.das_interp_method_combo.addItems([ "None", "linear", "cubic" ])
        self.das_interp_method_combo.setCurrentText("cubic")

        das_multiplier_label = QLabel("Multiplier:")
        self.das_interp_multiplier_combo = QComboBox( )
        self.das_interp_multiplier_combo.addItems([ "x1", "x2", "x3", "x5" ])
        self.das_interp_multiplier_combo.setCurrentText("x2")

        das_interp_layout.addWidget(das_method_label)
        das_interp_layout.addWidget(self.das_interp_method_combo)
        das_interp_layout.addWidget(das_multiplier_label)
        das_interp_layout.addWidget(self.das_interp_multiplier_combo)
        das_interp_layout.addStretch(1)

        das_interp_group.setLayout(das_interp_layout)
        input_layout.addWidget(das_interp_group)

        grid_plots = QGridLayout( )
        grid_plots.addWidget(QLabel(f"{self.x_axis_label} to Plot (comma-separated):"), 0, 0)
        self.probes_input = QLineEdit("1970, 2015, 1950, 1940")
        grid_plots.addWidget(self.probes_input, 0, 1)

        grid_plots.addWidget(QLabel("Plot Font Size:"), 1, 0)
        self.plot_font_size_input = QSpinBox( )
        self.plot_font_size_input.setRange(8, 24)
        self.plot_font_size_input.setValue(self.font_size)
        grid_plots.addWidget(self.plot_font_size_input, 1, 1)

        input_layout.addLayout(grid_plots)

        button_layout = QHBoxLayout( )
        self.run_button = QPushButton("Run Analysis")
        self.run_button.clicked.connect(self.run_analysis)
        button_layout.addWidget(self.run_button)

        self.export_button = QPushButton("Export Data")
        self.export_button.clicked.connect(self.export_fit_results)
        button_layout.addWidget(self.export_button)

        input_layout.addLayout(button_layout)

        main_layout.addWidget(input_panel, 0, 0, 1, 1, Qt.AlignTop)

        results_panel = QWidget( )
        results_layout = QVBoxLayout(results_panel)
        results_layout.addWidget(QLabel("Analysis Results:"))
        self.results_text_edit = QTextEdit( )
        self.results_text_edit.setReadOnly(True)
        results_layout.addWidget(self.results_text_edit)
        self.results_text_edit.setMinimumHeight(300)

        main_layout.addWidget(results_panel, 0, 1, 1, 1, Qt.AlignTop)

        self.plot_container = QWidget( )
        self.plot_layout = QHBoxLayout(self.plot_container)
        main_layout.addWidget(self.plot_container, 1, 0, 1, 2)
        main_layout.setRowStretch(1, 10)
        self.setLayout(main_layout)

    def toggle_convolved_options(self):
        enabled = self.convolved_checkbox.isChecked( )
        self.t0_label.setEnabled(enabled)
        self.t0_input.setEnabled(enabled)
        self.fix_t0_checkbox.setEnabled(enabled)
        self.fwhm_label.setEnabled(enabled)
        self.fwhm_input.setEnabled(enabled)
        self.fix_fwhm_checkbox.setEnabled(enabled)

    def export_fit_results(self):
        if self.best_fit_data is None or self.data_sliced_for_export is None:
            QMessageBox.warning(self, "No Data", "Please run the global fit analysis before exporting.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Global Fit Data", "", "CSV Files (*.csv)"
        )

        if not file_path:
            return

        dir_name = os.path.dirname(file_path)
        base_name = os.path.basename(file_path).split('.') [ 0 ]

        try:
            das_df = pd.DataFrame(
                self.A_final_data.T,
                index = self.probe_data,
                columns = [ f"τ = {tau:.2f}" for tau in self.best_taus_data ]
            )
            das_df.to_csv(os.path.join(dir_name, f"{base_name}_DAS.csv"), encoding = 'utf-8')

            fitted_data_df = pd.DataFrame(
                self.best_fit_data,
                index = self.time_data,
                columns = self.probe_data
            )
            fitted_data_df.to_csv(os.path.join(dir_name, f"{base_name}_2D_fitted_data.csv"), encoding = 'utf-8')

            traces_df = pd.DataFrame({
                "Time": self.time_data,
            })
            for probe_val in self.probes_to_plot_data:
                probe_idx = np.argmin(np.abs(self.probe_data - probe_val))

                traces_df [ f"Trace at {probe_val}" ] = self.data_sliced_for_export [ :, probe_idx ]

                traces_df [ f"{probe_val} (Fit)" ] = self.best_fit_data [ :, probe_idx ]

            traces_df.to_csv(os.path.join(dir_name, f"{base_name}_time_traces.csv"), encoding = 'utf-8')

            QMessageBox.information(
                self, "Export Successful",
                f"Fit results exported to multiple files in the directory:\n{dir_name}"
            )

        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"An error occurred during data export: {e}")

    def get_state(self):
        state = {
            'type': 'GlobalFitApp',
            'tab_title': self.objectName( ),
            'time_min_input': self.time_min_input.text( ),
            'time_max_input': self.time_max_input.text( ),
            'probe_min_input': self.probe_min_input.text( ),
            'probe_max_input': self.probe_max_input.text( ),
            'manual_num_components': self.num_components_input.text( ),
            'probes_input': self.probes_input.text( ),
            'use_convolved_model': self.convolved_checkbox.isChecked( ),
            'use_svd_initial_guess': self.svd_checkbox.isChecked( ),
            'manual_tau_guesses': self.manual_guess_input.text( ),
            'fixed_long_tau': self.fixed_long_tau_input.text( ),
            't0_input': self.t0_input.text( ),
            'fwhm_input': self.fwhm_input.text( ),
            'fix_t0': self.fix_t0_checkbox.isChecked( ),
            'fix_fwhm': self.fix_fwhm_checkbox.isChecked( ),
            'x_axis_label': self.x_axis_label,
            'y_axis_label': self.y_axis_label,
            'x_axis_unit': self.x_axis_unit,
            'y_axis_unit': self.y_axis_unit,
            'font_size': self.plot_font_size_input.value( ),
            'x_axis_data': self.x_axis_data.tolist( ) if self.x_axis_data is not None else None,
            'y_axis_data': self.y_axis_data.tolist( ) if self.y_axis_data is not None else None,
            'two_d_spectrum_data': self.two_d_spectrum_data.tolist( ) if self.two_d_spectrum_data is not None else None,
            'das_interp_method': self.das_interp_method_combo.currentText( ),
            'das_interp_multiplier': int(self.das_interp_multiplier_combo.currentText( ).replace('x', '')),
            'results_text': self.results_text_edit.toPlainText( ),
        }
        return state

    def set_fitter_state(self, state):
        self.setObjectName(state.get('tab_title', 'Global Fit'))
        self.setWindowTitle(state.get('tab_title', 'Global Fit'))
        self.x_axis_label = state.get('x_axis_label', 'Probe wavenumber')
        self.y_axis_label = state.get('y_axis_label', 'Time')
        self.x_axis_unit = state.get('x_axis_unit', 'cm\u207B\u00B9')
        self.y_axis_unit = state.get('y_axis_unit', 'ps')
        self.x_axis_data = np.array(state.get('x_axis_data')) if state.get('x_axis_data') is not None else None
        self.y_axis_data = np.array(state.get('y_axis_data')) if state.get('y_axis_data') is not None else None
        self.two_d_spectrum_data = np.array(state.get('two_d_spectrum_data')) if state.get(
            'two_d_spectrum_data') is not None else None

        self.time_min_input.setText(state.get('time_min_input', "1.0"))
        self.time_max_input.setText(state.get('time_max_input', "80.0"))
        self.probe_min_input.setText(state.get('probe_min_input', "1900"))
        self.probe_max_input.setText(state.get('probe_max_input', "2100"))
        self.num_components_input.setText(state.get('manual_num_components', "3"))
        self.probes_input.setText(state.get('probes_input', "1970, 2015, 1950, 1940"))
        self.convolved_checkbox.setChecked(state.get('use_convolved_model', True))
        self.svd_checkbox.setChecked(state.get('use_svd_initial_guess', False))
        self.manual_guess_input.setText(state.get('manual_tau_guesses', "5, 20"))
        self.manual_guess_input.setDisabled(self.svd_checkbox.isChecked( ))
        self.fixed_long_tau_input.setText(state.get('fixed_long_tau', '1000'))
        self.plot_font_size_input.setValue(state.get('font_size', 12))
        self.results_text_edit.setText(state.get('results_text', ''))

        self.t0_input.setText(state.get('t0_input', '0'))
        self.fwhm_input.setText(state.get('fwhm_input', '0.15'))
        self.fix_t0_checkbox.setChecked(state.get('fix_t0', False))
        self.fix_fwhm_checkbox.setChecked(state.get('fix_fwhm', False))

        self.das_interp_method = state.get('das_interp_method', 'None')
        self.das_interp_multiplier = state.get('das_interp_multiplier', 1)
        self.das_interp_method_combo.setCurrentText(self.das_interp_method)
        self.das_interp_multiplier_combo.setCurrentText(f"x{self.das_interp_multiplier}")

        self.toggle_convolved_options( )

        if self.two_d_spectrum_data is not None:
            self.run_analysis( )

    def run_analysis(self):
        try:
            params = {
                "time_min": float(self.time_min_input.text( )),
                "time_max": float(self.time_max_input.text( )),
                "probe_min": float(self.probe_min_input.text( )),
                "probe_max": float(self.probe_max_input.text( )),
                "manual_num_components": int(self.num_components_input.text( )),
                "probes_to_plot": [ float(p.strip( )) for p in self.probes_input.text( ).split(',') if p.strip( ) ],
                "use_convolved_model": self.convolved_checkbox.isChecked( ),
                "use_svd_initial_guess": self.svd_checkbox.isChecked( ),
                "manual_tau_guesses": [ float(t.strip( )) for t in self.manual_guess_input.text( ).split(
                    ',') if t.strip( ) ] if not self.svd_checkbox.isChecked( ) else None,
                "manual_t0_guess": float(self.t0_input.text( )) if self.convolved_checkbox.isChecked( ) else None,
                "manual_fwhm_guess": float(self.fwhm_input.text( )) if self.convolved_checkbox.isChecked( ) else None,
                "fix_t0": self.fix_t0_checkbox.isChecked( ) if self.convolved_checkbox.isChecked( ) else False,
                "fix_fwhm": self.fix_fwhm_checkbox.isChecked( ) if self.convolved_checkbox.isChecked( ) else False,
                "x_axis": self.x_axis_data,
                "y_axis": self.y_axis_data,
                "two_d_spectrum": self.two_d_spectrum_data,
                "fixed_tau_indices": [ ],
                "fixed_long_tau": float(self.fixed_long_tau_input.text( ))
            }

            self.results_text_edit.clear( )
            self.run_button.setDisabled(True)
            self.worker_thread = AnalysisWorker(params)
            self.worker_thread.finished.connect(self.plot_results)
            self.worker_thread.error.connect(self.handle_error)
            self.worker_thread.start( )

        except ValueError as ve:
            QMessageBox.critical(self, "Input Error",
                                 f"Invalid number format in one of the input fields. Please check your values and ensure they contain only numbers and a decimal point. Error: {ve}")
            self.run_button.setDisabled(False)
        except Exception as e:
            QMessageBox.critical(self, "Input Error", f"Please check your input values. Error: {e}")
            self.run_button.setDisabled(False)

    def plot_results(self, results):
        self.run_button.setDisabled(False)
        if not results:
            return

        (best_fit, A_final, best_taus, r_squared, probe, time, data_sliced, probes_to_plot, best_t0,
         best_fwhm, best_offset, use_convolved_model, initial_tau_guesses_for_print, fit_report_string,
         warnings) = results

        self.best_fit_data = best_fit
        self.A_final_data = A_final
        self.best_taus_data = best_taus
        self.r_squared_data = r_squared
        self.probe_data = probe
        self.time_data = time
        self.data_sliced_for_export = data_sliced
        self.probes_to_plot_data = probes_to_plot
        self.best_t0_data = best_t0
        self.best_fwhm_data = best_fwhm
        self.best_offset_data = best_offset
        self.use_convolved_model_data = use_convolved_model
        self.initial_tau_guesses_data = initial_tau_guesses_for_print
        self.fit_report_string_data = fit_report_string

        output_text = "--- Global Fit Results ---\n\n"
        output_text += f"Overall R-squared value: {r_squared:.4f}\n"
        output_text += f"Initial Tau Guesses: {initial_tau_guesses_for_print}\n\n"
        output_text += "--- Full lmfit Report ---\n"
        output_text += fit_report_string

        if warnings:
            output_text += "\n--- Warnings ---\n"
            for warning in warnings:
                output_text += f"- {warning}\n"

        self.results_text_edit.setText(output_text)

        for i in reversed(range(self.plot_layout.count( ))):
            self.plot_layout.itemAt(i).widget( ).setParent(None)

        plot1_widget = QWidget( )
        plot1_layout = QVBoxLayout(plot1_widget)
        plot2_widget = QWidget( )
        plot2_layout = QVBoxLayout(plot2_widget)
        plot3_widget = QWidget( )
        plot3_layout = QVBoxLayout(plot3_widget)

        font_size = self.plot_font_size_input.value( )
        plt.rcParams.update({ 'font.size': font_size })

        formatted_x_unit = self.x_axis_unit
        formatted_y_unit = self.y_axis_unit

        fig1 = Figure(figsize = (6, 6), dpi = 100)
        canvas1 = FigureCanvas(fig1)
        ax1 = fig1.add_subplot(111)
        toolbar1 = NavigationToolbar(canvas1, plot1_widget)
        plot1_layout.addWidget(toolbar1)
        plot1_layout.addWidget(canvas1)

        tau_string = ", ".join([ f'τ{i + 1} = {t:.2f} {formatted_y_unit}' for i, t in enumerate(best_taus) ])
        norm = TwoSlopeNorm(vmin = best_fit.min( ), vcenter = 0, vmax = best_fit.max( ))
        im1 = ax1.contourf(probe, time, best_fit, cmap = 'seismic', norm = norm, levels = 100)
        if use_convolved_model:
            ax1.set_title(
                f'Fitted data')
        else:
            ax1.set_title(f'Fit data')
        ax1.set_xlabel(f'{self.x_axis_label} ({formatted_x_unit})')
        ax1.set_ylabel(f'{self.y_axis_label} ({formatted_y_unit})')
        fig1.tight_layout( )
        ax1.minorticks_on( )
        canvas1.draw( )

        self.plot_layout.addWidget(plot1_widget)

        fig2 = Figure(figsize = (6, 6), dpi = 100)
        canvas2 = FigureCanvas(fig2)
        ax2 = fig2.add_subplot(111)
        toolbar2 = NavigationToolbar(canvas2, plot2_widget)
        plot2_layout.addWidget(toolbar2)
        plot2_layout.addWidget(canvas2)

        das_interp_method = self.das_interp_method_combo.currentText( )
        das_interp_multiplier = int(self.das_interp_multiplier_combo.currentText( ).replace('x', ''))

        if das_interp_method != "None" and das_interp_multiplier > 1:
            try:
                probe_interp = np.linspace(probe.min( ), probe.max( ), len(probe) * das_interp_multiplier)
                for i in range(A_final.shape [ 0 ]):
                    f_interp = interp1d(probe, A_final [ i, : ], kind = das_interp_method,
                                        fill_value = "extrapolate")
                    das_interp = f_interp(probe_interp)
                    ax2.plot(probe_interp, das_interp, linewidth = 2,
                             label = f'τ = {best_taus [ i ]:.2f} {formatted_y_unit}')
            except Exception as e:
                QMessageBox.warning(self, "Interpolation Error", f"DAS interpolation failed: {e}. Plotting raw DAS.")
                for i in range(A_final.shape [ 0 ]):
                    ax2.plot(probe, A_final [ i, : ], linewidth = 2,
                             label = f'τ = {best_taus [ i ]:.2f} {formatted_y_unit}')
        else:
            for i in range(A_final.shape [ 0 ]):
                ax2.plot(probe, A_final [ i, : ], linewidth = 2,
                         label = f'τ = {best_taus [ i ]:.2f} {formatted_y_unit}')

        ax2.set_title('DAS spectra')
        ax2.set_xlabel(f'{self.x_axis_label} ({formatted_x_unit})')
        ax2.set_ylabel('Absorbance change (mOD)')
        ax2.legend( )
        ax2.minorticks_on( )
        ax2.grid(True)
        fig2.tight_layout( )
        canvas2.draw( )

        self.plot_layout.addWidget(plot2_widget)

        fig3 = Figure(figsize = (6, 6), dpi = 100)
        canvas3 = FigureCanvas(fig3)
        ax3 = fig3.add_subplot(111)
        toolbar3 = NavigationToolbar(canvas3, plot3_widget)
        plot3_layout.addWidget(toolbar3)
        plot3_layout.addWidget(canvas3)

        for p_val in probes_to_plot:
            idx = find(probe, p_val)
            ax3.plot(time, data_sliced [ :, idx ], '-', linewidth = 2, label = f'{p_val} {formatted_x_unit}')
            ax3.plot(time, best_fit [ :, idx ], '-', linewidth = 1, color = 'k')
        ax3.set_title('Time Traces with Fits')
        ax3.set_xlabel(f'{self.y_axis_label} ({formatted_y_unit})')
        ax3.set_ylabel('Signal (mOD)')
        ax3.legend( )
        ax3.minorticks_on( )
        ax3.grid(True)
        fig3.tight_layout( )
        canvas3.draw( )

        self.plot_layout.addWidget(plot3_widget)

    def handle_error(self, message):
        self.run_button.setDisabled(False)
        QMessageBox.critical(self, "Analysis Error", message)


def gaussian(x, amp, pos, fwhm):
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return amp * np.exp(-(x - pos) ** 2 / (2 * sigma ** 2))


def multi_gaussian(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        amp, pos, fwhm = params [ i:i + 3 ]
        y += gaussian(x, amp, pos, fwhm)
    return y


def lorentzian(x, amplitude, mean, fwhm):
    gamma = fwhm / 2.0
    return amplitude * (gamma ** 2 / ((x - mean) ** 2 + gamma ** 2))


def multi_lorentzian(x, *params):
    y_sum = np.zeros_like(x, dtype = float)
    for i in range(0, len(params), 3):
        amplitude, mean, fwhm = params [ i:i + 3 ]
        y_sum += lorentzian(x, amplitude, mean, fwhm)
    return y_sum


def exponential(x, amp, tau):
    return amp * np.exp(-x / tau)


def multi_exponential(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 2):
        amp, tau = params [ i:i + 2 ]
        y += exponential(x, amp, tau)
    return y


class LineThicknessDialog(QDialog):
    def __init__(self, parent=None, initial_thickness=2):
        super( ).__init__(parent)
        self.setWindowTitle("Change Line Thickness")
        self.setGeometry(300, 300, 250, 100)
        layout = QVBoxLayout( )
        h_layout = QHBoxLayout( )
        h_layout.addWidget(QLabel("Line Thickness (px):"))
        self.thickness_input = QSpinBox(self)
        self.thickness_input.setRange(1, 10)
        self.thickness_input.setValue(initial_thickness)
        h_layout.addWidget(self.thickness_input)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
        self.setLayout(layout)

    def get_thickness(self):
        return self.thickness_input.value( )


class EditNamesDialog(QDialog):
    def __init__(self, parent=None, current_labels=None):
        super( ).__init__(parent)
        self.setWindowTitle("Edit Axis Names")
        self.setGeometry(200, 200, 1000, 550)
        self.labels = current_labels if current_labels else {
            'signal_bottom': 'Probe wavenumber [cm\u207B\u00B9]',
            'signal_left': 'Time [ps]',
            'x_slice_bottom': 'Time [ps]',
            'x_slice_left': 'ΔOD',
            'y_slice_bottom': 'Probe wavenumber [cm\u207B\u00B9]',
            'y_slice_left': 'ΔOD'
        }
        self.init_ui( )

    def init_ui(self):
        layout = QVBoxLayout( )
        form_layout = QGridLayout( )
        dialog_font = QFont("Times New Roman")
        dialog_font.setPointSize(12)

        label_signal_bottom = QLabel("2D Plot - X-axis:")
        label_signal_bottom.setFont(dialog_font)
        form_layout.addWidget(label_signal_bottom, 0, 0)
        self.signal_bottom_input = QLineEdit(self)
        self.signal_bottom_input.setText(self.labels [ 'signal_bottom' ])
        self.signal_bottom_input.setFont(dialog_font)
        form_layout.addWidget(self.signal_bottom_input, 0, 1)

        label_signal_left = QLabel("2D Plot - Y-axis:")
        label_signal_left.setFont(dialog_font)
        form_layout.addWidget(label_signal_left, 1, 0)
        self.signal_left_input = QLineEdit(self)
        self.signal_left_input.setText(self.labels [ 'signal_left' ])
        self.signal_left_input.setFont(dialog_font)
        form_layout.addWidget(self.signal_left_input, 1, 1)

        label_x_slice_bottom = QLabel("X-Slice Plot - X-axis:")
        label_x_slice_bottom.setFont(dialog_font)
        form_layout.addWidget(label_x_slice_bottom, 2, 0)
        self.x_slice_bottom_input = QLineEdit(self)
        self.x_slice_bottom_input.setText(self.labels [ 'x_slice_bottom' ])
        self.x_slice_bottom_input.setFont(dialog_font)
        form_layout.addWidget(self.x_slice_bottom_input, 2, 1)

        label_x_slice_left = QLabel("X-Slice Plot - Y-axis:")
        label_x_slice_left.setFont(dialog_font)
        form_layout.addWidget(label_x_slice_left, 3, 0)
        self.x_slice_left_input = QLineEdit(self)
        self.x_slice_left_input.setText(self.labels [ 'x_slice_left' ])
        self.x_slice_left_input.setFont(dialog_font)
        form_layout.addWidget(self.x_slice_left_input, 3, 1)

        label_y_slice_bottom = QLabel("Y-Slice Plot - X-axis:")
        label_y_slice_bottom.setFont(dialog_font)
        form_layout.addWidget(label_y_slice_bottom, 4, 0)
        self.y_slice_bottom_input = QLineEdit(self)
        self.y_slice_bottom_input.setText(self.labels [ 'y_slice_bottom' ])
        self.y_slice_bottom_input.setFont(dialog_font)
        form_layout.addWidget(self.y_slice_bottom_input, 4, 1)

        label_y_slice_left = QLabel("Y-Slice Plot - Y-axis:")
        label_y_slice_left.setFont(dialog_font)
        form_layout.addWidget(label_y_slice_left, 5, 0)
        self.y_slice_left_input = QLineEdit(self)
        self.y_slice_left_input.setText(self.labels [ 'y_slice_left' ])
        self.y_slice_left_input.setFont(dialog_font)
        form_layout.addWidget(self.y_slice_left_input, 5, 1)

        layout.addLayout(form_layout)
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)
        self.setLayout(layout)

    def get_names(self):
        return {
            'signal_bottom': self.signal_bottom_input.text( ),
            'signal_left': self.signal_left_input.text( ),
            'x_slice_bottom': self.x_slice_bottom_input.text( ),
            'x_slice_left': self.x_slice_left_input.text( ),
            'y_slice_bottom': self.y_slice_bottom_input.text( ),
            'y_slice_left': self.y_slice_left_input.text( )
        }


class ExponentialFitterApp(QWidget):
    def __init__(self, parent=None, x_data=None, y_data=None, xlabel="X-axis",
                 ylabel="Y-axis", slice_axis_name="", slice_value=None,
                 slice_unit="", is_spline_corrected=False):
        super( ).__init__(parent)
        self.x_data = x_data if x_data is not None else np.array([ ])
        self.y_data = y_data if y_data is not None else np.array([ ])
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.slice_axis_name = slice_axis_name
        self.slice_value = slice_value
        self.slice_unit = slice_unit
        self.is_spline_corrected = is_spline_corrected
        self.is_guessing_mode_active = False
        title_parts = [ "Exponential:" ]
        if self.slice_axis_name and self.slice_value is not None:
            title_parts.append(f"{self.slice_axis_name} = {self.slice_value:.1f}{self.slice_unit}")
        self.setWindowTitle(" ".join(title_parts))
        self.setObjectName(" ".join(title_parts))
        self.init_ui( )
        self.init_fitter_variables( )
        self.update_plot( )

    def init_ui(self):
        self.main_layout = QVBoxLayout(self)
        self.fig, self.ax = plt.subplots(figsize = (15, 10))
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.params_text_edit = QTextEdit( )
        self.params_text_edit.setReadOnly(True)
        self.params_text_edit.setMinimumHeight(100)
        self.params_text_edit.setStyleSheet("font-family: Consolas; font-size: 10pt;")
        self.start_guess_button = QPushButton("Start Initial Guess")
        self.fit_button = QPushButton("Fit Exponential")
        self.export_button = QPushButton("Export Fit Data")
        self.export_button.setDisabled(True)
        self.clear_button = QPushButton("Clear Guesses")
        self.close_button = QPushButton("Close Tab")
        self.info_label = QLabel("Use toolbar to zoom/pan. Click 'Start Initial Guess' to define components.")
        self.info_label.setWordWrap(True)
        self.splitter = QSplitter(Qt.Vertical)
        self.splitter.addWidget(self.canvas)
        self.splitter.addWidget(self.params_text_edit)
        self.splitter.setSizes([ 700, 300 ])
        self.control_layout = QHBoxLayout( )
        self.control_layout.addWidget(self.start_guess_button)
        self.control_layout.addWidget(self.fit_button)
        self.control_layout.addWidget(self.export_button)
        self.control_layout.addWidget(self.clear_button)
        self.control_layout.addWidget(self.info_label)
        self.control_layout.addWidget(self.close_button)
        self.control_layout.setSpacing(10)
        self.main_layout.addWidget(self.toolbar)
        self.main_layout.addLayout(self.control_layout)
        self.main_layout.addWidget(self.splitter)
        self.start_guess_button.clicked.connect(self._toggle_guessing_mode)
        self.fit_button.clicked.connect(self.on_fit)
        self.export_button.clicked.connect(self.export_fit_data)
        self.clear_button.clicked.connect(self.on_clear_guesses)
        self.close_button.clicked.connect(self._close_tab)
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def init_fitter_variables(self):
        self.current_component = None
        self.start_x = None
        self.fixed_components = [ ]
        self.fitted_params = None
        self.fitted_errors = None

    def _toggle_guessing_mode(self):
        self.is_guessing_mode_active = not self.is_guessing_mode_active
        if self.is_guessing_mode_active:
            self.current_component = None
            self.start_guess_button.setText("Stop Initial Guess")
            self.info_label.setText("Guessing Mode: ON. Click for amplitude, then drag to set decay and click again.")
        else:
            self.start_guess_button.setText("Start Initial Guess")
            self.info_label.setText("Guessing Mode: OFF. Use toolbar to zoom/pan.")
            self.current_component = None
        self.update_plot( )

    def on_click(self, event):
        if event.inaxes != self.ax or not self.is_guessing_mode_active:
            return
        if event.button == 1:
            if self.current_component is None:
                self.current_component = [ event.ydata, 1.0 ]
                self.start_x = event.xdata
                self.info_label.setText("Drag to adjust decay (tau), then click again to fix.")
            else:
                self.fixed_components.append(tuple(self.current_component))
                self.current_component = None
                self.info_label.setText(
                    f"Component {len(self.fixed_components)} fixed. Click for next, or 'Stop Guessing'.")
            self.update_plot( )

    def on_motion(self, event):
        if event.inaxes != self.ax or self.current_component is None or not self.is_guessing_mode_active:
            return
        self.current_component [ 1 ] = max(0.01, abs(event.xdata - self.start_x))
        self.update_plot( )

    def on_fit(self):
        if not self.fixed_components:
            QMessageBox.warning(self, "No Components", "Please fix at least one exponential component before fitting.")
            return
        initial_params = np.array(self.fixed_components).flatten( )
        try:
            self.fitted_params, pcov = curve_fit(multi_exponential, self.x_data, self.y_data, p0 = initial_params,
                                                 bounds = (-np.inf, np.inf))
            self.fitted_errors = np.sqrt(np.diag(pcov))
            self.display_fitted_parameters( )
            self.update_plot( )
            self.info_label.setText("Fitting complete. See fitted parameters below.")
            self.export_button.setDisabled(False)
        except RuntimeError as e:
            QMessageBox.critical(self, "Fitting Error", f"Failed to fit: {e}. Adjust initial guesses and try again.")
            self.info_label.setText("Fitting failed. Adjust guesses.")
            self.export_button.setDisabled(True)

    def on_clear_guesses(self):
        self.init_fitter_variables( )
        self.params_text_edit.clear( )
        self.info_label.setText("All guesses cleared. Click 'Start Initial Guess' to begin.")
        self.export_button.setDisabled(True)
        self.update_plot( )

    def export_fit_data(self):
        if self.fitted_params is None:
            QMessageBox.warning(self, "Export Error", "Please run a successful fit before exporting data.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Exponential Fit Data", "", "CSV Files (*.csv);;Text Files (*.txt)"
        )

        if not file_path:
            return

        try:
            fit_y = multi_exponential(self.x_data, *self.fitted_params)

            data_dict = {
                self.xlabel: self.x_data,
                f"Data ({self.ylabel})": self.y_data,
                f"Fit Total ({self.ylabel})": fit_y
            }

            for i in range(0, len(self.fitted_params), 2):
                amp, tau = self.fitted_params [ i:i + 2 ]
                comp_y = exponential(self.x_data, amp, tau)
                data_dict [ f"Component {i // 2 + 1} (tau={tau:.2g})" ] = comp_y

            df = pd.DataFrame(data_dict)

            report_lines = self.params_text_edit.toPlainText( ).split('\n')
            comment_lines = [ f"# {line}" for line in report_lines if line.strip( ) ]

            with open(file_path, 'w', encoding = 'utf-8') as f:
                f.write('\n'.join(comment_lines) + '\n')
                df.to_csv(f, index = False, lineterminator = '\n')

            QMessageBox.information(self, "Export Successful", f"Fit data exported to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"An error occurred during data export: {e}")

    def display_fitted_parameters(self):
        if self.fitted_params is None:
            self.params_text_edit.clear( )
            return
        output_text = "Fitted Exponential Parameters:\n"
        output_text += "----------------------------------\n"
        for i in range(0, len(self.fitted_params), 2):
            amp, tau = self.fitted_params [ i:i + 2 ]
            amp_err, tau_err = self.fitted_errors [ i:i + 2 ]
            output_text += (f"Component {i // 2 + 1}:\n"
                            f"  Amplitude (A): {amp:.4g} ± {amp_err:.2g}\n"
                            f"  Decay (τ):     {tau:.4g} ± {tau_err:.2g}\n"
                            f"----------------------------------\n")
        self.params_text_edit.setText(output_text)

    def update_plot(self):
        self.ax.clear( )
        self.ax.plot(self.x_data, self.y_data, 'b-', label = 'Data')
        if self.fitted_params is None and (len(self.fixed_components) > 0 or self.current_component):
            y_sum_guess = np.zeros_like(self.x_data)
            for amp, tau in self.fixed_components:
                y_sum_guess += exponential(self.x_data, amp, tau)
            if self.current_component:
                y_sum_guess += exponential(self.x_data, *self.current_component)
            self.ax.plot(self.x_data, y_sum_guess, 'r--', alpha = 0.7, label = 'Initial Guess Sum')
        if self.fitted_params is not None:
            y_fit_total = multi_exponential(self.x_data, *self.fitted_params)
            self.ax.plot(self.x_data, y_fit_total, 'g-', linewidth = 2, label = 'Fitted Sum')
            for i in range(0, len(self.fitted_params), 2):
                params = self.fitted_params [ i:i + 2 ]
                self.ax.plot(self.x_data, exponential(self.x_data, *params), '--', alpha = 0.6,
                             label = f'Component {i // 2 + 1}')
        self.ax.set_xlabel(self.xlabel, fontsize = 18)
        self.ax.set_ylabel(self.ylabel, fontsize = 18)
        self.ax.grid(True)
        self.ax.legend(loc = 'best', fontsize = 18)
        self.ax.tick_params(axis = 'both', which = 'major', labelsize = 18)
        self.canvas.draw( )

    def get_fitter_state(self):
        return {
            'type': 'ExponentialFitterApp',
            'tab_title': self.objectName( ),
            'x_data': self.x_data.tolist( ) if self.x_data is not None else None,
            'y_data': self.y_data.tolist( ) if self.y_data is not None else None,
            'xlabel': self.xlabel,
            'ylabel': self.ylabel,
            'slice_axis_name': self.slice_axis_name,
            'slice_value': self.slice_value,
            'slice_unit': self.slice_unit,
            'fixed_components': self.fixed_components,
            'fitted_params': self.fitted_params.tolist( ) if self.fitted_params is not None else None,
            'fitted_errors': self.fitted_errors.tolist( ) if self.fitted_errors is not None else None,
        }

    def set_fitter_state(self, state):
        self.setObjectName(state.get('tab_title', 'Exponential Fitter'))
        self.setWindowTitle(state.get('tab_title', 'Exponential Fitter'))

        self.x_data = np.array(state.get('x_data')) if state.get('x_data') is not None else np.array([ ])
        self.y_data = np.array(state.get('y_data')) if state.get('y_data') is not None else np.array([ ])
        self.xlabel = state.get('xlabel', "X-axis")
        self.ylabel = state.get('ylabel', "Y-axis")
        self.slice_axis_name = state.get('slice_axis_name', "")
        self.slice_value = state.get('slice_value')
        self.slice_unit = state.get('slice_unit', "")
        self.fixed_components = state.get('fixed_components', [ ])
        self.fitted_params = np.array(state.get('fitted_params')) if state.get('fitted_params') is not None else None
        self.fitted_errors = np.array(state.get('fitted_errors')) if state.get('fitted_errors') is not None else None

        self.export_button.setDisabled(self.fitted_params is None)

        self.display_fitted_parameters( )
        self.update_plot( )

    def _close_tab(self):
        parent_tab_widget = self.parent( )
        if parent_tab_widget and isinstance(parent_tab_widget, QTabWidget):
            tab_index = parent_tab_widget.indexOf(self)
            if tab_index != -1:
                parent_tab_widget.removeTab(tab_index)
        self.deleteLater( )


class GaussianFitterApp(QWidget):
    def __init__(self, parent=None, x_data=None, y_data=None, fitting_function_type="Gaussian", xlabel="X-axis",
                 ylabel="Y-axis",
                 slice_axis_name="", slice_value=None, slice_unit="", is_spline_corrected=False):
        super( ).__init__(parent)
        self.xga = x_data if x_data is not None else np.array([ ])
        self.yga = y_data if y_data is not None else np.array([ ])
        self._original_x_data = np.copy(x_data) if x_data is not None else np.array([ ])
        self._original_y_data = np.copy(y_data) if y_data is not None else np.array([ ])
        self.fitting_function_type = fitting_function_type.capitalize( )
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.slice_axis_name = slice_axis_name
        self.slice_value = slice_value
        self.slice_unit = slice_unit
        self.is_spline_corrected = is_spline_corrected
        self.is_guessing_mode_active = False
        title_parts = [ f"{self.fitting_function_type}:" ]
        if self.slice_axis_name and self.slice_value is not None:
            title_parts.append(f"{self.slice_axis_name} = {self.slice_value:.1f}{self.slice_unit}")
        self.setWindowTitle(" ".join(title_parts))
        self.setObjectName(" ".join(title_parts))
        self.init_ui( )
        self.init_fitter_variables( )
        self.update_plot( )

    def init_ui(self):
        self.main_layout = QVBoxLayout(self)
        self.fig, self.ax = plt.subplots(figsize = (15, 10))
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.params_text_edit = QTextEdit( )
        self.params_text_edit.setReadOnly(True)
        self.params_text_edit.setMinimumHeight(100)
        self.params_text_edit.setStyleSheet("font-family: Consolas; font-size: 10pt;")
        self.start_guess_button = QPushButton("Start Initial Guess")
        self.fit_button = QPushButton(f"Fit {self.fitting_function_type}")
        self.export_button = QPushButton("Export Fit Data")
        self.export_button.setDisabled(True)
        self.clear_button = QPushButton("Clear Guesses")
        self.close_button = QPushButton("Close Tab")
        self.info_label = QLabel(
            f"Use the toolbar for zoom/pan. Click 'Start Initial Guess' to define {self.fitting_function_type} parameters.")
        self.info_label.setWordWrap(True)
        self.splitter = QSplitter(Qt.Vertical)
        self.splitter.addWidget(self.canvas)
        self.splitter.addWidget(self.params_text_edit)
        self.splitter.setSizes([ 700, 300 ])
        self.control_layout = QHBoxLayout( )
        self.control_layout.addWidget(self.start_guess_button)
        self.control_layout.addWidget(self.fit_button)
        self.control_layout.addWidget(self.export_button)
        self.control_layout.addWidget(self.clear_button)
        self.control_layout.addWidget(self.info_label)
        self.control_layout.addWidget(self.close_button)
        self.control_layout.setSpacing(10)
        self.main_layout.addWidget(self.toolbar)
        self.main_layout.addLayout(self.control_layout)
        self.main_layout.addWidget(self.splitter)
        self.start_guess_button.clicked.connect(self._toggle_guessing_mode)
        self.fit_button.clicked.connect(self.on_fit)
        self.export_button.clicked.connect(self.export_fit_data)
        self.clear_button.clicked.connect(self.on_clear_guesses)
        self.close_button.clicked.connect(self._close_tab)
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def init_fitter_variables(self):
        self.amp = None
        self.pos = None
        self.fwhm = None
        self.temp_line = None
        self.fixed_peaks = [ ]
        self.fitted_params = None
        self.fitted_errors = None
        self.original_xlim = self.ax.get_xlim( )
        self.original_ylim = self.ax.get_ylim( )

    def _get_single_function(self):
        if self.fitting_function_type == "Gaussian":
            return gaussian
        elif self.fitting_function_type == "Lorentzian":
            return lorentzian
        else:
            raise ValueError("Invalid fitting function type selected.")

    def _get_multi_function(self):
        if self.fitting_function_type == "Gaussian":
            return multi_gaussian
        elif self.fitting_function_type == "Lorentzian":
            return multi_lorentzian
        else:
            raise ValueError("Invalid fitting function type selected.")

    def _toggle_guessing_mode(self):
        self.is_guessing_mode_active = not self.is_guessing_mode_active
        if self.is_guessing_mode_active:
            self.amp = None
            self.pos = None
            self.fwhm = None
            self.temp_line = None
            self.start_guess_button.setText("Stop Initial Guess")
            self.info_label.setText("Guessing Mode: ON. Click on the plot for peak position and amplitude.")
            self.update_plot( )
        else:
            self.start_guess_button.setText("Start Initial Guess")
            self.info_label.setText(
                "Guessing Mode: OFF. Use the toolbar for zoom/pan. Click 'Start Initial Guess' to define parameters.")
            self.amp = None
            self.pos = None
            self.fwhm = None
            self.temp_line = None
            self.update_plot( )

    def on_click(self, event):
        if event.inaxes != self.ax:
            return
        if not self.is_guessing_mode_active:
            return
        if event.button == 1:
            if self.amp is None:
                self.amp = event.ydata
                self.pos = event.xdata
                self.fwhm = (self.xga [ -1 ] - self.xga [ 0 ]) / 10 if (self.xga [ -1 ] - self.xga [ 0 ]) != 0 else 1.0
                self.start_x = event.xdata
                self.info_label.setText("Drag mouse to adjust FWHM, then click again to fix.")
            else:
                self.fixed_peaks.append((self.amp, self.pos, self.fwhm))
                self.amp = None
                self.pos = None
                self.fwhm = None
                self.temp_line = None
                self.info_label.setText(
                    f"{self.fitting_function_type} {len(self.fixed_peaks)} fixed. Click for next, or 'Stop Initial Guess'.")
            self.update_plot( )
        elif event.button == 3:
            if self.amp is not None:
                self.amp = None
                self.pos = None
                self.fwhm = None
                self.temp_line = None
                self.info_label.setText("Current peak selection cancelled. Click to start new guess.")
                self.update_plot( )

    def on_motion(self, event):
        if event.inaxes != self.ax or self.amp is None or not self.is_guessing_mode_active:
            return
        self.fwhm = 2 * abs(event.xdata - self.start_x)
        if self.fwhm < 0.001:
            self.fwhm = 0.001
        self.update_plot( )

    def apply_interpolation_settings(self, method, multiplier):
        self.on_clear_guesses( )
        if method == "None":
            self.xga = np.copy(self._original_x_data)
            self.yga = np.copy(self._original_y_data)
        else:
            try:
                target_n_points = int(len(self._original_x_data) * multiplier)
                if target_n_points < 2:
                    target_n_points = 2
                x_interp = np.linspace(self._original_x_data.min( ), self._original_x_data.max( ), target_n_points)
                f_interp = interp1d(self._original_x_data, self._original_y_data, kind = method,
                                    fill_value = "extrapolate")
                y_interp = f_interp(x_interp)
                self.xga = x_interp
                self.yga = y_interp
            except ValueError as e:
                QMessageBox.warning(self, "Interpolation Error",
                                    f"Could not apply {method} interpolation: {e}")
                self.xga = np.copy(self._original_x_data)
                self.yga = np.copy(self._original_y_data)
            except Exception as e:
                QMessageBox.warning(self, "Interpolation Error",
                                    f"An unexpected error occurred during interpolation: {e}")
                self.xga = np.copy(self._original_x_data)
                self.yga = np.copy(self._original_y_data)
        self.update_plot( )

    def on_fit(self):
        if not self.fixed_peaks:
            QMessageBox.warning(self, "No Peaks",
                                f"Please fix at least one {self.fitting_function_type} guess before fitting.")
            return
        if self.is_guessing_mode_active:
            self.is_guessing_mode_active = False
            self.start_guess_button.setText("Start Initial Guess")
            self.info_label.setText("Guessing Mode: OFF. Fitting in progress...")
        params = np.array(self.fixed_peaks).flatten( )
        lower_bounds = [ ]
        upper_bounds = [ ]
        for _ in range(len(self.fixed_peaks)):
            lower_bounds.extend([ -np.inf, self.xga.min( ), 0.001 ])
            upper_bounds.extend([ np.inf, self.xga.max( ), np.inf ])
        try:
            self.fitted_params, pcov = curve_fit(self._get_multi_function( ), self.xga, self.yga, p0 = params,
                                                 bounds = (lower_bounds, upper_bounds))
            self.fitted_errors = np.sqrt(np.diag(pcov))
            self.update_plot( )
            self.display_fitted_parameters( )
            self.info_label.setText("Fitting complete. See fitted parameters below.")
            self.export_button.setDisabled(False)
        except RuntimeError as e:
            QMessageBox.critical(self, "Fitting Error",
                                 f"Failed to fit {self.fitting_function_type}: {e}. Adjust peaks and try again.")
            self.info_label.setText("Fitting failed. Adjust guesses and try again.")
            self.export_button.setDisabled(True)
        except ValueError as e:
            QMessageBox.critical(self, "Fitting Error",
                                 f"Fitting input error: {e}. Check data and initial parameters.")
            self.info_label.setText("Fitting failed due to input error. Check data.")
            self.export_button.setDisabled(True)

    def on_clear_guesses(self):
        self.fixed_peaks = [ ]
        self.amp = None
        self.pos = None
        self.fwhm = None
        self.temp_line = None
        self.fitted_params = None
        self.fitted_errors = None
        self.params_text_edit.clear( )
        self.info_label.setText("All peak guesses cleared. Click 'Start Initial Guess' to define new ones.")
        self.export_button.setDisabled(True)
        self.update_plot( )

    def export_fit_data(self):
        if self.fitted_params is None:
            QMessageBox.warning(self, "Export Error", "Please run a successful fit before exporting data.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, f"Save {self.fitting_function_type} Fit Data", "", "CSV Files (*.csv);;Text Files (*.txt)"
        )

        if not file_path:
            return

        try:
            multi_func = self._get_multi_function( )
            single_func = self._get_single_function( )
            fit_y = multi_func(self.xga, *self.fitted_params)

            data_dict = {
                self.xlabel: self.xga,
                f"Data ({self.ylabel})": self.yga,
                f"Fit Total ({self.ylabel})": fit_y
            }

            for i in range(0, len(self.fitted_params), 3):
                amp, pos, fwhm = self.fitted_params [ i:i + 3 ]
                comp_y = single_func(self.xga, amp, pos, fwhm)

                if self.fitting_function_type == "Gaussian":
                    name = f"Component {i // 3 + 1} (Pos={pos:.2f}, FWHM={fwhm:.2f})"
                else:
                    name = f"Component {i // 3 + 1} (Mean={pos:.2f}, FWHM={fwhm:.2f})"

                data_dict [ name ] = comp_y

            df = pd.DataFrame(data_dict)

            report_lines = self.params_text_edit.toPlainText( ).split('\n')
            comment_lines = [ f"# {line}" for line in report_lines if line.strip( ) ]

            with open(file_path, 'w', encoding = 'utf-8') as f:
                f.write('\n'.join(comment_lines) + '\n')
                df.to_csv(f, index = False, lineterminator = '\n')

            QMessageBox.information(self, "Export Successful", f"Fit data exported to:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "Export Error", f"An error occurred during data export: {e}")

    def display_fitted_parameters(self):
        if self.fitted_params is None:
            self.params_text_edit.clear( )
            return
        output_text = f"Fitted {self.fitting_function_type} Parameters:\n"
        output_text += "--------------------------------------------------\n"
        for i in range(0, len(self.fitted_params), 3):
            amp, pos, fwhm = self.fitted_params [ i:i + 3 ]
            amp_err, pos_err, fwhm_err = (
                self.fitted_errors [ i ], self.fitted_errors [ i + 1 ], self.fitted_errors [ i + 2 ]
            ) if i + 2 < len(self.fitted_errors) else (0, 0, 0)
            output_text += (f"Peak {i // 3 + 1}:\n"
                            f"  Amplitude (Amp): {amp:.4g} ± {amp_err:.2g}\n"
                            f"  Position (Pos):  {pos:.4g} ± {pos_err:.2g}\n"
                            f"  FWHM:            {fwhm:.4g} ± {fwhm_err:.2g}\n")
            output_text += "--------------------------------------------------\n"
        self.params_text_edit.setText(output_text)

    def update_plot(self):
        self.ax.clear( )
        self.ax.plot(self.xga, self.yga, 'b-', label = 'Data')
        legend_entries = [ 'Data' ]
        if self.fitted_params is not None:
            self.ax.plot(self.xga, self._get_multi_function( )(self.xga, *self.fitted_params),
                         'g-', label = 'Fitted Curve', linewidth = 2)
            legend_entries.append('Fitted Curve')
            for i in range(0, len(self.fitted_params), 3):
                amp, pos, fwhm = self.fitted_params [ i:i + 3 ]
                amp_err, pos_err, fwhm_err = (
                    self.fitted_errors [ i ], self.fitted_errors [ i + 1 ], self.fitted_errors [ i + 2 ]
                ) if i + 2 < len(self.fitted_errors) else (0, 0, 0)
                label = (f'{self.fitting_function_type [ :5 ]} {i // 3 + 1}')
                self.ax.plot(self.xga, self._get_single_function( )(self.xga, amp, pos, fwhm),
                             '--', alpha = 0.7, linewidth = 1.5, label = label)
        else:
            for i, (amp, pos, fwhm) in enumerate(self.fixed_peaks):
                label = f'Initial Guess {i + 1}' if i == 0 else ""
                self.ax.plot(self.xga, self._get_single_function( )(self.xga, amp, pos, fwhm),
                             'r--', alpha = 0.5, label = label)
                if i == 0:
                    legend_entries.append('Initial Guesses')
        if self.amp is not None:
            if self.fwhm == 0: self.fwhm = 0.001
            self.temp_line, = self.ax.plot(self.xga,
                                           self._get_single_function( )(self.xga, self.amp, self.pos, self.fwhm),
                                           'g--', alpha = 0.5, label = 'Adjusting Width')
            if 'Adjusting Width' not in legend_entries:
                legend_entries.append('Adjusting Width')
        self.ax.set_ylabel(self.ylabel, fontsize = 18)
        self.ax.set_xlabel(self.xlabel, fontsize = 18)
        self.ax.grid(True)
        self.ax.legend(loc = 'best', fontsize = 18)
        self.ax.tick_params(axis = 'both', which = 'major', labelsize = 18)
        if self.fitted_params is None:
            pass
        self.canvas.draw( )

    def get_fitter_state(self):
        return {
            'type': 'GaussianFitterApp',
            'tab_title': self.objectName( ),
            'x_data': self.xga.tolist( ),
            'y_data': self.yga.tolist( ),
            'original_x_data': self._original_x_data.tolist( ),
            'original_y_data': self._original_y_data.tolist( ),
            'fitting_function_type': self.fitting_function_type,
            'xlabel': self.xlabel,
            'ylabel': self.ylabel,
            'slice_axis_name': self.slice_axis_name,
            'slice_value': self.slice_value,
            'slice_unit': self.slice_unit,
            'fixed_peaks': self.fixed_peaks,
            'fitted_params': self.fitted_params.tolist( ) if self.fitted_params is not None else None,
            'fitted_errors': self.fitted_errors.tolist( ) if self.fitted_errors is not None else None,
            'is_spline_corrected': self.is_spline_corrected
        }

    def set_fitter_state(self, state):
        self.setObjectName(state.get('tab_title', 'Gaussian Fitter'))
        self.setWindowTitle(state.get('tab_title', 'Gaussian Fitter'))

        self.xga = np.array(state.get('x_data'))
        self.yga = np.array(state.get('y_data'))
        self._original_x_data = np.array(
            state.get('original_x_data', self.xga))
        self._original_y_data = np.array(
            state.get('original_y_data', self.yga))
        self.fitting_function_type = state.get('fitting_function_type', 'Gaussian')
        self.xlabel = state.get('xlabel', 'X-axis')
        self.ylabel = state.get('ylabel', 'Y-axis')
        self.slice_axis_name = state.get('slice_axis_name', '')
        self.slice_value = state.get('slice_value')
        self.slice_unit = state.get('slice_unit', '')
        self.fixed_peaks = state.get('fixed_peaks', [ ])
        self.fitted_params = np.array(state.get('fitted_params')) if state.get('fitted_params') is not None else None
        self.fitted_errors = np.array(state.get('fitted_errors')) if state.get('fitted_errors') is not None else None
        self.is_spline_corrected = state.get('is_spline_corrected', False)

        self.export_button.setDisabled(self.fitted_params is None)

        self.display_fitted_parameters( )
        self.update_plot( )

    def _close_tab(self):
        parent_tab_widget = self.parent( )
        if parent_tab_widget and isinstance(parent_tab_widget, QTabWidget):
            tab_index = parent_tab_widget.indexOf(self)
            if tab_index != -1:
                parent_tab_widget.removeTab(tab_index)
        self.deleteLater( )


def signal_fitter_wrapper(parent, plot_data_item, is_x_slice, fitting_function_type, xlabel, ylabel, slice_axis_name,
                          slice_value, slice_unit, is_spline_corrected):
    x_data_full = plot_data_item.getData( ) [ 0 ]
    y_data_full = plot_data_item.getData( ) [ 1 ]
    view_range = plot_data_item.getViewBox( ).viewRange( )
    xlim_view = view_range [ 0 ]
    ylim_view = view_range [ 1 ]
    mask = np.logical_and(x_data_full >= xlim_view [ 0 ], x_data_full <= xlim_view [ 1 ])
    x_data_filtered = x_data_full [ mask ]
    y_data_filtered = y_data_full [ mask ]
    if len(x_data_filtered) < 3:
        QMessageBox.warning(None, "Fit Error",
                            "Not enough data points in visible range for fitting (need at least 3). Zoom in or adjust data.")
        return None
    fitter_app = GaussianFitterApp(parent, x_data_filtered, y_data_filtered, fitting_function_type, xlabel, ylabel,
                                   slice_axis_name, slice_value, slice_unit, is_spline_corrected)
    return fitter_app


def exponential_fitter_wrapper(parent, plot_data_item, xlabel, ylabel, slice_axis_name, slice_value, slice_unit,
                               is_spline_corrected):
    x_data_full = plot_data_item.getData( ) [ 0 ]
    y_data_full = plot_data_item.getData( ) [ 1 ]
    view_range = plot_data_item.getViewBox( ).viewRange( )
    xlim_view = view_range [ 0 ]
    mask = np.logical_and(x_data_full >= xlim_view [ 0 ], x_data_full <= xlim_view [ 1 ])
    x_data_filtered = x_data_full [ mask ]
    y_data_filtered = y_data_full [ mask ]
    if len(x_data_filtered) < 2:
        QMessageBox.warning(None, "Fit Error", "Not enough data points in visible range for fitting (need at least 2).")
        return None
    fitter_app = ExponentialFitterApp(parent, x_data_filtered, y_data_filtered, xlabel, ylabel,
                                      slice_axis_name, slice_value, slice_unit, is_spline_corrected)
    return fitter_app


class SignalPlotterApp(QMainWindow):
    def __init__(self):
        super( ).__init__( )
        self.base_title = "Kaalen-v2.0"
        self._current_project_file = None
        self._data_modified = False
        self._update_window_title( )

        self.setGeometry(100, 100, 1920, 1080)
        self.base_font_size = 12
        self.axis_label_font_size = 12
        self.axis_scale_font_size = 12
        self.central_widget = QWidget( )

        self.setCentralWidget(self.central_widget)
        self.main_layout = QGridLayout(self.central_widget)
        self.held_x_slices_count = 0
        self.held_y_slices_count = 0

        self.plot_colors = [
            (255, 0, 0),
            (0, 0, 255),
            (0, 200, 0),
            (255, 165, 0),
            (128, 0, 128),
            (0, 255, 255),
            (255, 0, 255),
            (150, 75, 0),
            (0, 128, 128),
            (255, 255, 0)
        ]
        self._initial_raw_x_values = None
        self._initial_raw_y_values = None
        self._initial_raw_signal_data = None
        self.current_x_values = None
        self.current_y_values = None
        self.current_signal_data = None
        self.is_spline_corrected = False

        self.x_values_interp = None
        self.y_values_interp = None
        self.signal_data_interp = None

        self.x_dim = 0
        self.y_dim = 0
        self.data_loaded = False

        self.x_legend_font_size = 14
        self.y_legend_font_size = 14

        self.axis_labels = {
            'signal_bottom': 'Probe wavenumber [cm\u207B\u00B9]',
            'signal_left': 'Time [ps]',
            'x_slice_bottom': 'Time [ps]',
            'x_slice_left': 'ΔOD',
            'y_slice_bottom': 'Probe wavenumber [cm\u207B\u00B9]',
            'y_slice_left': 'ΔOD'
        }

        self.current_slice_linewidth = 2
        self.x_slice_legend_unit = "cm^-1"
        self.y_slice_legend_unit = "ps"

        self._current_interp_method = "None"
        self._current_interp_multiplier = 1
        self._original_x_slice_data_main = None
        self._original_y_slice_data_main = None
        self._original_y_slice_data_main_y = None
        self._original_y_slice_data_main_x = None

        self.active_fitter_tabs = [ ]

        self.init_ui( )

        self.update_plots( )

        self.x_unit_input.textChanged.connect(self._set_data_modified)
        self.y_unit_input.textChanged.connect(self._set_data_modified)
        self.min_level_input.textChanged.connect(self._set_data_modified)
        self.max_level_input.textChanged.connect(self._set_data_modified)

    def _set_data_modified(self):
        if not self._data_modified:
            self._data_modified = True
            self._update_window_title( )

    def _update_window_title(self):
        title = self.base_title
        if self._current_project_file:
            project_name = os.path.basename(self._current_project_file)
            title += f" - {project_name}"
        else:
            title += " - Unsaved"
        if self._data_modified and self._current_project_file:
            title += " (Unsaved Changes)"
        elif self._data_modified and not self._current_project_file:
            title += " (Unsaved)"
        self.setWindowTitle(title)

    def _load_data_into_plots(self, x_raw, y_raw, z_raw):
        self._initial_raw_x_values = x_raw.copy( )
        self._initial_raw_y_values = y_raw.copy( )
        self._initial_raw_signal_data = z_raw.copy( )
        self.current_x_values = x_raw.copy( )
        self.current_y_values = y_raw.copy( )
        self.current_signal_data = z_raw.copy( )
        self.is_spline_corrected = False
        print(f"Original data shape: {self.current_signal_data.shape}")
        print(f"Original X range: {self.current_x_values.min( ):.2f} to {self.current_x_values.max( ):.2f}")
        print(f"Original Y range: {self.current_y_values.min( ):.2f} to {self.current_y_values.max( ):.2f}")
        self.data_loaded = True
        self._refresh_all_plots( )
        self._set_data_modified( )
        self.update_level_inputs( )

    def update_level_inputs(self):
        if self.data_loaded:
            self.min_level_input.setText(f"{self.signal_data_interp.min( ):.3f}")
            self.max_level_input.setText(f"{self.signal_data_interp.max( ):.3f}")
        else:
            self.min_level_input.setText("")
            self.max_level_input.setText("")

    def _refresh_all_plots(self, preserve_contour_levels=False, min_level=None, max_level=None,
                           preserve_plot_ranges=False, signal_xlim=None, signal_ylim=None,
                           x_slice_xlim=None, x_slice_ylim=None, y_slice_xlim=None, y_slice_ylim=None):
        if not self.data_loaded:
            self.update_plots( )
            return
        self.new_x_resolution = 1000
        self.new_y_resolution = 1000
        self.x_values_interp = np.linspace(self.current_x_values.min( ), self.current_x_values.max( ),
                                           self.new_x_resolution)
        self.y_values_interp = np.linspace(self.current_y_values.min( ), self.current_y_values.max( ),
                                           self.new_y_resolution)
        interp_func = RectBivariateSpline(self.current_y_values, self.current_x_values, self.current_signal_data)
        self.signal_data_interp = interp_func(self.y_values_interp, self.x_values_interp)
        self.x_dim = self.new_x_resolution
        self.y_dim = self.new_y_resolution
        self.x_slider.setRange(0, len(self.current_x_values) - 1)
        self.y_slider.setRange(0, len(self.current_y_values) - 1)
        self.image_item.setImage(self.signal_data_interp.T)
        self.image_item.setRect(pg.QtCore.QRectF(
            self.x_values_interp [ 0 ],
            self.y_values_interp [ 0 ],
            self.x_values_interp [ -1 ] - self.x_values_interp [ 0 ],
            self.y_values_interp [ -1 ] - self.y_values_interp [ 0 ]
        ))
        if preserve_contour_levels and min_level is not None and max_level is not None:
            self.min_level_input.setText(f"{min_level:.2f}")
            self.max_level_input.setText(f"{max_level:.2f}")
            self.image_item.setLevels((min_level, max_level))
        else:
            self.min_level_input.setText(f"{self.signal_data_interp.min( ):.3f}")
            self.max_level_input.setText(f"{self.signal_data_interp.max( ):.3f}")
            self.image_item.setLevels((self.signal_data_interp.min( ), self.signal_data_interp.max( )))
        if preserve_plot_ranges and signal_xlim is not None and signal_ylim is not None:
            self.signal_plot_widget.setXRange(signal_xlim [ 0 ], signal_xlim [ 1 ], padding = 0)
            self.signal_plot_widget.setYRange(signal_ylim [ 0 ], signal_ylim [ 1 ], padding = 0)
        else:
            self.signal_plot_widget.setXRange(self.x_values_interp.min( ), self.x_values_interp.max( ), padding = 0)
            self.signal_plot_widget.setYRange(self.y_values_interp.min( ), self.y_values_interp.max( ), padding = 0)
        if preserve_plot_ranges and x_slice_xlim is not None and x_slice_ylim is not None:
            self.x_slice_plot_widget.setXRange(x_slice_xlim [ 0 ], x_slice_xlim [ 1 ], padding = 0.05)
            self.x_slice_plot_widget.setYRange(x_slice_ylim [ 0 ], x_slice_ylim [ 1 ], padding = 0.05)
        else:
            self.x_slice_plot_widget.setXRange(self.current_y_values.min( ), self.current_y_values.max( ),
                                               padding = 0.05)
            self.x_slice_plot_widget.setYRange(self.current_signal_data.min( ), self.current_signal_data.max( ),
                                               padding = 0.05)
        if preserve_plot_ranges and y_slice_xlim is not None and y_slice_ylim is not None:
            self.y_slice_plot_widget.setXRange(y_slice_xlim [ 0 ], y_slice_xlim [ 1 ], padding = 0.05)
            self.y_slice_plot_widget.setYRange(y_slice_ylim [ 0 ], y_slice_ylim [ 1 ], padding = 0.05)
        else:
            self.y_slice_plot_widget.setXRange(self.current_x_values.min( ), self.current_x_values.max( ),
                                               padding = 0.05)
        if preserve_plot_ranges and y_slice_xlim is not None and y_slice_ylim is not None:
            self.y_slice_plot_widget.setXRange(y_slice_xlim [ 0 ], y_slice_xlim [ 1 ], padding = 0.05)
            self.y_slice_plot_widget.setYRange(y_slice_ylim [ 0 ], y_slice_ylim [ 1 ], padding = 0.05)
        else:
            self.y_slice_plot_widget.setXRange(self.current_x_values.min( ), self.current_x_values.max( ),
                                               padding = 0.05)
            self.y_slice_plot_widget.setYRange(self.current_signal_data.min( ), self.current_signal_data.max( ),
                                               padding = 0.05)
        x_idx_raw = self.x_slider.value( )
        y_idx_raw = self.y_slider.value( )
        x_pos_val_raw = self.current_x_values [ x_idx_raw ]
        y_pos_val_raw = self.current_y_values [ y_idx_raw ]
        self._original_y_slice_data_main_x = np.copy(self.current_x_values)
        self._original_y_slice_data_main_y = np.copy(self.current_signal_data [ y_idx_raw, : ])
        self._original_x_slice_data_main = np.copy(self.current_y_values)
        self._original_x_slice_data_main_y = np.copy(self.current_signal_data [ :, x_idx_raw ])
        self.update_plots( )
        self._update_spline_button_text( )

    def init_ui(self):
        self.menu_bar = self.menuBar( )
        self.file_menu = self.menu_bar.addMenu("&File")
        import_data_action = QAction("&Import Data...", self)
        import_data_action.setShortcut("Ctrl+I")
        import_data_action.setStatusTip("Import data from a file")
        import_data_action.triggered.connect(self.on_import_data_action_triggered)
        self.file_menu.addAction(import_data_action)
        save_project_action = QAction("&Save Project", self)
        save_project_action.setShortcut("Ctrl+S")
        save_project_action.setStatusTip("Save current project state")
        save_project_action.triggered.connect(self._save_project)
        self.file_menu.addAction(save_project_action)
        save_project_as_action = QAction("Save Project &As...", self)
        save_project_as_action.setShortcut("Ctrl+Shift+S")
        save_project_as_action.setStatusTip("Save current project with a new name")
        save_project_as_action.triggered.connect(self._save_project_as)
        self.file_menu.addAction(save_project_as_action)
        load_project_action = QAction("&Load Project...", self)
        load_project_action.setShortcut("Ctrl+L")
        load_project_action.setStatusTip("Load a previously saved project")
        load_project_action.triggered.connect(self._load_project)
        self.file_menu.addAction(load_project_action)
        exit_action = QAction("&Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        self.file_menu.addAction(exit_action)
        self.edit_menu = self.menu_bar.addMenu("&Edit")
        edit_names_action = QAction("Edit &Names...", self)
        edit_names_action.setShortcut("Ctrl+N")
        edit_names_action.setStatusTip("Edit axis labels for all plots")
        edit_names_action.triggered.connect(self._show_edit_names_dialog)
        self.edit_menu.addAction(edit_names_action)
        label_font = QFont("Times New Roman")
        label_font.setPointSize(self.base_font_size)
        axis_scale_font = QFont("Times New Roman")
        axis_scale_font.setPointSize(self.axis_scale_font_size)
        axis_label_font = QFont("Times New Roman")
        axis_label_font.setPointSize(self.axis_label_font_size)
        tick_length = 10
        tick_text_offset = 5
        self.tab_widget = QTabWidget( )
        self.main_plot_tab = QWidget( )
        self.tab_widget.setTabsClosable(True)
        self.tab_widget.addTab(self.main_plot_tab, "Main Plots")
        self.tab_widget.tabCloseRequested.connect(self.close_tab)
        self.tab_widget.tabBar( ).tabBarDoubleClicked.connect(self.rename_tab)
        self.main_plot_tab_layout = QGridLayout(self.main_plot_tab)
        self.signal_plot_widget = pg.PlotWidget( )
        self.signal_plot_widget.setLabel('bottom', self.axis_labels [ 'signal_bottom' ],
                                         **{ 'font-size': f'{self.axis_label_font_size}pt' })
        self.signal_plot_widget.setLabel('left', self.axis_labels [ 'signal_left' ],
                                         **{ 'font-size': f'{self.axis_label_font_size}pt' })
        self.signal_plot_widget.setBackground('w')
        self.signal_plot_widget.getViewBox( ).enableAutoRange(axis = pg.ViewBox.XYAxes, enable = False)
        self.signal_plot_widget.getViewBox( ).setMouseMode(pg.ViewBox.RectMode)
        self.signal_plot_widget.getViewBox( ).setMouseEnabled(x = True, y = True)
        bottom_axis = self.signal_plot_widget.getAxis('bottom')
        left_axis = self.signal_plot_widget.getAxis('left')
        bottom_axis.setTickFont(axis_scale_font)
        left_axis.setTickFont(axis_scale_font)
        bottom_axis.showMinorTicks = True
        left_axis.showMinorTicks = True
        bottom_axis.tickLength = tick_length
        left_axis.tickLength = tick_length
        bottom_axis.setStyle(tickTextOffset = tick_text_offset)
        left_axis.setStyle(tickTextOffset = tick_text_offset)
        self.image_item = pg.ImageItem( )
        self.signal_plot_widget.addItem(self.image_item)
        colors_seismic_base = [
            (0, 0, 128, 255),
            (0, 0, 255, 255),
            (0, 255, 255, 255),
            (255, 255, 255, 255),
            (255, 255, 0, 255),
            (255, 0, 0, 255),
            (128, 0, 0, 255)
        ]
        positions_seismic_base = np.linspace(0.0, 0.99, len(colors_seismic_base))
        final_colors = colors_seismic_base + [ (255, 255, 255, 255) ]
        final_positions = np.append(positions_seismic_base, 1.0)
        self.seismic_colormap = pg.ColorMap(final_positions, final_colors)
        self.image_item.setLookupTable(self.seismic_colormap.getLookupTable( ))
        self.cursor_x_line = pg.InfiniteLine(angle = 90, movable = False,
                                             pen = pg.mkPen('k', width = 1, style = Qt.DotLine))
        self.cursor_y_line = pg.InfiniteLine(angle = 0, movable = False,
                                             pen = pg.mkPen('k', width = 1, style = Qt.DotLine))
        self.signal_plot_widget.addItem(self.cursor_x_line)
        self.signal_plot_widget.addItem(self.cursor_y_line)
        self.cursor_x_line.setVisible(False)
        self.cursor_y_line.setVisible(False)
        self.signal_plot_widget.scene( ).sigMouseMoved.connect(self.update_signal_cursor_pos)
        vb_menu = self.signal_plot_widget.getPlotItem( ).getViewBox( ).menu
        self.x_slice_plot_widget = pg.PlotWidget( )
        self.x_slice_plot_widget.setLabel('bottom', self.axis_labels [ 'x_slice_bottom' ],
                                          **{ 'font-size': f'{self.axis_label_font_size}pt' })
        self.x_slice_plot_widget.setLabel('left', self.axis_labels [ 'x_slice_left' ],
                                          **{ 'font-size': f'{self.axis_label_font_size}pt' })
        self.x_slice_plot_widget.setBackground('w')
        self.x_slice_plot_widget.getViewBox( ).enableAutoRange(axis = pg.ViewBox.XYAxes, enable = False)
        self.x_slice_plot_widget.getViewBox( ).setMouseMode(pg.ViewBox.RectMode)
        self.x_slice_plot_widget.getViewBox( ).setMouseEnabled(x = True, y = True)
        x_slice_bottom_axis = self.x_slice_plot_widget.getAxis('bottom')
        x_slice_left_axis = self.x_slice_plot_widget.getAxis('left')
        x_slice_bottom_axis.setTickFont(axis_scale_font)
        x_slice_left_axis.setTickFont(axis_scale_font)
        x_slice_bottom_axis.showMinorTicks = True
        x_slice_left_axis.showMinorTicks = True
        x_slice_bottom_axis.tickLength = tick_length
        x_slice_left_axis.tickLength = tick_text_offset
        x_slice_bottom_axis.setStyle(tickTextOffset = tick_text_offset)
        x_slice_left_axis.setStyle(tickTextOffset = tick_text_offset)
        self.x_slice_curve = self.x_slice_plot_widget.plot(pen = pg.mkPen('b', width = self.current_slice_linewidth))
        self.x_slice_legend = self.x_slice_plot_widget.addLegend( )
        self._apply_legend_font_size(self.x_slice_legend, self.x_legend_font_size)
        x_vb_menu = self.x_slice_plot_widget.getPlotItem( ).getViewBox( ).menu
        self.change_x_linewidth_action = x_vb_menu.addAction("Change Line Thickness")
        self.change_x_linewidth_action.triggered.connect(lambda: self._show_linewidth_dialog(self.x_slice_plot_widget))
        self.x_slice_plot_widget.scene( ).sigMouseMoved.connect(self.update_x_slice_cursor_pos)
        self.x_hold_button = QPushButton("Hold")
        self.x_hold_button.setFont(label_font)
        self.x_hold_button.clicked.connect(self.hold_x_slice_plot)
        self.x_hold_button.setMinimumHeight(30)
        self.x_unit_input = QLineEdit(self)
        self.x_unit_input.setFont(label_font)
        self.x_unit_input.setPlaceholderText("Enter X-slice unit (e.g., ps)")
        self.x_unit_input.setText(self.x_slice_legend_unit)
        self.x_unit_input.setFixedWidth(150)
        self.x_clear_button = QPushButton("Clear")
        self.x_clear_button.setFont(label_font)
        self.x_clear_button.clicked.connect(self.clear_x_slice_plots)
        self.x_clear_button.setMinimumHeight(30)
        self.x_fit_button = QPushButton("Fit Exponential")
        self.x_fit_button.setFont(label_font)
        self.x_fit_button.setMinimumHeight(30)
        self.x_fit_button.clicked.connect(self._open_x_fitter_tab)
        self.y_slice_plot_widget = pg.PlotWidget( )
        self.y_slice_plot_widget.setLabel('bottom', self.axis_labels [ 'y_slice_bottom' ],
                                          **{ 'font-size': f'{self.axis_label_font_size}pt' })
        self.y_slice_plot_widget.setLabel('left', self.axis_labels [ 'y_slice_left' ],
                                          **{ 'font-size': f'{self.axis_label_font_size}pt' })
        self.y_slice_plot_widget.setBackground('w')
        self.y_slice_plot_widget.getViewBox( ).enableAutoRange(axis = pg.ViewBox.XYAxes, enable = False)
        self.y_slice_plot_widget.getViewBox( ).setMouseMode(pg.ViewBox.RectMode)
        self.y_slice_plot_widget.getViewBox( ).setMouseEnabled(x = True, y = True)
        y_slice_bottom_axis = self.y_slice_plot_widget.getAxis('bottom')
        y_slice_left_axis = self.y_slice_plot_widget.getAxis('left')
        y_slice_bottom_axis.setTickFont(axis_scale_font)
        y_slice_left_axis.setTickFont(axis_scale_font)
        y_slice_bottom_axis.showMinorTicks = True
        y_slice_left_axis.showMinorTicks = True
        y_slice_bottom_axis.tickLength = tick_length
        y_slice_left_axis.tickLength = tick_length
        y_slice_bottom_axis.setStyle(tickTextOffset = tick_text_offset)
        y_slice_left_axis.setStyle(tickTextOffset = tick_text_offset)
        self.y_slice_curve = self.y_slice_plot_widget.plot(pen = pg.mkPen('r', width = self.current_slice_linewidth))
        self.y_slice_legend = self.y_slice_plot_widget.addLegend( )
        self._apply_legend_font_size(self.y_slice_legend, self.y_legend_font_size)
        y_vb_menu = self.y_slice_plot_widget.getPlotItem( ).getViewBox( ).menu
        self.change_y_linewidth_action = y_vb_menu.addAction("Change Line Thickness...")
        self.change_y_linewidth_action.triggered.connect(lambda: self._show_linewidth_dialog(self.y_slice_plot_widget))
        self.y_slice_plot_widget.scene( ).sigMouseMoved.connect(self.update_y_slice_cursor_pos)
        self.y_hold_button = QPushButton("Hold")
        self.y_hold_button.setFont(label_font)
        self.y_hold_button.clicked.connect(self.hold_y_slice_plot)
        self.y_hold_button.setMinimumHeight(30)
        self.y_unit_input = QLineEdit(self)
        self.y_unit_input.setFont(label_font)
        self.y_unit_input.setPlaceholderText("Enter Y-slice unit (e.g., cm^-1)")
        self.y_unit_input.setText(self.y_slice_legend_unit)
        self.y_unit_input.setFixedWidth(150)
        self.y_clear_button = QPushButton("Clear")
        self.y_clear_button.setFont(label_font)
        self.y_clear_button.clicked.connect(self.clear_y_slice_plots)
        self.y_clear_button.setMinimumHeight(30)
        self.y_fit_function_selector = QComboBox(self)
        self.y_fit_function_selector.addItem("Gaussian")
        self.y_fit_function_selector.addItem("Lorentzian")
        self.y_fit_function_selector.setMinimumHeight(30)
        self.y_fit_function_selector.setFont(label_font)
        self.y_fit_button = QPushButton("Fit")
        self.y_fit_button.setFont(label_font)
        self.y_fit_button.setMinimumHeight(30)
        self.y_fit_button.clicked.connect(
            lambda: self._open_fitter_tab(
                self.y_slice_curve,
                False,
                self.y_fit_function_selector.currentText( ),
                self.axis_labels [ 'y_slice_bottom' ],
                self.axis_labels [ 'y_slice_left' ],
                self._strip_html_tags(self.axis_labels [ 'signal_left' ]),
                self.current_y_values [ self.y_slider.value( ) ],
                self.y_unit_input.text( ),
                self.is_spline_corrected
            )
        )
        self.x_slider = QSlider(Qt.Horizontal)
        self.x_slider.setTickPosition(QSlider.TicksBelow)
        self.x_slider.valueChanged.connect(self.update_plots)
        self.y_slider = QSlider(Qt.Horizontal)
        self.y_slider.setTickPosition(QSlider.TicksBelow)
        self.y_slider.valueChanged.connect(self.update_plots)
        self.min_level_input = QLineEdit(self)
        self.min_level_input.setPlaceholderText("-0.5")
        self.min_level_input.setFont(label_font)
        self.min_level_input.setValidator(QDoubleValidator( ))
        self.min_level_input.editingFinished.connect(self.update_plots)
        self.min_level_input.setFixedWidth(100)
        self.max_level_input = QLineEdit(self)
        self.max_level_input.setFont(label_font)
        self.max_level_input.setValidator(QDoubleValidator( ))
        self.max_level_input.editingFinished.connect(self.update_plots)
        self.max_level_input.setFixedWidth(100)
        self.x_label = QLabel(f"Probe:")
        self.x_label.setFont(label_font)
        self.x_input = QLineEdit(self)
        self.x_input.setValidator(QDoubleValidator( ))
        self.x_input.setFixedWidth(100)
        self.x_input.editingFinished.connect(self.update_x_slider_from_input)
        self.y_label = QLabel(f"Delay time:")
        self.y_label.setFont(label_font)
        self.y_input = QLineEdit(self)
        self.y_input.setValidator(QDoubleValidator( ))
        self.y_input.setFixedWidth(100)
        self.y_input.editingFinished.connect(self.update_y_slider_from_input)
        self.min_level_label = QLabel("Contour min:")
        self.min_level_label.setFont(label_font)
        self.max_level_label = QLabel("Contour max:")
        self.max_level_label.setFont(label_font)
        self.cursor_pos_label = QLabel("Cursor: (X: -, Y: -)")
        self.cursor_pos_label.setFont(label_font)
        self.spline_baseline_button = QPushButton("Apply Spline Baseline")
        self.spline_baseline_button.setFont(label_font)
        self.spline_baseline_button.setMinimumHeight(30)
        self.spline_baseline_button.clicked.connect(self._toggle_spline_correction)
        self.interp_method_label = QLabel("Interpolate probe axis:")
        self.interp_method_combo = QComboBox(self)
        self.interp_method_combo.addItems([ "None", "linear", "cubic", "nearest" ])
        self.interp_method_combo.setFont(label_font)
        self.interp_method_combo.setMinimumHeight(30)
        self.interp_method_combo.currentIndexChanged.connect(self._apply_interpolation_to_all_plots)
        self.interp_multiplier_label = QLabel("")
        self.interp_multiplier_combo = QComboBox(self)
        self.interp_multiplier_combo.addItems([ "x1", "x2", "x3", "x5" ])
        self.interp_multiplier_combo.setFont(label_font)
        self.interp_multiplier_combo.setMinimumHeight(30)
        self.interp_multiplier_combo.currentIndexChanged.connect(self._apply_interpolation_to_all_plots)

        self.global_fit_button = QPushButton("Open Global Fit Tab")
        self.global_fit_button.setFont(label_font)
        self.global_fit_button.setMinimumHeight(30)
        self.global_fit_button.clicked.connect(self._launch_global_fit_tab)

        self.pfid_fit_button = QPushButton("Open PFID Fit Tab")
        self.pfid_fit_button.setFont(label_font)
        self.pfid_fit_button.setMinimumHeight(30)
        self.pfid_fit_button.clicked.connect(self._launch_pfid_fit_tab)

        self.main_plot_tab_layout.addWidget(self.signal_plot_widget, 0, 0, 2, 2)
        x_slice_controls_layout = QHBoxLayout( )
        x_slice_controls_layout.addWidget(self.x_hold_button)
        x_slice_controls_layout.addWidget(self.x_unit_input)
        x_slice_controls_layout.addWidget(self.x_clear_button)
        x_slice_controls_layout.addWidget(self.x_fit_button)
        x_slice_controls_layout.setSpacing(10)
        x_slice_plot_and_buttons_layout = QVBoxLayout( )
        x_slice_plot_and_buttons_layout.addWidget(self.x_slice_plot_widget)
        x_slice_plot_and_buttons_layout.addLayout(x_slice_controls_layout)
        self.main_plot_tab_layout.addLayout(x_slice_plot_and_buttons_layout, 0, 2, 1, 1)
        y_slice_controls_layout = QHBoxLayout( )
        y_slice_controls_layout.addWidget(self.y_hold_button)
        y_slice_controls_layout.addWidget(self.y_unit_input)
        y_slice_controls_layout.addWidget(self.y_clear_button)
        y_slice_controls_layout.addWidget(QLabel("Fit Type:"))
        y_slice_controls_layout.addWidget(self.y_fit_function_selector)
        y_slice_controls_layout.addWidget(self.y_fit_button)
        y_slice_controls_layout.setSpacing(10)
        y_slice_plot_and_buttons_layout = QVBoxLayout( )
        y_slice_plot_and_buttons_layout.addWidget(self.y_slice_plot_widget)
        y_slice_plot_and_buttons_layout.addLayout(y_slice_controls_layout)
        self.main_plot_tab_layout.addLayout(y_slice_plot_and_buttons_layout, 1, 2, 1, 1)
        slider_controls_layout = QVBoxLayout( )

        x_slider_layout = QHBoxLayout( )
        x_slider_layout.addWidget(self.x_label)
        x_slider_layout.addWidget(self.x_input)
        x_slider_layout.addWidget(self.x_slider)
        slider_controls_layout.addLayout(x_slider_layout)

        y_slider_layout = QHBoxLayout( )
        y_slider_layout.addWidget(self.y_label)
        y_slider_layout.addWidget(self.y_input)
        y_slider_layout.addWidget(self.y_slider)
        slider_controls_layout.addLayout(y_slider_layout)

        level_inputs_layout = QHBoxLayout( )
        level_inputs_layout.addWidget(self.min_level_label)
        level_inputs_layout.addWidget(self.min_level_input)
        level_inputs_layout.addSpacing(20)
        level_inputs_layout.addWidget(self.max_level_label)
        level_inputs_layout.addWidget(self.max_level_input)
        level_inputs_layout.addStretch(1)
        level_inputs_layout.addWidget(self.global_fit_button)
        level_inputs_layout.addWidget(self.pfid_fit_button)
        level_inputs_layout.addSpacing(5)
        level_inputs_layout.addWidget(self.spline_baseline_button)
        level_inputs_layout.addWidget(self.interp_method_label)
        level_inputs_layout.addWidget(self.interp_method_combo)
        level_inputs_layout.addWidget(self.interp_multiplier_label)
        level_inputs_layout.addWidget(self.interp_multiplier_combo)
        slider_controls_layout.addLayout(level_inputs_layout)
        slider_controls_layout.addWidget(self.cursor_pos_label)
        self.main_plot_tab_layout.addLayout(slider_controls_layout, 2, 0, 1, 3)
        self.main_plot_tab_layout.setColumnStretch(0, 2)
        self.main_plot_tab_layout.setColumnStretch(1, 2)
        self.main_plot_tab_layout.setColumnStretch(2, 3)
        self.main_plot_tab_layout.setRowStretch(0, 3)
        self.main_plot_tab_layout.setRowStretch(1, 3)
        self.main_plot_tab_layout.setRowStretch(2, 1)
        self.main_plot_tab_layout.setRowStretch(3, 0)
        self.main_layout.addWidget(self.tab_widget, 0, 0, 1, 1)
        self._update_spline_button_text( )

    def update_x_slider_from_input(self):
        if not self.data_loaded:
            return
        try:
            val = float(self.x_input.text( ))
            closest_index = np.argmin(np.abs(self.current_x_values - val))
            self.x_slider.setValue(closest_index)
        except ValueError:
            pass

    def update_y_slider_from_input(self):
        if not self.data_loaded:
            return
        try:
            val = float(self.y_input.text( ))
            closest_index = np.argmin(np.abs(self.current_y_values - val))
            self.y_slider.setValue(closest_index)
        except ValueError:
            pass

    def _launch_global_fit_tab(self):
        if not self.data_loaded or self.current_signal_data is None:
            QMessageBox.warning(self, "No Data", "Please import a dataset before running a global fit.")
            return
        x_axis_label, x_axis_unit = _parse_label_and_unit(self.axis_labels [ 'signal_bottom' ])
        y_axis_label, y_axis_unit = _parse_label_and_unit(self.axis_labels [ 'signal_left' ])

        global_fit_tab = GlobalFitApp(
            x_axis_data = self.current_x_values,
            y_axis_data = self.current_y_values,
            two_d_spectrum_data = self.current_signal_data,
            parent = self.tab_widget,
            x_axis_label = x_axis_label,
            y_axis_label = y_axis_label,
            x_axis_unit = x_axis_unit,
            y_axis_unit = y_axis_unit,
            font_size = self.axis_label_font_size
        )

        tab_index = self.tab_widget.addTab(global_fit_tab, "Global Fit")
        self.tab_widget.setCurrentIndex(tab_index)
        self.active_fitter_tabs.append({ 'widget': global_fit_tab, 'tab_index': tab_index })
        global_fit_tab.destroyed.connect(partial(self._remove_fitter_tab, global_fit_tab))
        self._set_data_modified( )

    def _launch_pfid_fit_tab(self):
        if not self.data_loaded or self.current_signal_data is None:
            QMessageBox.warning(self, "No Data", "Please import a dataset before running the PFID fit.")
            return

        x_axis_label, x_axis_unit = _parse_label_and_unit(self.axis_labels [ 'signal_bottom' ])
        y_axis_label, y_axis_unit = _parse_label_and_unit(self.axis_labels [ 'signal_left' ])

        pfid_fit_tab = PFIDFitterApp(
            main_window = self,
            x_axis_label = x_axis_label,
            y_axis_label = y_axis_label,
            x_axis_unit = x_axis_unit,
            y_axis_unit = y_axis_unit,
            font_size = self.axis_label_font_size
        )

        tab_index = self.tab_widget.addTab(pfid_fit_tab, "PFID Fit")
        self.tab_widget.setCurrentIndex(tab_index)
        self.active_fitter_tabs.append({ 'widget': pfid_fit_tab, 'tab_index': tab_index })
        pfid_fit_tab.destroyed.connect(partial(self._remove_fitter_tab, pfid_fit_tab))
        self._set_data_modified( )

    def _get_interpolated_1d_data(self, original_x, original_y, method, multiplier):
        if method == "None" or len(original_x) < 2:
            return np.copy(original_x), np.copy(original_y)
        try:
            target_n_points = int(len(original_x) * multiplier)
            if target_n_points < 2:
                target_n_points = 2
            x_interp = np.linspace(original_x.min( ), original_x.max( ), target_n_points)
            f_interp = interp1d(original_x, original_y, kind = method,
                                fill_value = "extrapolate")
            y_interp = f_interp(x_interp)
            return x_interp, y_interp
        except Exception as e:
            print(f"Error during 1D interpolation ({method}, x{multiplier}): {e}")
            return np.copy(original_x), np.copy(original_y)

    def _apply_interpolation_to_all_plots(self):
        self._current_interp_method = self.interp_method_combo.currentText( )
        self._current_interp_multiplier = int(self.interp_multiplier_combo.currentText( ).replace('x', ''))
        self.update_plots( )
        print(
            f"Applying interpolation to fitter tabs: Method={self._current_interp_method}, Multiplier={self._current_interp_multiplier}")
        for tab_info in self.active_fitter_tabs:
            fitter_app_instance = tab_info [ 'widget' ]
            if isinstance(fitter_app_instance, GaussianFitterApp):
                fitter_app_instance.apply_interpolation_settings(self._current_interp_method,
                                                                 self._current_interp_multiplier)
        self._set_data_modified( )

    def rename_tab(self, index):
        if index == 0:
            return
        old_name = self.tab_widget.tabText(index)
        new_name, ok = QInputDialog.getText(self, "Rename Tab", "Enter new tab name:", QLineEdit.Normal, old_name)
        if ok and new_name:
            self.tab_widget.setTabText(index, new_name)
            self.tab_widget.widget(index).setObjectName(new_name)
            self._set_data_modified( )

    def close_tab(self, index):
        if index == 0:
            QMessageBox.information(self, "Cannot Close Tab", "The main plot window cannot be closed.")
            return

        widget = self.tab_widget.widget(index)
        if widget is not None:
            is_fitter_tab = False
            for i, tab_info in enumerate(self.active_fitter_tabs):
                if tab_info [ 'widget' ] == widget:
                    self.active_fitter_tabs.pop(i)
                    is_fitter_tab = True
                    break
            self.tab_widget.removeTab(index)
            widget.deleteLater( )
            self._set_data_modified( )

    def _apply_legend_font_size(self, legend_item, new_size):
        if legend_item:
            for sample, label_item in legend_item.items:
                font = label_item.font( )
                font.setPointSize(new_size)
                label_item.setFont(font)

    def _show_linewidth_dialog(self, plot_widget):
        dialog = LineThicknessDialog(self, initial_thickness = self.current_slice_linewidth)
        if dialog.exec_( ) == QDialog.Accepted:
            new_thickness = dialog.get_thickness( )
            self.current_slice_linewidth = new_thickness
            if plot_widget == self.x_slice_plot_widget:
                current_pen = self.x_slice_curve.opts [ 'pen' ]
                self.x_slice_curve.setPen(
                    pg.mkPen(current_pen.color( ), width = new_thickness, style = Qt.SolidLine))
            elif plot_widget == self.y_slice_plot_widget:
                current_pen = self.y_slice_curve.opts [ 'pen' ]
                self.y_slice_curve.setPen(
                    pg.mkPen(current_pen.color( ), width = new_thickness, style = old_pen.style( )))
            for item in plot_widget.listDataItems( ):
                if (plot_widget == self.x_slice_plot_widget and item == self.x_slice_curve) or \
                        (plot_widget == self.y_slice_plot_widget and item == self.y_slice_curve):
                    continue
                if item.opts.get('pen') is not None:
                    old_pen = item.opts [ 'pen' ]
                    new_pen = pg.mkPen(old_pen.color( ), width = new_thickness, style = old_pen.style( ))
                    item.setPen(new_pen)
            self._set_data_modified( )

    def update_plots(self):
        if not self.data_loaded:
            self.image_item.setImage(np.zeros((2, 2)).T)
            self.image_item.setRect(pg.QtCore.QRectF(0, 0, 1, 1))
            self.signal_plot_widget.setXRange(0, 1, padding = 0)
            self.signal_plot_widget.setYRange(0, 1, padding = 0)
            self.x_slice_curve.setData([ ], [ ])
            self.y_slice_curve.setData([ ], [ ])
            self.x_slice_plot_widget.setXRange(0, 1, padding = 0.05)
            self.x_slice_plot_widget.setYRange(0, 1, padding = 0.05)
            self.y_slice_plot_widget.setXRange(0, 1, padding = 0.05)
            self.y_slice_plot_widget.setYRange(0, 1, padding = 0.05)
            self.x_label.setText(f"Probe:")
            self.x_input.setText("-")
            self.y_label.setText(f"Delay time:")
            self.y_input.setText("-")
            self.cursor_x_line.setVisible(False)
            self.cursor_y_line.setVisible(False)
            self.min_level_input.setText("")
            self.max_level_input.setText("")
            self.x_slider.setEnabled(False)
            self.y_slider.setEnabled(False)
            self.x_input.setEnabled(False)
            self.y_input.setEnabled(False)
            self.min_level_input.setEnabled(False)
            self.max_level_input.setEnabled(False)
            self.x_hold_button.setEnabled(False)
            self.x_unit_input.setEnabled(False)
            self.x_clear_button.setEnabled(False)
            self.x_fit_button.setEnabled(False)
            self.y_hold_button.setEnabled(False)
            self.y_unit_input.setEnabled(False)
            self.y_clear_button.setEnabled(False)
            self.y_fit_button.setEnabled(False)
            self.y_fit_function_selector.setEnabled(False)
            self.spline_baseline_button.setEnabled(False)
            self.interp_method_combo.setEnabled(False)
            self.interp_multiplier_combo.setEnabled(False)
            self.global_fit_button.setEnabled(False)
            self.pfid_fit_button.setEnabled(False)
            return
        self.x_slider.setEnabled(True)
        self.y_slider.setEnabled(True)
        self.x_input.setEnabled(True)
        self.y_input.setEnabled(True)
        self.min_level_input.setEnabled(True)
        self.max_level_input.setEnabled(True)
        self.x_hold_button.setEnabled(True)
        self.x_unit_input.setEnabled(True)
        self.x_clear_button.setEnabled(True)
        self.x_fit_button.setEnabled(True)
        self.y_hold_button.setEnabled(True)
        self.y_unit_input.setEnabled(True)
        self.y_clear_button.setEnabled(True)
        self.y_fit_button.setEnabled(True)
        self.y_fit_function_selector.setEnabled(True)
        self.spline_baseline_button.setEnabled(True)
        self.interp_method_combo.setEnabled(True)
        self.interp_multiplier_combo.setEnabled(True)
        self.global_fit_button.setEnabled(True)
        self.pfid_fit_button.setEnabled(True)
        x_idx_raw = self.x_slider.value( )
        y_idx_raw = self.y_slider.value( )
        x_pos_val_raw = self.current_x_values [ x_idx_raw ]
        y_pos_val_raw = self.current_y_values [ y_idx_raw ]

        x_val_display = x_pos_val_raw
        y_val_display = y_pos_val_raw

        self.x_label.setText(f"Probe:")
        self.x_input.setText(f"{x_pos_val_raw:.1f}")
        self.y_label.setText(f"Delay time:")
        self.y_input.setText(f"{y_pos_val_raw:.1f}")
        original_y_slice_x_data = self.current_x_values
        original_y_slice_y_data = self.current_signal_data [ y_idx_raw, : ]
        original_x_slice_x_data = self.current_y_values
        original_x_slice_y_data = self.current_signal_data [ :, x_idx_raw ]
        interp_x_slice_x, interp_x_slice_y = original_x_slice_x_data, original_x_slice_y_data
        interp_y_slice_x, interp_y_slice_y = self._get_interpolated_1d_data(
            original_y_slice_x_data, original_y_slice_y_data,
            self._current_interp_method, self._current_interp_multiplier
        )
        self.y_slice_curve.setData(interp_y_slice_x, interp_y_slice_y)
        self.x_slice_curve.setData(interp_x_slice_x, interp_x_slice_y)
        vb_y_slice = self.y_slice_plot_widget.getViewBox( )
        x_auto_y_slice, y_auto_y_slice = vb_y_slice.autoRangeEnabled( )
        if x_auto_y_slice:
            self.y_slice_plot_widget.setXRange(interp_y_slice_x.min( ), interp_y_slice_x.max( ),
                                               padding = 0.05)
        if y_auto_y_slice:
            self.y_slice_plot_widget.setYRange(interp_y_slice_y.min( ), interp_y_slice_y.max( ),
                                               padding = 0.05)
        vb_x_slice = self.x_slice_plot_widget.getViewBox( )
        x_auto_x_slice, y_auto_x_slice = vb_x_slice.autoRangeEnabled( )
        if x_auto_x_slice:
            self.x_slice_plot_widget.setXRange(interp_x_slice_x.min( ), interp_x_slice_x.max( ),
                                               padding = 0.05)
        if y_auto_x_slice:
            self.x_slice_plot_widget.setYRange(interp_x_slice_y.min( ), interp_x_slice_y.max( ),
                                               padding = 0.05)
        try:
            min_level = float(self.min_level_input.text( ))
        except ValueError:
            min_level = self.signal_data_interp.min( )
        try:
            max_level = float(self.max_level_input.text( ))
        except ValueError:
            max_level = self.signal_data_interp.max( )
        if min_level > max_level:
            min_level, max_level = max_level, min_level
        self.image_item.setLevels((min_level, max_level))
        self.cursor_x_line.setPos(x_pos_val_raw)
        self.cursor_y_line.setPos(y_pos_val_raw)
        self.cursor_x_line.setVisible(True)
        self.cursor_y_line.setVisible(True)
        self.cursor_pos_label.setText(f"Cursor: (X: {x_val_display:.2f}, Y: {y_val_display:.2f})")

    def update_signal_cursor_pos(self, evt):
        if not self.data_loaded:
            self.cursor_pos_label.setText("Cursor: (X: -, Y: -)")
            return
        pos = evt
        if self.signal_plot_widget.sceneBoundingRect( ).contains(pos):
            mousePoint = self.signal_plot_widget.plotItem.vb.mapSceneToView(pos)
            x_val_display = mousePoint.x( )
            y_val_display = mousePoint.y( )
            self.cursor_x_line.setVisible(True)
            self.cursor_y_line.setVisible(True)
            self.cursor_x_line.setPos(x_val_display)
            self.cursor_y_line.setPos(y_val_display)
            self.cursor_pos_label.setText(f"Cursor: (X: {x_val_display:.2f}, Y: {y_val_display:.2f})")
        else:
            self.cursor_x_line.setVisible(False)
            self.cursor_y_line.setVisible(False)
            self.cursor_pos_label.setText("Cursor: (X: -, Y: -)")

    def update_x_slice_cursor_pos(self, evt):
        if not self.data_loaded: return
        pos = evt
        if self.x_slice_plot_widget.sceneBoundingRect( ).contains(pos):
            mousePoint = self.x_slice_plot_widget.plotItem.vb.mapSceneToView(pos)
            y_val_current = mousePoint.x( )
            amp_val = mousePoint.y( )
            self.cursor_pos_label.setText(f"X-Slice Cursor: ({y_val_current:.2f},{amp_val:.2g})")
        else:
            self.cursor_pos_label.setText("Cursor: (X: -, Y: -)")

    def update_y_slice_cursor_pos(self, evt):
        if not self.data_loaded: return
        pos = evt
        if self.y_slice_plot_widget.sceneBoundingRect( ).contains(pos):
            mousePoint = self.y_slice_plot_widget.plotItem.vb.mapSceneToView(pos)
            x_val_current = mousePoint.x( )
            amp_val = mousePoint.y( )
            self.cursor_pos_label.setText(f"Y-Slice Cursor: ({x_val_current:.2f}, {amp_val:.2g})")
        else:
            self.cursor_pos_label.setText("Cursor: (X: -, Y: -)")

    def _format_unit_for_display(self, unit_string):
        unit_string = unit_string.replace("^-1", "\u207B\u00B9")
        unit_string = unit_string.replace("^2", "\u00B2")
        unit_string = unit_string.replace("^3", "\u00B3")
        unit_string = unit_string.replace("^-3", "\u207B\u00B3")
        unit_string = unit_string.replace("^-4", "\u207B\u2074")
        unit_string = unit_string.replace("_1", "\u2081")
        unit_string = unit_string.replace("_2", "\u2082")
        unit_string = unit_string.replace("_3", "\u2083")
        return unit_string

    def _strip_html_tags(self, text):
        text = re.sub(r'<sup[^>]*>.*?</sup>', '', text)
        text = re.sub(r'<sub[^>]*>.*?</sub>', '', text)
        text = text.replace('<sup>', '').replace('</sup>', '')
        text = text.replace('<sub>', '').replace('</sub>', '')
        text = re.sub(r'\[.*?\]', '', text).strip( )
        return text

    def hold_x_slice_plot(self):
        if not self.data_loaded: return
        x_idx_raw = self.x_slider.value( )
        x_pos_val = self.current_x_values [ x_idx_raw ]
        self.held_x_slices_count += 1
        x_data = self.x_slice_curve.getData( ) [ 0 ]
        y_data = self.x_slice_curve.getData( ) [ 1 ]
        color_index = (self.held_x_slices_count - 1) % len(self.plot_colors)
        color = self.plot_colors [ color_index ]
        pen = pg.mkPen(color, width = self.current_slice_linewidth, style = Qt.SolidLine)
        unit_text = self.x_unit_input.text( ).strip( )
        formatted_unit_text = self._format_unit_for_display(unit_text)
        if formatted_unit_text:
            name = f'{x_pos_val:.0f} {formatted_unit_text}'
        else:
            name = f'{x_pos_val:.0f}'
        self.x_slice_plot_widget.plot(x_data, y_data, pen = pen, name = name)
        self._set_data_modified( )

    def clear_x_slice_plots(self):
        items_to_remove = [ item for item in self.x_slice_plot_widget.listDataItems( ) if item != self.x_slice_curve ]
        if items_to_remove:
            for item in items_to_remove:
                self.x_slice_plot_widget.removeItem(item)
            self.held_x_slices_count = 0
            self._set_data_modified( )

    def hold_y_slice_plot(self):
        if not self.data_loaded: return
        y_idx_raw = self.y_slider.value( )
        y_pos_val = self.current_y_values [ y_idx_raw ]
        self.held_y_slices_count += 1
        x_data = self.y_slice_curve.getData( ) [ 0 ]
        y_data = self.y_slice_curve.getData( ) [ 1 ]
        color_index = (self.held_y_slices_count - 1) % len(self.plot_colors)
        color = self.plot_colors [ color_index ]
        pen = pg.mkPen(color, width = self.current_slice_linewidth, style = Qt.SolidLine)
        unit_text = self.y_unit_input.text( ).strip( )
        formatted_unit_text = self._format_unit_for_display(unit_text)
        if formatted_unit_text:
            name = f'{y_pos_val:.0f} {formatted_unit_text}'
        else:
            name = f'{y_pos_val:.0f}'
        self.y_slice_plot_widget.plot(x_data, y_data, pen = pen, name = name)
        self._set_data_modified( )

    def clear_y_slice_plots(self):
        items_to_remove = [ item for item in self.y_slice_plot_widget.listDataItems( ) if item != self.y_slice_curve ]
        if items_to_remove:
            for item in items_to_remove:
                self.y_slice_plot_widget.removeItem(item)
            self.held_y_slices_count = 0
            self._set_data_modified( )

    def spline_baseline_correction(self, data, probe_wn):
        baseline = np.zeros_like(data, dtype = float)
        for i in range(data.shape [ 0 ]):
            if len(probe_wn) >= 2:
                try:
                    valid_mask = np.isfinite(data [ i, : ])
                    if np.sum(valid_mask) >= 2:
                        spline = UnivariateSpline(probe_wn [ valid_mask ], data [ i, valid_mask ])
                        baseline [ i, : ] = spline(probe_wn)
                    else:
                        print(f"Not enough valid data points for spline in row {i}. Skipping spline for this row.")
                        baseline [ i, : ] = 0.0
                except Exception as e:
                    print(f"Error during spline calculation for row {i}: {e}")
                    baseline [ i, : ] = 0.0
            else:
                print(f"Not enough data points for spline in row {i}. Skipping spline for this row.")
                baseline [ i, : ] = 0.0
        corrected_data = data - baseline
        return corrected_data

    def _toggle_spline_correction(self):
        if not self.data_loaded:
            QMessageBox.warning(self, "No Data", "Please import data before applying spline baseline.")
            return
        current_min_level = float(self.min_level_input.text( )) if self.min_level_input.text( ) else None
        current_max_level = float(self.max_level_input.text( )) if self.max_level_input.text( ) else None
        current_signal_xlim = None
        current_signal_ylim = None
        if self.signal_plot_widget.plotItem.vb:
            signal_view_range = self.signal_plot_widget.plotItem.vb.viewRange( )
            if signal_view_range and len(signal_view_range) == 2:
                current_signal_xlim = signal_view_range [ 0 ]
                current_signal_ylim = signal_view_range [ 1 ]
        current_x_slice_xlim = None
        current_x_slice_ylim = None
        if self.x_slice_plot_widget.plotItem.vb:
            x_slice_view_range = self.x_slice_plot_widget.plotItem.vb.viewRange( )
            if x_slice_view_range and len(x_slice_view_range) == 2:
                current_x_slice_xlim = x_slice_view_range [ 0 ]
                current_x_slice_ylim = x_slice_view_range [ 1 ]
        current_y_slice_xlim = None
        current_y_slice_ylim = None
        if self.y_slice_plot_widget.plotItem.vb:
            y_slice_view_range = self.y_slice_plot_widget.plotItem.vb.viewRange( )
            if y_slice_view_range and len(y_slice_view_range) == 2:
                current_y_slice_xlim = y_slice_view_range [ 0 ]
                current_y_slice_ylim = y_slice_view_range [ 1 ]
        if not self.is_spline_corrected:
            try:
                corrected_data = self.spline_baseline_correction(
                    self._initial_raw_signal_data, self._initial_raw_x_values
                )
                self.current_signal_data = corrected_data
                self.is_spline_corrected = True
                self._set_data_modified( )
                self._refresh_all_plots(preserve_contour_levels = True,
                                        min_level = current_min_level,
                                        max_level = current_max_level,
                                        preserve_plot_ranges = True,
                                        signal_xlim = current_signal_xlim,
                                        signal_ylim = current_signal_ylim,
                                        x_slice_xlim = current_x_slice_xlim,
                                        x_slice_ylim = current_x_slice_ylim,
                                        y_slice_xlim = current_y_slice_xlim,
                                        y_slice_ylim = current_y_slice_ylim)
            except Exception as e:
                print(f"Failed to apply spline baseline: {e}")
                self.current_signal_data = self._initial_raw_signal_data.copy( )
                self.is_spline_corrected = False
                self._refresh_all_plots( )
        else:
            self.current_signal_data = self._initial_raw_signal_data.copy( )
            self.is_spline_corrected = False
            self._set_data_modified( )
            self._refresh_all_plots(preserve_contour_levels = True,
                                    min_level = current_min_level,
                                    max_level = current_max_level,
                                    preserve_plot_ranges = True,
                                    signal_xlim = current_signal_xlim,
                                    signal_ylim = current_signal_ylim,
                                    x_slice_xlim = current_x_slice_xlim,
                                    x_slice_ylim = current_x_slice_ylim,
                                    y_slice_xlim = current_y_slice_xlim,
                                    y_slice_ylim = current_y_slice_ylim)
        self._update_spline_button_text( )

    def _update_spline_button_text(self):
        if self.is_spline_corrected:
            self.spline_baseline_button.setText("Revert to Original")
        else:
            self.spline_baseline_button.setText("Use Spline Baseline")

    def _show_edit_names_dialog(self):
        dialog = EditNamesDialog(self, current_labels = self.axis_labels)
        if dialog.exec_( ) == QDialog.Accepted:
            new_labels = dialog.get_names( )
            if any(self.axis_labels [ key ] != new_labels [ key ] for key in self.axis_labels):
                self.axis_labels.update(new_labels)
                self._apply_axis_labels( )
                self._set_data_modified( )

    def _apply_axis_labels(self):
        self.signal_plot_widget.setLabel('bottom', self.axis_labels [ 'signal_bottom' ],
                                         **{ 'font-size': f'{self.axis_label_font_size}pt' })
        self.signal_plot_widget.setLabel('left', self.axis_labels [ 'signal_left' ],
                                         **{ 'font-size': f'{self.axis_label_font_size}pt' })
        self.x_slice_plot_widget.setLabel('bottom', self.axis_labels [ 'x_slice_bottom' ],
                                          **{ 'font-size': f'{self.axis_label_font_size}pt' })
        self.x_slice_plot_widget.setLabel('left', self.axis_labels [ 'x_slice_left' ],
                                          **{ 'font-size': f'{self.axis_label_font_size}pt' })
        self.y_slice_plot_widget.setLabel('bottom', self.axis_labels [ 'y_slice_bottom' ],
                                          **{ 'font-size': f'{self.axis_label_font_size}pt' })
        self.y_slice_plot_widget.setLabel('left', self.axis_labels [ 'y_slice_left' ],
                                          **{ 'font-size': f'{self.axis_label_font_size}pt' })

    def on_import_data_action_triggered(self):
        if self._data_modified:
            reply = QMessageBox.question(self, 'Save Changes',
                                         "You have unsaved changes. Do you want to save them before importing new data?",
                                         QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel,
                                         QMessageBox.Save)
            if reply == QMessageBox.Save:
                save_successful = self._save_project( )
                if not save_successful:
                    return
            elif reply == QMessageBox.Cancel:
                return
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Data File",
            "",
            "CSV Files (*.csv);;Text Files (*.txt);;All Files (*)"
        )
        if file_path:
            try:
                df = pd.read_csv(file_path, header = None)
                if df.shape [ 0 ] < 2 or df.shape [ 1 ] < 2:
                    raise ValueError("Data file must have at least 2 rows and 2 columns for X, Y, Z extraction.")
                y_values = df.iloc [ 1:, 0 ].values.astype(float)
                x_values = df.iloc [ 0, 1: ].values.astype(float)
                z_data = df.iloc [ 1:, 1: ].values.astype(float)
                self._load_data_into_plots(x_values, y_values, z_data)
                QMessageBox.information(self, "File Processed",
                                        f"File '{file_path}' loaded and data parsed successfully.\n"
                                        f"X-values shape: {x_values.shape}\n"
                                        f"Y-values shape: {y_values.shape}\n"
                                        f"Z-data shape: {z_data.shape}")
                print(f"Data file selected and parsed: {file_path}")
                self._current_project_file = None
                self._data_modified = True
                self._update_window_title( )
            except ValueError as ve:
                QMessageBox.critical(self, "Data Error", f"Error in data structure: {ve}")
                print(f"Error in data structure for {file_path}: {ve}")
                self.data_loaded = False
                self.update_plots( )
                self._data_modified = False
                self._update_window_title( )
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to process file: {e}")
                print(f"Error processing file {file_path}: {e}")
                self.data_loaded = False
                self.update_plots( )
                self._data_modified = False
                self._update_window_title( )
        else:
            QMessageBox.information(self, "Action", "No data file selected.")
            print("No data file selected.")

    def _open_x_fitter_tab(self):
        plot_data_item = self.x_slice_curve
        xlabel = self.axis_labels [ 'x_slice_bottom' ]
        ylabel = self.axis_labels [ 'x_slice_left' ]
        slice_axis_name = self._strip_html_tags(self.axis_labels [ 'signal_bottom' ])
        slice_value = self.current_x_values [ self.x_slider.value( ) ]
        slice_unit = self.x_unit_input.text( )
        is_spline_corrected = self.is_spline_corrected
        fitter_widget = exponential_fitter_wrapper(self, plot_data_item, xlabel, ylabel, slice_axis_name,
                                                   slice_value, slice_unit, is_spline_corrected)
        if fitter_widget:
            tab_index = self.tab_widget.addTab(fitter_widget, fitter_widget.objectName( ))
            self.tab_widget.setCurrentIndex(tab_index)
            self.active_fitter_tabs.append({ 'widget': fitter_widget, 'tab_index': tab_index })
            fitter_widget.destroyed.connect(partial(self._remove_fitter_tab, fitter_widget))
            self._set_data_modified( )

    def _open_fitter_tab(self, plot_data_item, is_x_slice, fitting_function_type, xlabel, ylabel, slice_axis_name,
                         slice_value, slice_unit, is_spline_corrected):
        fitter_widget = signal_fitter_wrapper(self, plot_data_item, is_x_slice, fitting_function_type, xlabel, ylabel,
                                              slice_axis_name, slice_value, slice_unit, is_spline_corrected)
        if fitter_widget:
            tab_index = self.tab_widget.addTab(fitter_widget, fitter_widget.objectName( ))
            self.tab_widget.setCurrentIndex(tab_index)
            self.active_fitter_tabs.append({ 'widget': fitter_widget, 'tab_index': tab_index })
            fitter_widget.destroyed.connect(
                partial(self._remove_fitter_tab, fitter_widget)
            )
            self._set_data_modified( )

    def _remove_fitter_tab(self, fitter_widget_to_remove):
        try:
            _ = self.objectName( )
        except RuntimeError as e:
            if "wrapped C/C++ object" in str(e):
                return

        for i, tab_info in enumerate(list(self.active_fitter_tabs)):
            if tab_info [ 'widget' ] == fitter_widget_to_remove:
                try:
                    index = self.tab_widget.indexOf(tab_info [ 'widget' ])
                    if index != -1:
                        self.tab_widget.removeTab(index)
                except RuntimeError as e:
                    if "wrapped C/C++ object" in str(e):
                        pass
                    else:
                        raise e

                self.active_fitter_tabs.pop(i)
                self._set_data_modified( )
                break

    def _save_project_data(self, file_path):
        project_state = {
            'data_loaded': self.data_loaded,
            'axis_labels': self.axis_labels,
            'x_legend_font_size': self.x_legend_font_size,
            'y_legend_font_size': self.y_legend_font_size,
            'current_slice_linewidth': self.current_slice_linewidth,
            'x_slice_legend_unit': self.x_unit_input.text( ),
            'y_slice_legend_unit': self.y_unit_input.text( ),
            'is_spline_corrected': self.is_spline_corrected,
            'current_interp_method': self._current_interp_method,
            'current_interp_multiplier': self._current_interp_multiplier,
            'active_fitter_tabs_states': [ ],
            'initial_raw_x_values': self._initial_raw_x_values.tolist( ) if self._initial_raw_x_values is not None else None,
            'initial_raw_y_values': self._initial_raw_y_values.tolist( ) if self._initial_raw_y_values is not None else None,
            'initial_raw_signal_data': self._initial_raw_signal_data.tolist( ) if self._initial_raw_signal_data is not None else None,
            'x_slider_value': self.x_slider.value( ),
            'y_slider_value': self.y_slider.value( ),
            'min_level_input': self.min_level_input.text( ),
            'max_level_input': self.max_level_input.text( ),
        }
        if self.data_loaded:
            project_state.update({
                'initial_raw_x_values': self._initial_raw_x_values.tolist( ),
                'initial_raw_y_values': self._initial_raw_y_values.tolist( ),
                'initial_raw_signal_data': self._initial_raw_signal_data.tolist( ),
                'x_slider_value': self.x_slider.value( ),
                'y_slider_value': self.y_slider.value( ),
                'min_level_input': self.min_level_input.text( ),
                'max_level_input': self.max_level_input.text( )
            })

        def save_plot_item(item):
            x, y = item.getData( )
            pen_color_rgb = (0.0, 0.0, 0.0, 1.0)
            pen_width = 1
            pen_style = str(Qt.SolidLine)
            if item.opts.get('pen') is not None:
                pen_obj = item.opts [ 'pen' ]
                pen_color_rgb = pen_obj.color( ).getRgbF( )
                pen_width = pen_obj.width( )
                pen_style = str(pen_obj.style( ))
            return {
                'x_data': x.tolist( ),
                'y_data': y.tolist( ),
                'pen_color_rgb': pen_color_rgb,
                'pen_width': pen_width,
                'pen_style': pen_style,
                'name': item.name( ),
                'z_value': item.zValue( )
            }

        project_state [ 'held_x_plots' ] = [
            save_plot_item(item) for item in self.x_slice_plot_widget.listDataItems( )
            if item != self.x_slice_curve
        ]
        project_state [ 'held_y_plots' ] = [
            save_plot_item(item) for item in self.y_slice_plot_widget.listDataItems( )
            if item != self.y_slice_curve
        ]
        for tab_info in self.active_fitter_tabs:
            widget = tab_info [ 'widget' ]
            if hasattr(widget, 'get_state'):
                project_state [ 'active_fitter_tabs_states' ].append(widget.get_state( ))
            elif hasattr(widget, 'get_fitter_state'):
                project_state [ 'active_fitter_tabs_states' ].append(widget.get_fitter_state( ))

        try:
            with open(file_path, 'w') as f:
                json.dump(project_state, f, indent = 4)
            QMessageBox.information(self, "Success", f"Project saved to {file_path}")
            self._current_project_file = file_path
            self._data_modified = False
            self._update_window_title( )
            return True
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Save failed: {str(e)}")
            return False

    def _save_project(self):
        if not self._current_project_file:
            return self._save_project_as( )
        else:
            return self._save_project_data(self._current_project_file)

    def _save_project_as(self):
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Project As", "", "Project Files (*.specdatpp);;JSON Files (*.json);;All Files (*)"
        )
        if not file_path:
            QMessageBox.information(self, "Save Project As", "Project save cancelled.")
            return False
        return self._save_project_data(file_path)

    def _load_project_from_path(self, file_path):
        try:
            with open(file_path, 'r') as f:
                project_state = json.load(f)
            self.x_unit_input.textChanged.disconnect(self._set_data_modified)
            self.y_unit_input.textChanged.disconnect(self._set_data_modified)
            self.min_level_input.textChanged.disconnect(self._set_data_modified)
            self.max_level_input.textChanged.disconnect(self._set_data_modified)
            self.clear_x_slice_plots( )
            self.clear_y_slice_plots( )
            for tab_info in list(self.active_fitter_tabs):
                self.tab_widget.removeTab(self.tab_widget.indexOf(tab_info [ 'widget' ]))
                tab_info [ 'widget' ].deleteLater( )
            self.active_fitter_tabs.clear( )
            if project_state.get('data_loaded', False):
                self._initial_raw_x_values = np.array(project_state [ 'initial_raw_x_values' ]) if project_state.get(
                    'initial_raw_x_values') is not None else None
                self._initial_raw_y_values = np.array(project_state [ 'initial_raw_y_values' ]) if project_state.get(
                    'initial_raw_y_values') is not None else None
                self._initial_raw_signal_data = np.array(
                    project_state [ 'initial_raw_signal_data' ]) if project_state.get(
                    'initial_raw_signal_data') is not None else None

                self.current_x_values = self._initial_raw_x_values.copy( )
                self.current_y_values = self._initial_raw_y_values.copy( )
                self.current_signal_data = self._initial_raw_signal_data.copy( )
                self.data_loaded = True

                self.is_spline_corrected = project_state.get('is_spline_corrected', False)

                if self.is_spline_corrected:
                    corrected_data = self.spline_baseline_correction(
                        self.current_signal_data, self.current_x_values
                    )
                    self.current_signal_data = corrected_data
                    print("Spline correction re-applied during load.")

                self._refresh_all_plots( )

                self.x_slider.setValue(project_state.get('x_slider_value', self.x_dim // 2 if self.x_dim else 0))
                self.y_slider.setValue(project_state.get('y_slider_value', self.y_dim // 2 if self.y_dim else 0))
                self.min_level_input.setText(project_state.get('min_level_input', ''))
                self.max_level_input.setText(project_state.get('max_level_input', ''))
            else:
                self.data_loaded = False
                self.update_plots( )
            self.axis_labels.update(project_state.get('axis_labels', { }))
            self._apply_axis_labels( )
            self.x_legend_font_size = project_state.get('x_legend_font_size', 14)
            self.y_legend_font_size = project_state.get('y_legend_font_size', 14)
            self._apply_legend_font_size(self.x_slice_legend, self.x_legend_font_size)
            self._apply_legend_font_size(self.y_slice_legend, self.y_legend_font_size)
            self.current_slice_linewidth = project_state.get('current_slice_linewidth', 2)
            self.x_slice_curve.setPen(pg.mkPen('b', width = self.current_slice_linewidth))
            self.y_slice_curve.setPen(pg.mkPen('r', width = self.current_slice_linewidth))
            self.x_slice_legend_unit = project_state.get('x_slice_legend_unit', "cm^-1")
            self.y_slice_legend_unit = project_state.get('y_slice_legend_unit', "ps")
            self.x_unit_input.setText(self.x_slice_legend_unit)
            self.y_unit_input.setText(self.y_slice_legend_unit)
            self._current_interp_method = project_state.get('current_interp_method', "None")
            self._current_interp_multiplier = project_state.get('current_interp_multiplier', 1)
            self.interp_method_combo.setCurrentText(self._current_interp_method)
            self.interp_multiplier_combo.setCurrentText(f"x{self._current_interp_multiplier}")
            self.held_x_slices_count = 0
            for plot_data in project_state.get('held_x_plots', [ ]):
                x_data = np.array(plot_data [ 'x_data' ])
                y_data = np.array(plot_data [ 'y_data' ])
                pen_color = QColor.fromRgbF(*plot_data [ 'pen_color_rgb' ])
                name = plot_data [ 'name' ]
                pen_width = plot_data.get('pen_width', self.current_slice_linewidth)
                pen_style_str = plot_data.get('pen_style', str(Qt.SolidLine))
                pen_style = Qt.SolidLine
                if pen_style_str == str(Qt.DotLine):
                    pen_style = Qt.DotLine
                elif pen_style_str == str(Qt.DashLine):
                    pen_style = Qt.DashLine
                elif pen_style_str == str(Qt.DashDotLine):
                    pen_style = Qt.DashDotLine
                elif pen_style_str == str(Qt.SolidLine):
                    pen_style = Qt.SolidLine
                pen = pg.mkPen(pen_color, width = pen_width, style = pen_style)
                self.x_slice_plot_widget.plot(x_data, y_data, pen = pen, name = name)
                self.held_x_slices_count += 1
            self.held_y_slices_count = 0
            for plot_data in project_state.get('held_y_plots', [ ]):
                x_data = np.array(plot_data [ 'x_data' ])
                y_data = np.array(plot_data [ 'y_data' ])
                pen_color = QColor.fromRgbF(*plot_data [ 'pen_color_rgb' ])
                name = plot_data [ 'name' ]
                pen_width = plot_data.get('pen_width', self.current_slice_linewidth)
                pen_style_str = plot_data.get('pen_style', str(Qt.SolidLine))
                pen_style = Qt.SolidLine
                if pen_style_str == str(Qt.DotLine):
                    pen_style = Qt.DotLine
                elif pen_style_str == str(Qt.DashLine):
                    pen_style = Qt.DashLine
                elif pen_style_str == str(Qt.DashDotLine):
                    pen_style = Qt.DashDotLine
                elif pen_style_str == str(Qt.SolidLine):
                    pen_style = Qt.SolidLine
                pen = pg.mkPen(pen_color, width = pen_width, style = pen_style)
                self.y_slice_plot_widget.plot(x_data, y_data, pen = pen, name = name)
                self.held_y_slices_count += 1
            fitter_states = project_state.get('active_fitter_tabs_states', [ ])
            for state in fitter_states:
                widget_type = state.get('type')
                new_widget = None
                if widget_type == 'GlobalFitApp':
                    new_widget = GlobalFitApp(parent = self)
                elif widget_type == 'GaussianFitterApp':
                    new_widget = GaussianFitterApp(parent = self)
                elif widget_type == 'ExponentialFitterApp':
                    new_widget = ExponentialFitterApp(parent = self)
                elif widget_type == 'PFIDFitterApp':
                    new_widget = PFIDFitterApp(main_window = self)

                if new_widget:
                    if hasattr(new_widget, 'set_fitter_state'):
                        new_widget.set_fitter_state(state)
                    else:
                        new_widget.set_state(state)
                    tab_name = state.get('tab_title', 'Restored Tab')
                    tab_index = self.tab_widget.addTab(new_widget, tab_name)
                    self.active_fitter_tabs.append({ 'widget': new_widget, 'tab_index': tab_index })
                    new_widget.destroyed.connect(partial(self._remove_fitter_tab, new_widget))

            QMessageBox.information(self, "Load Project", f"Project loaded successfully from '{file_path}'")
            self._current_project_file = file_path
            self._data_modified = False
            self._update_window_title( )
            return True
        except json.JSONDecodeError as jde:
            QMessageBox.critical(self, "Load Error", f"Failed to parse JSON file. Invalid format: {jde}")
            self._data_modified = False
            self._update_window_title( )
            return False
        except KeyError as ke:

            QMessageBox.critical(self, "Load Error",
                                 f"Missing data in project file: {ke}. File might be corrupted or from an incompatible version.")
            self._data_modified = False
            self._update_window_title( )
            return False
        except Exception as e:
            QMessageBox.critical(self, "Load Error", f"Failed to load project: {e}")
            self._data_modified = False
            self._update_window_title( )
            return False
        finally:
            self.x_unit_input.textChanged.connect(self._set_data_modified)
            self.y_unit_input.textChanged.connect(self._set_data_modified)
            self.min_level_input.textChanged.connect(self._set_data_modified)
            self.max_level_input.textChanged.connect(self._set_data_modified)

    def _load_project(self):
        if self._data_modified:
            reply = QMessageBox.question(self, 'Save Changes',
                                         "You have unsaved changes. Do you want to save them before loading a new project?",
                                         QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel,
                                         QMessageBox.Save)
            if reply == QMessageBox.Save:
                save_successful = self._save_project( )
                if not save_successful:
                    return
            elif reply == QMessageBox.Cancel:
                return
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Load Project",
            "",
            "Project Files (*.specdatpp);;JSON Files (*.json);;All Files (*)"
        )
        if file_path:
            self._load_project_from_path(file_path)
        else:
            QMessageBox.information(self, "Load Project", "Project load cancelled.")

    def closeEvent(self, event):
        """
        Overrides the standard close event to ensure the entire application,
        including the main event loop, terminates cleanly when the user clicks 'X'.
        This fixes the 'zombie process' bug by ensuring the application releases
        the lock on the executable file.
        """
        if self._data_modified:
            reply = QMessageBox.question(self, 'Save Changes',
                                         "You have unsaved changes. Do you want to save them before quitting?",
                                         QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel,
                                         QMessageBox.Save)
            if reply == QMessageBox.Save:
                save_successful = self._save_project( )
                if not save_successful:
                    event.ignore( )
                    return
            elif reply == QMessageBox.Cancel:
                event.ignore( )
                return

        event.accept( )

        QApplication.instance( ).quit( )


if __name__ == '__main__':
    if sys.platform == 'win32':
        try:
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID('dataviewer2D.Kaalen_app.1.0')
        except AttributeError:
            pass
    app = QApplication(sys.argv)
    app.setStyleSheet("QWidget { font-size: 10pt; }")
    app.setWindowIcon(QIcon(':/icons/icon.ico'))
    window = SignalPlotterApp( )
    window.show( )
    if len(sys.argv) > 1:
        file_to_open = sys.argv [ 1 ]
        if os.path.exists(file_to_open) and file_to_open.lower( ).endswith(('.specdatpp', '.json')):
            try:
                window._load_project_from_path(file_to_open)
            except Exception as e:
                QMessageBox.critical(window, "Error Opening Project",
                                     f"Failed to load project from '{file_to_open}': {e}")
        else:
            QMessageBox.warning(window, "Unsupported File",
                                f"The file '{file_to_open}' is not a recognized project file (.specdatpp or .json).")

    sys.exit(app.exec_( ))
