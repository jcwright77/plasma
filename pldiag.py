"""
TORIC Diagnostic Plotting Tool - Complete Version
Converted from IDL pldiag.pro to Python

This module reads diagnostic data from toric.asc files and creates
comprehensive plots of all plasma configuration and wave properties.

USAGE
python pldiag.py toric.asc --plots 4
python pldiag.py toric.asc --output eps --plots 2
python pldiag.py toric.asc --linetype 1 --plots 1 --output pdf

from pldiag import ToricDiagnostics
import matplotlib.pyplot as plt

# Configure plotter
diag = ToricDiagnostics(
    filename='toric.asc',
    output_type='screen',
    linetype=0,
    plots_per_page=2
)

# Run analysis
diag.run()

# All plots are now displayed or saved

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import sys
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime

class ToricDiagnostics:
    """Class to handle reading and plotting TORIC diagnostic data."""
    
    def __init__(self, filename: str, output_type: str = 'screen', 
                 linetype: int = 0, plots_per_page: int = 4):
        """
        Initialize the diagnostics plotter.
        
        Parameters
        ----------
        filename : str
            Path to the toric.asc file
        output_type : str
            'screen' (default), 'eps' (Encapsulated PostScript), 'pdf', or 'pdf-single'
        linetype : int
            0 for colored lines (default), 1 for solid/dotted lines
        plots_per_page : int
            Number of plots per page: 1, 2, or 4 (default: 4)
        """
        self.filename = Path(filename)
        self.output_type = output_type
        self.linetype = linetype
        self.plots_per_page = plots_per_page
        self.data = {}
        self.fig_counter = 0
        self.plot_counter = 0
        self.current_fig = None
        self.current_axes = None
        self.total_plots = 0

        # PDF single file handling
        self.pdf_pages = None
        if output_type == 'pdf-single':
            pdf_filename = Path(str(self.filename).replace('.asc', '_complete.pdf'))
            self.pdf_pages = PdfPages(str(pdf_filename))
            print(f"Creating combined PDF: {pdf_filename}")
        
        # Determine grid layout
        if plots_per_page == 1:
            self.grid_shape = (1, 1)
        elif plots_per_page == 2:
            self.grid_shape = (2, 1)
        else:
            self.grid_shape = (2, 2)
        
        self.setup_plot_params()
        self.setup_colors()
        
    def setup_plot_params(self):
        """Set up matplotlib plotting parameters."""
        plt.rcParams['font.size'] = 8
        plt.rcParams['lines.linewidth'] = 1.2
        plt.rcParams['axes.linewidth'] = 0.8
        if self.output_type == 'screen':
            plt.rcParams['figure.figsize'] = (12, 9)
        else:
            plt.rcParams['figure.figsize'] = (8.5, 11)
    
    def setup_colors(self):
        """Set up color scheme."""
        if self.output_type == 'screen':
            self.color_raw = '#1f77b4'
            self.color_interp = '#ff7f0e'
            self.color_real = '#2ca02c'
            self.color_imag = '#d62728'
            self.linestyle_raw = '-'
            self.linestyle_interp = '--'
        else:
            self.color_raw = 'black'
            self.color_interp = 'gray'
            self.color_real = 'black'
            self.color_imag = 'gray'
            self.linestyle_raw = '-'
            self.linestyle_interp = '--'
    
    def read_ascii_file(self) -> bool:
        """Read data from toric.asc file."""
        try:
            with open(self.filename, 'r') as f:
                lines = f.readlines()
            self.lines = lines
            self.line_idx = 0
            return True
        except FileNotFoundError:
            print(f"Error: File {self.filename} not found")
            return False
    
    def read_line(self) -> str:
        """Read next non-empty line from file."""
        while self.line_idx < len(self.lines):
            line = self.lines[self.line_idx].strip()
            self.line_idx += 1
            if line:
                return line
        return ""
    
    def peek_line(self) -> str:
        """Peek at next line without consuming it."""
        idx = self.line_idx
        while idx < len(self.lines):
            line = self.lines[idx].strip()
            if line:
                return line
            idx += 1
        return ""
    
    def read_int_line(self, n: int) -> List[int]:
        """Read a line with integers."""
        line = self.read_line()
        values = []
        for x in line.split():
            try:
                values.append(int(x))
            except ValueError:
                continue
        return values[:n] if n > 0 else values
    
    def read_float_line(self, n: int) -> np.ndarray:
        """Read floats spanning multiple lines."""
        values = []
        while len(values) < n:
            line = self.read_line()
            if not line:
                break
            try:
                vals = [float(x) for x in line.split()]
                values.extend(vals)
            except ValueError:
                continue
        return np.array(values[:n])
    
    def read_float_array(self, n: int) -> np.ndarray:
        """Read array of n floats."""
        return self.read_float_line(n)
    
    def get_next_ax(self, title: str = ""):
        """Get next subplot axis, creating new figure if needed."""
        if self.current_fig is None or self.plot_counter >= self.plots_per_page:
            self.finalize_page()
            self.fig_counter += 1
            self.current_fig, self.current_axes = plt.subplots(
                self.grid_shape[0], self.grid_shape[1],
                figsize=(12, 9) if self.plots_per_page == 4 else (12, 5.5)
            )
            
            if self.plots_per_page == 1:
                self.current_axes = np.array([self.current_axes])
            elif isinstance(self.current_axes, np.ndarray) and self.current_axes.ndim == 2:
                self.current_axes = self.current_axes.flatten()
            elif not isinstance(self.current_axes, np.ndarray):
                self.current_axes = np.array([self.current_axes])
            
            self.plot_counter = 0
        
        ax = self.current_axes[self.plot_counter]
        if title:
            ax.set_title(title, fontsize=9, fontweight='bold')
        self.plot_counter += 1
        self.total_plots += 1
        
        return ax
    
    def finalize_page(self):
        """Finalize and save/show current page."""
        if self.current_fig is not None:
            # Hide unused subplots
            for idx in range(self.plot_counter, len(self.current_axes)):
                self.current_axes[idx].set_visible(False)
            
            plt.tight_layout()

            # Add page number if saving to single PDF
            if self.output_type == 'pdf-single':
                page_num = self.fig_counter
                self.current_fig.text(0.99, 0.01, f'Page {page_num}', 
                                ha='right', va='bottom', fontsize=7, alpha=0.5)
            
            self.save_or_show(self.current_fig)
            self.current_fig = None
            self.current_axes = None
            self.plot_counter = 0

    def save_or_show(self, fig=None):
        """Save or show figure."""
        if fig is None:
            fig = plt.gcf()
    
        if self.output_type == 'screen':
            plt.show()
        elif self.output_type == 'pdf':
            filename = f"pldiag_{self.fig_counter:02d}.pdf"
            fig.savefig(filename, format='pdf', dpi=150, bbox_inches='tight')
            print(f"  Saved {filename}")
            plt.close(fig)
        elif self.output_type == 'pdf-single':
            self.pdf_pages.savefig(fig, dpi=150, bbox_inches='tight')
            plt.close(fig)
        elif self.output_type == 'ps':
            filename = f"pldiag_{self.fig_counter:02d}.ps"
            fig.savefig(filename, format='ps', dpi=150, bbox_inches='tight')
            print(f"  Saved {filename}")
            plt.close(fig)
        elif self.output_type == 'eps':
            filename = f"pldiag_{self.fig_counter:02d}.eps"
            fig.savefig(filename, format='eps', dpi=150, bbox_inches='tight')
            print(f"  Saved {filename}")
            plt.close(fig)
    
    def plot_coefficient_pair(self, title: str, srad: np.ndarray, 
                             irad: np.ndarray, rawcf: np.ndarray, 
                             intcf: np.ndarray):
        """Plot raw and interpolated coefficients."""
        ax = self.get_next_ax(title=title)
        
        ax.plot(srad, rawcf, 'o', color=self.color_raw, markersize=3,
               linestyle=self.linestyle_raw, linewidth=0.8, alpha=0.6,
               label='Raw')
        ax.plot(irad, intcf, 's', color=self.color_interp, markersize=2,
               linestyle=self.linestyle_interp, linewidth=0.8, alpha=0.6,
               label='Interp')
        
        ax.set_xlabel('psi', fontsize=8)
        ax.set_ylabel('Coefficient', fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
        ax.tick_params(labelsize=7)
    
    def plot_real_imag_pair(self, title: str, x: np.ndarray, real: np.ndarray,
                           imag: np.ndarray):
        """Plot real and imaginary components."""
        ax = self.get_next_ax(title=title)
        
        ax.plot(x, real, 's-', color=self.color_real, markersize=2,
               linewidth=0.8, alpha=0.7, label='Real')
        ax.plot(x, imag, '^--', color=self.color_imag, markersize=2,
               linewidth=0.8, alpha=0.7, label='Imag')
        
        ax.set_xlabel('x (cm)', fontsize=8)
        ax.set_ylabel('Value', fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)
        ax.tick_params(labelsize=7)
    
    def plot_single_curve(self, title: str, x: np.ndarray, y: np.ndarray,
                         xlabel: str = 'x', ylabel: str = 'y'):
        """Plot single curve."""
        ax = self.get_next_ax(title=title)
        
        ax.plot(x, y, 'o-', color=self.color_real, markersize=3,
               linewidth=0.8, alpha=0.7)
        
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)
    
    def plot_multiline(self, title: str, x: np.ndarray, curves: Dict[str, np.ndarray]):
        """Plot multiple lines."""
        ax = self.get_next_ax(title=title)
        
        colors = plt.cm.tab10(np.linspace(0, 1, len(curves)))
        for (label, y), color in zip(curves.items(), colors):
            ax.plot(x, y, 'o-', color=color, markersize=2, linewidth=0.8,
                   alpha=0.7, label=label)
        
        ax.set_xlabel('x', fontsize=8)
        ax.set_ylabel('Value', fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7, loc='best')
        ax.tick_params(labelsize=7)
    
    def plot_interpolation_tests(self):
        """Plot all interpolation coefficients."""
        print("\n" + "=" * 70)
        print("INTERPOLATION COEFFICIENT TESTS")
        print("=" * 70)
        
        # Read test parameters
        line = self.read_line()
        iqtest_val = int(line.split()[0]) if line else 0
        print(f"iqtest = {iqtest_val}")
        
        plot_name = self.read_line()
        print(f"Section: {plot_name}")
        
        dims = self.read_int_line(5)
        igsmhd, iudsym, npsi, ipsi, modmhd = dims
        print(f"Grid: npsi={npsi}, ipsi={ipsi}, modmhd={modmhd}")
        
        srad = self.read_float_array(npsi)
        irad = self.read_float_array(ipsi)
        self.data['srad'] = srad
        self.data['irad'] = irad
        self.data['npsi'] = npsi
        self.data['ipsi'] = ipsi
        
        # xc(0)
        plot_name = self.read_line()
        rawcf = self.read_float_array(npsi)
        intcf = self.read_float_array(ipsi)
        self.plot_coefficient_pair(f"xc(0): {plot_name}", srad, irad, rawcf, intcf)
        
        # zc(0) if asymmetric
        if iudsym == 0:
            plot_name = self.read_line()
            rawcf = self.read_float_array(npsi)
            intcf = self.read_float_array(ipsi)
            self.plot_coefficient_pair(f"zc(0): {plot_name}", srad, irad, rawcf, intcf)
        
        # Mode-dependent coefficients
        for m in range(1, modmhd + 1):
            plot_name = self.read_line()
            rawcf = self.read_float_array(npsi)
            intcf = self.read_float_array(ipsi)
            self.plot_coefficient_pair(f"xc({m}): {plot_name}", srad, irad, rawcf, intcf)
            
            plot_name = self.read_line()
            rawcf = self.read_float_array(npsi)
            intcf = self.read_float_array(ipsi)
            self.plot_coefficient_pair(f"zs({m}): {plot_name}", srad, irad, rawcf, intcf)
            
            if iudsym == 0:
                plot_name = self.read_line()
                rawcf = self.read_float_array(npsi)
                intcf = self.read_float_array(ipsi)
                self.plot_coefficient_pair(f"xs({m}): {plot_name}", srad, irad, rawcf, intcf)
                
                plot_name = self.read_line()
                rawcf = self.read_float_array(npsi)
                intcf = self.read_float_array(ipsi)
                self.plot_coefficient_pair(f"zc({m}): {plot_name}", srad, irad, rawcf, intcf)
        
        print(f"Plotted {m*2 + (2 if iudsym==0 else 1)} interpolation coefficient plots")
        return iqtest_val == 0
    
    def plot_safety_factor_and_current(self):
        """Plot safety factor and current profiles."""
        print("\nSAFETY FACTOR AND CURRENT PROFILES")
        print("-" * 70)
        
        srad = self.data.get('srad', np.linspace(0, 1, 151))
        npsi = len(srad)
        
        plot_name = self.read_line()
        qq = self.read_float_array(npsi)
        self.plot_single_curve(f"Safety Factor: {plot_name}", srad, qq, 'psi', 'q')
        
        plot_name = self.read_line()
        aj = self.read_float_array(npsi)
        self.plot_single_curve(f"Current Density: {plot_name}", srad, aj, 'psi', 'A/cm²')
        
        plot_name = self.read_line()
        ai = self.read_float_array(npsi)
        self.plot_single_curve(f"Integrated Current: {plot_name}", srad, ai, 'psi', 'kA')
        
        plot_name = self.read_line()
        vi = self.read_float_array(npsi)
        plot_name = self.read_line()
        ai_area = self.read_float_array(npsi)
        
        self.plot_multiline(f"Volume and Area", srad, 
                           {'Volume': vi, 'Area': ai_area})
        
        print("Plotted 4 equilibrium profile plots")
    
    def plot_magnetic_configuration(self):
        """Plot magnetic configuration."""
        print("\nMAGNETIC CONFIGURATION")
        print("-" * 70)
        
        # Read singularities
        plot_name = self.read_line()
        if plot_name == "    0":
            print("Magnetic configuration not available (idlout=0)")
            return
        
        dims = self.read_int_line(5)
        ncy1, ncy2, nres, ncof, npvert = dims
        print(f"Config: ncy1={ncy1}, ncy2={ncy2}, nres={nres}, ncof={ncof}, npvert={npvert}")
        
        # Read singularity data
        singularities = {'ncy1': [], 'ncy2': [], 'nres': [], 'ncof': []}
        
        if ncy1 > 0:
            plot_name = self.read_line()
            for _ in range(ncy1):
                k, npv = self.read_int_line(2)
                if npv > 0:
                    aux_x = self.read_float_array(npv)
                    aux_z = self.read_float_array(npv)
                    singularities['ncy1'].append((aux_x, aux_z))
        
        if ncy2 > 0:
            plot_name = self.read_line()
            for _ in range(ncy2):
                k, npv = self.read_int_line(2)
                if npv > 0:
                    aux_x = self.read_float_array(npv)
                    aux_z = self.read_float_array(npv)
                    singularities['ncy2'].append((aux_x, aux_z))
        
        if nres > 0:
            plot_name = self.read_line()
            for _ in range(nres):
                k, npv = self.read_int_line(2)
                if npv > 0:
                    aux_x = self.read_float_array(npv)
                    aux_z = self.read_float_array(npv)
                    singularities['nres'].append((aux_x, aux_z))
        
        if ncof > 0:
            plot_name = self.read_line()
            for _ in range(ncof):
                k, npv = self.read_int_line(2)
                if npv > 0:
                    aux_x = self.read_float_array(npv)
                    aux_z = self.read_float_array(npv)
                    singularities['ncof'].append((aux_x, aux_z))
        
        # Read configuration
        plot_name = self.read_line()
        igsmhd, iqtest = self.read_int_line(2)
        ncopsi, jptheta = self.read_int_line(2)
        
        # Read and plot flux contours
        splotx = self.read_float_array(jptheta)
        sploty = self.read_float_array(jptheta)
        
        ax = self.get_next_ax(title="Flux Contours and theta lines")
        ax.plot(splotx, sploty, 'b-', linewidth=0.8)
        ax.set_aspect('equal')
        ax.set_xlabel('x (cm)', fontsize=8)
        ax.set_ylabel('z (cm)', fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)
        
        # Separatrix and other contours
        plot_name = self.read_line()
        splotx = self.read_float_array(jptheta)
        sploty = self.read_float_array(jptheta)
        
        #ax = self.get_next_ax(title="Flux Contours - Separatrix")
        ax.plot(splotx, sploty, 'r--', linewidth=0.8)
        #ax.set_xlabel('x (cm)', fontsize=8)
        #ax.set_ylabel('z (cm)', fontsize=8)
        #ax.grid(True, alpha=0.3)
        #ax.tick_params(labelsize=7)
        
        # Additional contours
        for i in range(ncopsi): # range(min(2, ncopsi)):
            splotx = self.read_float_array(jptheta)
            sploty = self.read_float_array(jptheta)
            
            #ax = self.get_next_ax(title=f"Flux Contours - {i+1}")
            ax.plot(splotx, sploty, 'g-', linewidth=0.8, alpha=0.7)
        #    ax.set_xlabel('x (cm)', fontsize=8)
        #    ax.set_ylabel('z (cm)', fontsize=8)
        #    ax.grid(True, alpha=0.3)
        #    ax.tick_params(labelsize=7)
        
        # Skip remaining contours
        #for i in range(2, ncopsi):
        #    splotx = self.read_float_array(jptheta)
        #    sploty = self.read_float_array(jptheta)
        
        # Plot theta lines
        plot_name = self.read_line()
        ntt, lpl = self.read_int_line(2)
        
        #ax = self.get_next_ax(title=f"Constant Theta Lines ({ntt} lines)")
        for k in range(ntt): #range(min(5, ntt)):
            ssx = self.read_float_array(lpl)
            ssy = self.read_float_array(lpl)
            ax.plot(ssx, ssy, linewidth=0.6, alpha=0.7)
        
        # Skip remaining theta lines
        #for k in range(5, nAtt):
        #    ssx = self.read_float_array(lpl)
        #    ssy = self.read_float_array(lpl)
        
        #ax.set_xlabel('x (cm)', fontsize=8)
        #ax.set_ylabel('z (cm)', fontsize=8)
        #ax.grid(True, alpha=0.3)
        #ax.tick_params(labelsize=7)

        for res in singularities.keys():
            if len(singularities[res])>0:
                lres=singularities[res][0]
                for i in range(len(lres)):
                    ax.plot(lres[i][0],lrespi[[1],'-.')
        print(f"Plotted magnetic configuration plots")
        
    def plot_metric_elements(self):
        """Plot metric elements on magnetic surfaces."""
        print("\nMETRIC ELEMENTS ON MAGNETIC SURFACES")
        print("-" * 70)
        
        plot_name = self.read_line()
        nsfs, ntt0 = self.read_int_line(2)
        ntt = ntt0 + 1
        
        plot_name = self.read_line()
        thdeg = self.read_float_array(ntt)
        
        metrics_names = ['Jpol', 'Ntau', 'Gpol', 'J_p/R', 'Gpol/N²', 'tanQ']
        plot_count = 0
        
        for i in range(nsfs):
            # Read 6 metric elements for this surface
            metric_data = {}
            for j, metric_name in enumerate(metrics_names):
                plot_name = self.read_line()
                tqaus = self.read_float_array(ntt)
                metric_data[metric_name] = tqaus
            if i==0:
                for metric_name, values in metric_data.items():
                    ax = self.get_next_ax(title=f"{metric_name} vs theta (surf {i})")
                    ax.plot(thdeg, values, 'o-', markersize=2, linewidth=0.8)
                    ax.set_xlabel('θ/2π', fontsize=8)
                    ax.set_ylabel(metric_name, fontsize=8)
                    ax.grid(True, alpha=0.3)
                    ax.tick_params(labelsize=7)
                    plot_count += 1
        
        print(f"Plotted {plot_count} metric element plots")
    
    def plot_equatorial_metrics(self):
        """Plot equatorial plane metrics."""
        print("\nEQUATORIAL PLANE METRICS")
        print("-" * 70)
        
        plot_name = self.read_line()
        nupx = int(self.read_line())
        
        plot_name = self.read_line()
        stbx = self.read_float_array(nupx)
        plot_name = self.read_line()
        stbz = self.read_float_array(nupx)
        
        # Psi mapping
        plot_name = self.read_line()
        nmhd = int(self.read_line())
        
        plot_name = self.read_line()
        psipol = self.read_float_array(nmhd)
        plot_name = self.read_line()
        psitor = self.read_float_array(nmhd)
        plot_name = self.read_line()
        xovera = self.read_float_array(nmhd)
        
        # Plot psi mappings
        ax = self.get_next_ax(title="Mapping psi_pol to X")
        ax.plot(psipol, xovera, 'b-', linewidth=0.8, label='X/a (outer)')
        ax.plot(psipol, psipol, 'r--', linewidth=0.8, label='psi')
        ax.set_xlabel('ψ_pol', fontsize=8)
        ax.set_ylabel('Value', fontsize=8)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)
        
        ax = self.get_next_ax(title="Mapping psi_tor to X")
        ax.plot(psitor, xovera, 'b-', linewidth=0.8, label='X/a (outer)')
        ax.plot(psitor, psitor, 'r--', linewidth=0.8, label='psi')
        ax.set_xlabel('ψ_tor', fontsize=8)
        ax.set_ylabel('Value', fontsize=8)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)
        
        # Plot equatorial metrics
        plot_count = 2
        for ipl in range(7):
            plot_name = self.read_line()
            stby = self.read_float_array(nupx)
            
            ax = self.get_next_ax(title=f"{plot_name}")
            ax.plot(stbz, stby, 's-', markersize=2, linewidth=0.8, color='blue')
            ax.set_xlabel('X (cm)', fontsize=8)
            ax.set_ylabel('Value', fontsize=8)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
        
        print(f"Plotted {plot_count} equatorial metric plots")
    
    def plot_density_temperature(self):
        """Plot density and temperature profiles."""
        print("\nDENSITY AND TEMPERATURE PROFILES")
        print("-" * 70)
        
        plot_name = self.read_line()
        npxeq_line = self.read_line()
        npxeq = int(npxeq_line)
        
        plot_name = self.read_line()
        xeqpl = self.read_float_array(npxeq)
        
        # Density
        plot_name = self.read_line()
        tbden = self.read_float_array(npxeq)
        
        ax = self.get_next_ax(title="Density Profile")
        ax.plot(xeqpl, tbden, 'b-', linewidth=0.8, marker='o', markersize=2)
        ax.set_xlabel('x (cm)', fontsize=8)
        ax.set_ylabel('n (cm⁻³)', fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)
        
        # Temperature
        plot_name = self.read_line()
        tbte = self.read_float_array(npxeq)
        plot_name = self.read_line()
        tbti = self.read_float_array(npxeq)
        
        ax = self.get_next_ax(title="Temperature Profile")
        ax.plot(xeqpl, tbte, 'r-', linewidth=0.8, marker='o', markersize=2, label='Te')
        ax.plot(xeqpl, tbti, 'b-', linewidth=0.8, marker='s', markersize=2, label='Ti')
        ax.set_xlabel('x (cm)', fontsize=8)
        ax.set_ylabel('T (keV)', fontsize=8)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)
        
        self.data['xeqpl'] = xeqpl
        self.data['npxeq'] = npxeq
        
        print("Plotted 2 density/temperature plots")
        return npxeq
    
    def plot_dielectric_tensor(self, npxeq: int):
        """Plot dielectric tensor components."""
        print("\nDIELECTRIC TENSOR COMPONENTS")
        print("-" * 70)
        
        xeqpl = self.data.get('xeqpl', np.linspace(0, 1, npxeq))
        plot_count = 0
        
        # R (real, 1 curve)
        plot_name = self.read_line()
        auxr = self.read_float_array(npxeq)
        
        ax = self.get_next_ax(title=f"R (Zero Larmor Radius)")
        ax.plot(xeqpl, auxr, 'b-', linewidth=0.8, marker='o', markersize=2)
        ax.set_xlabel('x (cm)', fontsize=8)
        ax.set_ylabel('R', fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)
        plot_count += 1
        
        # L, S, P components (3 curves each with real and imaginary)
        components = ['L', 'S', 'P']
        for comp_name in components:
            plot_name = self.read_line()
            real_data = []
            imag_data = []
            
            for k in range(3):
                auxr = self.read_float_array(npxeq)
                auxi = self.read_float_array(npxeq)
                real_data.append(auxr)
                imag_data.append(auxi)
            
            # Plot real parts
            ax = self.get_next_ax(title=f"{comp_name} - Real Parts")
            colors = ['b', 'g', 'r']
            for k, (data, color) in enumerate(zip(real_data, colors)):
                ax.plot(xeqpl, data, linestyle='-', color=color, linewidth=0.8,
                       marker='o', markersize=1, alpha=0.7, label=f'Mode {k}')
            ax.set_xlabel('x (cm)', fontsize=8)
            ax.set_ylabel(f'{comp_name} (Real)', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
            
            # Plot imaginary parts
            ax = self.get_next_ax(title=f"{comp_name} - Imaginary Parts")
            for k, (data, color) in enumerate(zip(imag_data, colors)):
                ax.plot(xeqpl, data, linestyle='--', color=color, linewidth=0.8,
                       marker='s', markersize=1, alpha=0.7, label=f'Mode {k}')
            ax.set_xlabel('x (cm)', fontsize=8)
            ax.set_ylabel(f'{comp_name} (Imag)', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
        
        # Finite Larmor Radius elements
        # rho_i (real, 1 curve)
        plot_name = self.read_line()
        auxr = self.read_float_array(npxeq)
        
        ax = self.get_next_ax(title="ρ_i (Larmor Radius)")
        ax.plot(xeqpl, auxr, 'b-', linewidth=0.8, marker='o', markersize=2)
        ax.set_xlabel('x (cm)', fontsize=8)
        ax.set_ylabel('ρ_i', fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)
        plot_count += 1
        
        # lambda_i, lambda_e, xi_e (each with real and imaginary)
        flr_names = ['λ_i', 'λ_e', 'ξ_e']
        for flr_name in flr_names:
            plot_name = self.read_line()
            real_data = []
            imag_data = []
            
            for k in range(3):
                auxr = self.read_float_array(npxeq)
                auxi = self.read_float_array(npxeq)
                real_data.append(auxr)
                imag_data.append(auxi)
            
            # Plot real parts
            ax = self.get_next_ax(title=f"{flr_name} - Real Parts")
            colors = ['b', 'g', 'r']
            for k, (data, color) in enumerate(zip(real_data, colors)):
                ax.plot(xeqpl, data, linestyle='-', color=color, linewidth=0.8,
                       marker='o', markersize=1, alpha=0.7, label=f'Mode {k}')
            ax.set_xlabel('x (cm)', fontsize=8)
            ax.set_ylabel(f'{flr_name} (Real)', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
            
            # Plot imaginary parts
            ax = self.get_next_ax(title=f"{flr_name} - Imaginary Parts")
            for k, (data, color) in enumerate(zip(imag_data, colors)):
                ax.plot(xeqpl, data, linestyle='--', color=color, linewidth=0.8,
                       marker='s', markersize=1, alpha=0.7, label=f'Mode {k}')
            ax.set_xlabel('x (cm)', fontsize=8)
            ax.set_ylabel(f'{flr_name} (Imag)', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
        
        print(f"Plotted {plot_count} dielectric tensor plots")
        return plot_count
    
    def plot_dispersion_indices(self, npxeq: int):
        """Plot wave mode dispersion indices."""
        print("\nDISPERSION RELATION INDICES")
        print("-" * 70)
        
        xeqpl = self.data.get('xeqpl', np.linspace(0, 1, npxeq))
        plot_count = 0
        
        # Fast Wave indices
        wave_modes = [
            ('Fast Wave (Full FLR)', 'N²_FW Full'),
            ('Fast Wave k²ρ_i/2', 'k²ρ_i/2 FW'),
            ('Fast Wave (Approx)', 'N²_FW Approx'),
            ('Fast Wave (Hot)', 'N²_FW Hot'),
            ('IBW (Full FLR)', 'N²_IBW Full'),
            ('IBW k²ρ_i/2', 'k²ρ_i/2 IBW'),
            ('IBW (Approx)', 'N²_IBW Approx'),
            ('IBW (Hot)', 'N²_IBW Hot'),
            ('Shear Alfvén (Full)', 'N²_SA Full'),
            ('Shear Alfvén (Approx)', 'N²_SA Approx'),
            ('Shear Alfvén (Hot)', 'N²_SA Hot'),
        ]
        
        for mode_label, plot_label in wave_modes:
            plot_name = self.read_line()
            if not plot_name:
                break
            
            # Check if this is k²ρ_i/2 type (only real, 3 curves)
            if 'k²ρ' in mode_label or 'k²ρ' in plot_label:
                real_data = []
                for k in range(3):
                    auxr = self.read_float_array(npxeq)
                    real_data.append(auxr)
                
                ax = self.get_next_ax(title=mode_label)
                colors = ['b', 'g', 'r']
                for k, (data, color) in enumerate(zip(real_data, colors)):
                    ax.plot(xeqpl, data, linestyle='-', color=color, linewidth=0.8,
                           marker='o', markersize=1, alpha=0.7, label=f'Mode {k}')
                ax.set_xlabel('x (cm)', fontsize=8)
                ax.set_ylabel('Value', fontsize=8)
                ax.legend(fontsize=7)
                ax.grid(True, alpha=0.3)
                ax.tick_params(labelsize=7)
                plot_count += 1
            else:
                # Real and imaginary parts
                real_data = []
                imag_data = []
                
                for k in range(3):
                    auxr = self.read_float_array(npxeq)
                    auxi = self.read_float_array(npxeq)
                    real_data.append(auxr)
                    imag_data.append(auxi)
                
                # Plot real parts
                ax = self.get_next_ax(title=f"{mode_label} - Real")
                colors = ['b', 'g', 'r']
                for k, (data, color) in enumerate(zip(real_data, colors)):
                    ax.plot(xeqpl, data, linestyle='-', color=color, linewidth=0.8,
                           marker='o', markersize=1, alpha=0.7, label=f'Mode {k}')
                ax.set_xlabel('x (cm)', fontsize=8)
                ax.set_ylabel('Real(N²)', fontsize=8)
                ax.legend(fontsize=7)
                ax.grid(True, alpha=0.3)
                ax.tick_params(labelsize=7)
                plot_count += 1
                
                # Plot imaginary parts
                ax = self.get_next_ax(title=f"{mode_label} - Imag")
                for k, (data, color) in enumerate(zip(imag_data, colors)):
                    ax.plot(xeqpl, data, linestyle='--', color=color, linewidth=0.8,
                           marker='s', markersize=1, alpha=0.7, label=f'Mode {k}')
                ax.set_xlabel('x (cm)', fontsize=8)
                ax.set_ylabel('Imag(N²)', fontsize=8)
                ax.legend(fontsize=7)
                ax.grid(True, alpha=0.3)
                ax.tick_params(labelsize=7)
                plot_count += 1
        
        print(f"Plotted {plot_count} dispersion index plots")
        return plot_count
    
    def plot_parallel_index_and_damping(self, npxeq: int):
        """Plot parallel index, phase velocity, and damping."""
        print("\nPARALLEL INDEX AND DAMPING")
        print("-" * 70)
        
        xeqpl = self.data.get('xeqpl', np.linspace(0, 1, npxeq))
        plot_count = 0
        
        # Parallel index
        plot_name = self.read_line()
        parallel_data = []
        for k in range(3):
            auxr = self.read_float_array(npxeq)
            parallel_data.append(auxr)
        
        ax = self.get_next_ax(title="Parallel Index")
        colors = ['b', 'g', 'r']
        for k, (data, color) in enumerate(zip(parallel_data, colors)):
            ax.plot(xeqpl, data, linestyle='-', color=color, linewidth=0.8,
                   marker='o', markersize=1, alpha=0.7, label=f'Mode {k}')
        ax.set_xlabel('x (cm)', fontsize=8)
        ax.set_ylabel(r'$N_\|$', fontsize=8)
        ax.set_ylim( -100,100 ) #JCW
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)
        plot_count += 1
        
        # Parallel phase velocity
        plot_name = self.read_line()
        phase_vel_data = []
        for k in range(3):
            auxr = self.read_float_array(npxeq)
            phase_vel_data.append(auxr)
        
        ax = self.get_next_ax(title="Normalized Phase Velocity")
        for k, (data, color) in enumerate(zip(phase_vel_data, colors)):
            ax.plot(xeqpl, data, linestyle='-', color=color, linewidth=0.8,
                   marker='s', markersize=1, alpha=0.7, label=f'Mode {k}')
        ax.set_xlabel('x (cm)', fontsize=8)
        ax.set_ylabel('ω/(k_z·v_th)', fontsize=8)
        ax.set_ylim( -1,10 ) #JCW
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)
        plot_count += 1
        
        # Electron Landau damping
        ibweld_line = self.read_line()
        try:
            ibweld = int(ibweld_line)
        except:
            ibweld = 0
        
        if ibweld > 0:
            plot_name = self.read_line()
            damping_data = []
            for k in range(3):
                auxr = self.read_float_array(npxeq)
                damping_data.append(auxr)
            
            ax = self.get_next_ax(title="Electron Landau Damping (IBW)")
            for k, (data, color) in enumerate(zip(damping_data, colors)):
                ax.plot(xeqpl, data, linestyle='-', color=color, linewidth=0.8,
                       marker='d', markersize=1, alpha=0.7, label=f'Mode {k}')
            ax.set_xlabel('x (cm)', fontsize=8)
            ax.set_ylabel('Damping', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
        
        print(f"Plotted {plot_count} parallel index/damping plots")
        return plot_count
    
    def plot_theta_dependent_quantities(self):
        """Plot theta-dependent quantities and Fourier transforms."""
        print("\nTHETA-DEPENDENT QUANTITIES AND FOURIER TRANSFORMS")
        print("-" * 70)
        
        iplth_line = self.read_line()
        try:
            iplth = int(iplth_line)
        except:
            iplth = 0
        
        if iplth == 0:
            print("No theta-dependent data available")
            return 0
        
        plot_count = 0
        
        for mp in range(1, 4):
            nnpl, npthet = self.read_int_line(2)
            npl = nnpl - 1
            
            print(f"  Processing mode set {mp}: {nnpl} poloidal positions, {npthet} theta points")
            
            # Read theta values
            plot_name = self.read_line()
            thval = self.read_float_array(npthet)
            
            # Create color array for different surfaces
            colors = plt.cm.tab10(np.linspace(0, 1, nnpl))
            
            # L vs theta
            plot_name = self.read_line()
            realpart = []
            imagpart = []
            for i in range(npl + 1):
                auxpr = self.read_float_array(npthet)
                auxpi = self.read_float_array(npthet)
                realpart.append(auxpr)
                imagpart.append(auxpi)
            
            ax = self.get_next_ax(title=f"L vs θ - Mode Set {mp} (Real)")
            for i, (rp, color) in enumerate(zip(realpart[:3], colors[:3])):
                ax.plot(thval, rp, linewidth=0.8, color=color, alpha=0.7, label=f'Surf {i}')
            ax.set_xlabel('θ/2π', fontsize=8)
            ax.set_ylabel('L (Real)', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
            
            ax = self.get_next_ax(title=f"L vs θ - Mode Set {mp} (Imag)")
            for i, (ip, color) in enumerate(zip(imagpart[:3], colors[:3])):
                ax.plot(thval, ip, linewidth=0.8, color=color, linestyle='--', alpha=0.7, label=f'Surf {i}')
            ax.set_xlabel('θ/2π', fontsize=8)
            ax.set_ylabel('L (Imag)', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
            
            # N²_FW vs theta
            plot_name = self.read_line()
            realpart = []
            imagpart = []
            for i in range(npl + 1):
                auxpr = self.read_float_array(npthet)
                auxpi = self.read_float_array(npthet)
                realpart.append(auxpr)
                imagpart.append(auxpi)
            
            ax = self.get_next_ax(title=f"N²_FW vs θ - Mode Set {mp} (Real)")
            for i, (rp, color) in enumerate(zip(realpart[:3], colors[:3])):
                ax.plot(thval, rp, linewidth=0.8, color=color, alpha=0.7, label=f'Surf {i}')
            ax.set_xlabel('θ/2π', fontsize=8)
            ax.set_ylabel('N²_FW (Real)', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
            
            # Lambda_i vs theta
            plot_name = self.read_line()
            realpart = []
            imagpart = []
            for i in range(npl + 1):
                auxpr = self.read_float_array(npthet)
                auxpi = self.read_float_array(npthet)
                realpart.append(auxpr)
                imagpart.append(auxpi)
            
            ax = self.get_next_ax(title=f"λ_i vs θ - Mode Set {mp}")
            for i, (rp, color) in enumerate(zip(realpart[:3], colors[:3])):
                ax.plot(thval, rp, linewidth=0.8, color=color, alpha=0.7, label=f'Surf {i}')
            ax.set_xlabel('θ/2π', fontsize=8)
            ax.set_ylabel('λ_i', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
            
            # N²_B vs theta
            plot_name = self.read_line()
            realpart = []
            imagpart = []
            for i in range(npl + 1):
                auxpr = self.read_float_array(npthet)
                auxpi = self.read_float_array(npthet)
                realpart.append(auxpr)
                imagpart.append(auxpi)
            
            ax = self.get_next_ax(title=f"N²_B vs θ - Mode Set {mp}")
            for i, (rp, color) in enumerate(zip(realpart[:3], colors[:3])):
                ax.plot(thval, rp, linewidth=0.8, color=color, alpha=0.7, label=f'Surf {i}')
            ax.set_xlabel('θ/2π', fontsize=8)
            ax.set_ylabel('N²_B', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
            
            # Fourier transform data
            plot_name = self.read_line()
            m_values = self.read_float_array(npthet)
            
            # L vs m
            plot_name = self.read_line()
            realpart_m = []
            imagpart_m = []
            for i in range(npl + 1):
                auxpr = self.read_float_array(npthet)
                auxpi = self.read_float_array(npthet)
                realpart_m.append(auxpr)
                imagpart_m.append(auxpi)
            
            ax = self.get_next_ax(title=f"L vs m - Mode Set {mp}")
            for i, (rp, color) in enumerate(zip(realpart_m[:3], colors[:3])):
                ax.plot(m_values, rp, 'o-', markersize=2, linewidth=0.8,
                       color=color, alpha=0.7, label=f'Surf {i}')
            ax.set_xlabel('m (poloidal mode)', fontsize=8)
            ax.set_ylabel('L', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
            
            # N²_FW vs m
            plot_name = self.read_line()
            realpart_m = []
            imagpart_m = []
            for i in range(npl + 1):
                auxpr = self.read_float_array(npthet)
                auxpi = self.read_float_array(npthet)
                realpart_m.append(auxpr)
                imagpart_m.append(auxpi)
            
            ax = self.get_next_ax(title=f"N²_FW vs m - Mode Set {mp}")
            for i, (rp, color) in enumerate(zip(realpart_m[:3], colors[:3])):
                ax.plot(m_values, rp, 'o-', markersize=2, linewidth=0.8,
                       color=color, alpha=0.7, label=f'Surf {i}')
            ax.set_xlabel('m (poloidal mode)', fontsize=8)
            ax.set_ylabel('N²_FW', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
            
            # Lambda_i vs m
            plot_name = self.read_line()
            realpart_m = []
            imagpart_m = []
            for i in range(npl + 1):
                auxpr = self.read_float_array(npthet)
                auxpi = self.read_float_array(npthet)
                realpart_m.append(auxpr)
                imagpart_m.append(auxpi)
            
            ax = self.get_next_ax(title=f"λ_i vs m - Mode Set {mp}")
            for i, (rp, color) in enumerate(zip(realpart_m[:3], colors[:3])):
                ax.plot(m_values, rp, 'o-', markersize=2, linewidth=0.8,
                       color=color, alpha=0.7, label=f'Surf {i}')
            ax.set_xlabel('m (poloidal mode)', fontsize=8)
            ax.set_ylabel('λ_i', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
            
            # N²_B vs m
            plot_name = self.read_line()
            realpart_m = []
            imagpart_m = []
            for i in range(npl + 1):
                auxpr = self.read_float_array(npthet)
                auxpi = self.read_float_array(npthet)
                realpart_m.append(auxpr)
                imagpart_m.append(auxpi)
            
            ax = self.get_next_ax(title=f"N²_B vs m - Mode Set {mp}")
            for i, (rp, color) in enumerate(zip(realpart_m[:3], colors[:3])):
                ax.plot(m_values, rp, 'o-', markersize=2, linewidth=0.8,
                       color=color, alpha=0.7, label=f'Surf {i}')
            ax.set_xlabel('m (poloidal mode)', fontsize=8)
            ax.set_ylabel('N²_B', fontsize=8)
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            plot_count += 1
        
        print(f"Plotted {plot_count} theta-dependent and Fourier plots")
        return plot_count

    def finalize_pdf(self):
        """Finalize and close combined PDF."""
        if self.pdf_pages is not None:
            # Add metadata
            d = self.pdf_pages.infodict()
            d['Title'] = 'TORIC Diagnostic Analysis'
            d['Author'] = 'TORIC Plotter'
            d['Subject'] = 'Plasma Diagnostic Plots'
            d['Keywords'] = 'TORIC, Plasma, Diagnostics'
            d['CreationDate'] = datetime.now()
        
            self.pdf_pages.close()
            print(f"\n✓ Combined PDF file created successfully")
    
    def run(self):
        """Run complete diagnostic analysis."""
        if not self.read_ascii_file():
            return False
        
        print(f"\n{'='*70}")
        print(f"TORIC COMPLETE DIAGNOSTIC PLOTTER")
        print(f"{'='*70}")
        print(f"File: {self.filename}")
        print(f"Output: {self.output_type}")
        print(f"Plots per page: {self.plots_per_page}")
        print(f"Layout: {self.grid_shape[0]}×{self.grid_shape[1]}")
        
        try:
            # Plot all data sections
            if not self.plot_interpolation_tests():
                print("iqtest != 0, stopping...")
                return False
            
            self.plot_safety_factor_and_current()
            self.plot_magnetic_configuration()
            self.plot_metric_elements()
            self.plot_equatorial_metrics()
            npxeq = self.plot_density_temperature()
            self.plot_dielectric_tensor(npxeq)
            self.plot_dispersion_indices(npxeq)
            self.plot_parallel_index_and_damping(npxeq)
            self.plot_theta_dependent_quantities()
            
            # Finalize last page
            self.finalize_page()
            
        except Exception as e:
            print(f"\nError during plotting: {e}")
            import traceback
            traceback.print_exc()
            return False

         # Finalize output
        if self.output_type == 'pdf-single':
            self.finalize_pdf()
        else:
            self.finalize_page()
        
        print(f"\n{'='*70}")
        print(f"PLOTTING COMPLETE")
        print(f"{'='*70}")
        print(f"Output type: {self.output_type}")
        print(f"Total figures: {self.fig_counter}")
        print(f"Total plots: {self.total_plots}")


        if self.output_type == 'pdf':
            print(f"Generated PDF files: pldiag_01.pdf through pldiag_{self.fig_counter:02d}.pdf")
        elif self.output_type == 'pdf-single':
            pdf_name = str(self.filename).replace('.asc', '_complete.pdf')
            print(f"Single combined PDF: {pdf_name}")
        elif self.output_type == 'ps':
            print(f"Generated PS files: pldiag_01.ps through pldiag_{self.fig_counter:02d}.ps")
        elif self.output_type == 'eps':
            print(f"Generated EPS files: pldiag_01.eps through pldiag_{self.fig_counter:02d}.eps")

        return True


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="TORIC Complete Diagnostic Plotter",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python pldiag.py toric.asc
  python pldiag.py toric.asc --output eps --plots 2
  python pldiag.py toric.asc --output pdf --plots 2
  python pldiag.py toric.asc --output pdf-single
        """
    )
    parser.add_argument('filename', nargs='?', default='toric.asc',
                       help='Input toric.asc file')
    parser.add_argument('--output', '-o',
                        choices=['screen', 'pdf', 'pdf-single', 'ps', 'eps'],
                       default='screen', help='Output type (default: screen)')
    parser.add_argument('--linetype', '-l', type=int, choices=[0, 1],
                       default=0, help='0: colored (default), 1: solid/dotted')
    parser.add_argument('--plots', '-p', type=int, choices=[1, 2, 4],
                       default=4, help='Plots per page (default: 4)')
    
    args = parser.parse_args()

    # Validate pdf-single option (should not use with --plots option for layout control)
    if args.output == 'pdf-single' and args.plots != 4:
        print(f"Note: pdf-single uses plots_per_page={args.plots}")
        
    diag = ToricDiagnostics(
        filename=args.filename,
        output_type=args.output,
        linetype=args.linetype,
        plots_per_page=args.plots
    )
    
    success = diag.run()
    return 0 if success else 1


if __name__ == '__main__':
    sys.exit(main())        
