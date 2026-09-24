import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog

from matplotlib.legend_handler import HandlerPatch
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.lines import Line2D
from cmcrameri import cm
import numpy as np
import tkinter.messagebox as messagebox

import rheolopy as re


class HandlerRect(HandlerPatch):
    def create_artists(
        self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
    ):
        rect = Rectangle(
            [xdescent, ydescent],
            width,
            height,
            facecolor=orig_handle.get_color(),
            edgecolor="none",
            transform=trans,
        )
        return [rect]


class RheologyExplorerApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Rheology Explorer")
        self.geometry("1400x1000")
        self.model = None  # Cache for the loaded background model
        self.bg_sigma_plot = None
        self.bg_z_plot_shifted = None
        self.bg_top_z = None
        self.config = re.load_config()

        # Accessing parameters from the General section
        self.strain_rate = self.config.getfloat("General", "strain_rate")
        self.geotherm = re.Geotherm()
        self.mat_dbase = re.materials()
        self.bg_artists = []  # Store background rectangles and labels
        self.bg_model_visible = False

        self.configure(bg="#2b2b2b")  # App background

        # Apply dark theme style to ttk
        style = ttk.Style(self)
        style.theme_use("default")

        style.configure("TFrame", background="#2b2b2b")
        style.configure(
            "TLabel", background="#2b2b2b", foreground="#e0e0e0", font=("Segoe UI", 10)
        )
        style.configure(
            "TButton",
            background="#3c3f41",
            foreground="#ffffff",
            font=("Segoe UI", 10),
            padding=6,
        )
        style.map("TButton", background=[("active", "#505354")])

        style.configure(
            "TEntry",
            fieldbackground="#3c3f41",
            foreground="#ffffff",
            insertcolor="#ffffff",
        )
        self.grid_on = False
        self.setup_data()
        self.setup_plot()
        self.create_controls()
        self.draw_plot()

    def setup_data(self):

        self.zs = np.linspace(0, 100000, num=1001)
        self.T = self.geotherm.interpolate(self.zs)

        self.xmin = -2.5
        self.xmax = 1.5
        self.ymin = 0
        self.ymax = 100

        self.lines = []
        self.labeld = {}
        self.lined = {}
        self.matlist = []

    def setup_plot(self):
        self.fig = plt.Figure(figsize=(8, 8))
        self.ax = self.fig.add_subplot(111)
        self.ax2 = self.ax.twiny()

        # Create the canvas early so it exists for all later calls
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().grid(
            row=0, column=1, sticky="nsew", padx=20, pady=10
        )

        self.ax.plot([0, 0], [self.ymin, self.ymax], lw=1, c="black", alpha=0.5)

        nmax = len(self.mat_dbase)
        for n, mat in enumerate(self.mat_dbase):
            self.matlist.append(mat)
            label = f"{mat.type}, {mat.source}"
            self.labeld[label] = n
            color = cm.romaO(n / nmax)
            sigma_plot, z_plot = re.compute_dsigma(
                mat, self.zs, self.T, self.strain_rate
            )
            (line,) = self.ax.plot(
                sigma_plot / 1e9,
                z_plot / 1000,
                label=label,
                c=color,
                zorder=10,
                linewidth=2,
            )
            self.lines.append(line)

        self.ax.set_xlabel("Differential stress in GPa")
        self.ax.set_ylabel("Depth in km")
        self.ax.set_anchor("C")  # Center the plot within the axes
        self.ax.set_aspect("auto")  # Allow plot to fill axes area
        self.fig.subplots_adjust(left=0.18, right=0.82)  # Add equal margins left/right

        all_x = []
        for line in self.lines:
            xdata = line.get_xdata()
            if len(xdata) > 0:
                all_x.extend(xdata)

        if all_x:
            min_x, max_x = min(all_x), max(all_x)
            padding = 0.1 * (max_x - min_x) if max_x != min_x else 0.1
            self.ax.set_xlim(min_x - padding, max_x + padding)

        self.ax.set_ylim(self.ymax, self.ymin)

        (self.temp_line,) = self.ax2.plot(
            self.T - 273.15,
            self.zs / 1000,
            linestyle="--",
            color="orange",
            label="Temperature",
            visible=False,  # Make temperature line invisible by default
        )
        self.ax2.set_xlabel("Temperature in °C")
        self.ax2.set_xlim(0, 1350)
        self.ax2.get_xaxis().set_visible(False)

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)

        # --- Custom scrollable legend panel (left) ---
        self.legend_panel = tk.Frame(self, bg="#2b2b2b")
        self.legend_panel.grid(row=0, column=0, sticky="nsw", padx=10, pady=10)
        self.legend_panel.grid_rowconfigure(0, weight=1)
        self.legend_panel.grid_columnconfigure(0, weight=1)

        legend_frame_container = tk.Frame(self.legend_panel, bg="#2b2b2b")
        legend_frame_container.grid(row=0, column=0, sticky="nsew")
        legend_frame_container.grid_rowconfigure(0, weight=1)
        legend_frame_container.grid_columnconfigure(0, weight=1)

        legend_canvas = tk.Canvas(
            legend_frame_container,
            bg="#2b2b2b",
            highlightthickness=0,
            borderwidth=0,
            height=300,
        )
        legend_canvas.grid(row=0, column=0, sticky="nsew")
        legend_scroll = ttk.Scrollbar(
            legend_frame_container, orient="vertical", command=legend_canvas.yview
        )
        legend_scroll.grid(row=0, column=1, sticky="ns")
        legend_canvas.configure(yscrollcommand=legend_scroll.set)

        legend_frame = tk.Frame(legend_canvas, bg="#2b2b2b")
        legend_canvas.create_window((0, 0), window=legend_frame, anchor="nw")

        def on_frame_configure(event):
            legend_canvas.configure(scrollregion=legend_canvas.bbox("all"))

        legend_frame.bind("<Configure>", on_frame_configure)

        # Add mouse wheel scrolling support
        def on_mousewheel(event):
            legend_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        # Bind mouse wheel to legend canvas and frame
        legend_canvas.bind("<MouseWheel>", on_mousewheel)
        legend_frame.bind("<MouseWheel>", on_mousewheel)

        # Add material labels with color swatches
        self.legend_labels = []
        for i, line in enumerate(self.lines):
            color = line.get_color()
            label = line.get_label()
            row = tk.Frame(legend_frame, bg="#2b2b2b")
            # Set pady=0 to eliminate vertical gap between entries
            row.pack(fill="x", pady=0, padx=2)

            def rgba_to_hex(c):
                if isinstance(c, str):
                    return c
                r, g, b = [int(255 * x) for x in c[:3]]
                return f"#{r:02x}{g:02x}{b:02x}"

            swatch = tk.Label(row, width=2, bg=rgba_to_hex(color))
            swatch.pack(side="left", padx=(0, 5))
            # Increase ipady for more internal vertical padding
            lbl = tk.Label(row, text=label, fg="#e0e0e0", bg="#2b2b2b", anchor="w")
            lbl.pack(side="left", fill="x", expand=True, ipady=8)

            def toggle_line(event, idx=i, label_widget=lbl):
                line = self.lines[idx]
                # Get the current toggle state from the label color, not the line visibility
                # (since line visibility might be temporarily changed by hover)
                label_color = label_widget.cget("fg")
                is_currently_visible = label_color == "#e0e0e0"

                # Toggle the state
                new_visible = not is_currently_visible
                line.set_visible(new_visible)
                label_widget.config(fg="#e0e0e0" if new_visible else "#888888")

                # Reset alpha and ensure proper state
                line.set_alpha(1.0)

                self.canvas.draw()
                show_material_info(self.matlist[idx] if new_visible else None)

            def on_legend_enter(event, idx=i):
                highlight_legend(idx)
                show_material_info(self.matlist[idx])
                # Always show the hovered line, regardless of its toggle state
                # Other lines: only show if they were originally visible (toggled on)
                for j, line in enumerate(self.lines):
                    if j == idx:
                        # Always make the hovered line visible at full opacity
                        line.set_visible(True)
                        line.set_alpha(1.0)
                    else:
                        # Check if this line is originally visible (toggled on)
                        label_color = self.legend_labels[j].cget("fg")
                        is_originally_visible = label_color == "#e0e0e0"
                        if is_originally_visible:
                            line.set_visible(True)
                            line.set_alpha(0.2)
                        else:
                            # Keep hidden lines hidden during hover
                            line.set_visible(False)
                self.canvas.draw()

            def on_legend_leave(event):
                unhighlight_legend()
                # Restore original visibility states based on legend label colors
                for j, line in enumerate(self.lines):
                    # Check if the line should be visible based on legend label color
                    label_color = self.legend_labels[j].cget("fg")
                    is_originally_visible = label_color == "#e0e0e0"  # visible color
                    line.set_visible(is_originally_visible)
                    line.set_alpha(1.0)
                self.canvas.draw()

            lbl.bind("<Button-1>", toggle_line)
            lbl.bind("<Enter>", on_legend_enter)
            lbl.bind("<Leave>", on_legend_leave)
            lbl.bind("<MouseWheel>", on_mousewheel)  # Add mouse wheel support to labels
            self.legend_labels.append(lbl)

        # Material info panel (fixed width, below legend)
        self.material_info_frame = tk.Frame(
            self.legend_panel,
            bg="#232323",
            bd=1,
            relief="groove",
            width=220,
            height=500,
        )
        self.material_info_frame.grid(
            row=1, column=0, sticky="ew", padx=5, pady=(10, 0)
        )
        self.material_info_frame.grid_propagate(False)
        # --- Refactored show_material_info: only update values, not labels ---
        self._material_info_labels = {}
        self._material_info_values = {}

        def show_material_info(mat):
            # Only create labels once
            if not self._material_info_labels:
                # Clear frame
                for widget in self.material_info_frame.winfo_children():
                    widget.destroy()
                row = 0
                # Title
                tk.Label(
                    self.material_info_frame,
                    text="Material Info",
                    fg="#e0e0e0",
                    bg="#232323",
                    font=("Segoe UI", 10, "bold"),
                ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(0, 4))
                row += 1
                # General attributes
                for key in ["id", "source", "type"]:
                    lbl = tk.Label(
                        self.material_info_frame,
                        text=f"{key.title()}:",
                        fg="#e0e0e0",
                        bg="#232323",
                        font=("Segoe UI", 9, "bold"),
                    )
                    lbl.grid(row=row, column=0, sticky="w")
                    val = tk.Label(
                        self.material_info_frame,
                        text="",
                        fg="#e0e0e0",
                        bg="#232323",
                        font=("Segoe UI", 11),
                    )
                    val.grid(row=row, column=1, sticky="w")
                    self._material_info_labels[key] = lbl
                    self._material_info_values[key] = val
                    row += 1
                # Byerlee
                tk.Label(
                    self.material_info_frame,
                    text="Byerlee:",
                    fg="#e0e0e0",
                    bg="#232323",
                    font=("Segoe UI", 9, "underline"),
                ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(4, 0))
                row += 1
                for key in ["fc_e", "fc_c", "lambda_pore", "rho_b"]:
                    lbl = tk.Label(
                        self.material_info_frame,
                        text=f"  {key}:",
                        fg="#e0e0e0",
                        bg="#232323",
                        font=("Segoe UI", 9),
                    )
                    lbl.grid(row=row, column=0, sticky="w")
                    val = tk.Label(
                        self.material_info_frame,
                        text="",
                        fg="#e0e0e0",
                        bg="#232323",
                        font=("Segoe UI", 11),
                    )
                    val.grid(row=row, column=1, sticky="w")
                    self._material_info_labels[key] = lbl
                    self._material_info_values[key] = val
                    row += 1
                # Dislocation
                tk.Label(
                    self.material_info_frame,
                    text="Dislocation:",
                    fg="#e0e0e0",
                    bg="#232323",
                    font=("Segoe UI", 9, "underline"),
                ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(4, 0))
                row += 1
                for key in ["a_disloc", "q_disloc", "n"]:
                    lbl = tk.Label(
                        self.material_info_frame,
                        text=f"  {key}:",
                        fg="#e0e0e0",
                        bg="#232323",
                        font=("Segoe UI", 9),
                    )
                    lbl.grid(row=row, column=0, sticky="w")
                    val = tk.Label(
                        self.material_info_frame,
                        text="",
                        fg="#e0e0e0",
                        bg="#232323",
                        font=("Segoe UI", 11),
                    )
                    val.grid(row=row, column=1, sticky="w")
                    self._material_info_labels[key] = lbl
                    self._material_info_values[key] = val
                    row += 1
                # Diffusion
                tk.Label(
                    self.material_info_frame,
                    text="Diffusion:",
                    fg="#e0e0e0",
                    bg="#232323",
                    font=("Segoe UI", 9, "underline"),
                ).grid(row=row, column=0, columnspan=2, sticky="w", pady=(4, 0))
                row += 1
                for key in ["a_diff", "q_diff", "d", "m"]:
                    lbl = tk.Label(
                        self.material_info_frame,
                        text=f"  {key}:",
                        fg="#e0e0e0",
                        bg="#232323",
                        font=("Segoe UI", 9),
                    )
                    lbl.grid(row=row, column=0, sticky="w")
                    val = tk.Label(
                        self.material_info_frame,
                        text="",
                        fg="#e0e0e0",
                        bg="#232323",
                        font=("Segoe UI", 11),
                    )
                    val.grid(row=row, column=1, sticky="w")
                    self._material_info_labels[key] = lbl
                    self._material_info_values[key] = val
                    row += 1
                # Other

            # Update values only
            def fmt(val):
                if isinstance(val, float):
                    if abs(val) > 1e4 or (abs(val) < 1e-2 and val != 0):
                        return f"{val:.2e}"
                    else:
                        return f"{val:.4g}"
                return str(val)

            if mat is None:
                for key, val_lbl in self._material_info_values.items():
                    val_lbl.config(text="")
                return
            for key in [
                "id",
                "source",
                "type",
                "fc_e",
                "fc_c",
                "lambda_pore",
                "rho_b",
                "a_disloc",
                "q_disloc",
                "n",
                "a_diff",
                "q_diff",
                "d",
                "m",
            ]:
                value = getattr(mat, key, "")
                self._material_info_values[key].config(text=fmt(value))

        show_material_info(None)

        # Track which legend label is highlighted
        self.highlighted_legend_idx = None

        def highlight_legend(idx):
            for j, lbl in enumerate(self.legend_labels):
                if j == idx:
                    lbl.config(bg="#444444")
                else:
                    lbl.config(bg="#2b2b2b")
            self.highlighted_legend_idx = idx

        def unhighlight_legend():
            for lbl in self.legend_labels:
                lbl.config(bg="#2b2b2b")
            self.highlighted_legend_idx = None

        self.canvas.draw()

    def update_legend_hitboxes(self, event):
        """Update rectangle positions after drawing"""
        for legtext in self.legend.get_texts():
            if hasattr(legtext, "_click_rect"):
                # Get text position in figure coordinates
                bbox = legtext.get_window_extent()
                fig_bbox = bbox.transformed(self.fig.transFigure.inverted())

                # Expand rectangle to full legend width
                leg_bbox = self.legend.get_window_extent()
                fig_leg_bbox = leg_bbox.transformed(self.fig.transFigure.inverted())

                # Update rectangle position
                rect = legtext._click_rect
                rect.set_bounds(
                    fig_leg_bbox.x0, fig_bbox.y0, fig_leg_bbox.width, fig_bbox.height
                )

    def on_legend_pick(self, event):
        # Handle rectangle clicks
        for rect, index in self.legend_rects:
            if event.artist == rect:
                line = self.lines[index]
                visible = not line.get_visible()
                line.set_visible(visible)

                # Update legend appearance
                legend_lines = self.legend.get_lines()
                legend_texts = self.legend.get_texts()
                alpha = 1.0 if visible else 0.2
                legend_lines[index].set_alpha(alpha)
                legend_texts[index].set_alpha(alpha)

                self.canvas.draw()
                break

    def export_figure(self):
        file_path = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[
                ("PNG Image", "*.png"),
                ("PDF File", "*.pdf"),
                ("All Files", "*.*"),
            ],
            title="Save Figure As",
        )
        if file_path:
            self.fig.savefig(
                file_path,
                dpi=300,
                # bbox_inches="tight",
                facecolor=self.fig.get_facecolor(),
            )
            print(f"Figure saved to {file_path}")

    def toggle_grid(self):
        self.grid_on = not self.grid_on
        self.ax.grid(self.grid_on)
        self.canvas.draw()
        self.grid_button.config(text="Grid Off" if self.grid_on else "Grid On")

    def toggle_temperature_line(self):
        if hasattr(self, "temp_line") and hasattr(self, "ax2"):
            visible = not self.temp_line.get_visible()
            self.temp_line.set_visible(visible)

            # Toggle the entire top axis visibility
            self.ax2.get_xaxis().set_visible(visible)
            self.ax2.get_yaxis().set_visible(
                False
            )  # Just in case, disable Y on twin axis

            self.canvas.draw()

    def create_controls(self):
        control_frame = ttk.Frame(self)
        control_frame.grid(row=0, column=2, sticky="nsw", padx=40, pady=10)
        control_frame.grid_rowconfigure(0, weight=1)
        control_frame.grid_columnconfigure(0, weight=1)

        ttk.Label(control_frame, text="Strain rate:").pack(pady=5)
        self.strain_entry = ttk.Entry(control_frame, width=15)
        self.strain_entry.insert(0, str(self.strain_rate))
        self.strain_entry.pack(pady=5)
        self.strain_entry.bind("<Return>", self.on_strain_rate_change)
        self.strain_entry.bind("<FocusOut>", self.on_strain_rate_change)

        self.all_visible = True  # Track current state for toggle
        self.toggle_all_button = ttk.Button(
            control_frame, text="Hide All", command=self.toggle_all_lines
        )
        self.toggle_all_button.pack(pady=5)

        self.grid_button = ttk.Button(
            control_frame, text="Grid Off", command=self.toggle_grid
        )

        self.grid_button.pack(pady=5)

        ttk.Button(
            control_frame, text="Toggle Temp Line", command=self.toggle_temperature_line
        ).pack(pady=5)

        ttk.Button(
            control_frame, text="Comp. & Ext.", command=lambda: self.set_xlim()
        ).pack(pady=5)

        ttk.Button(
            control_frame, text="Only Comp.", command=lambda: self.set_xlim(xmax=0)
        ).pack(pady=5)

        ttk.Button(
            control_frame, text="Only Ext.", command=lambda: self.set_xlim(xmin=0)
        ).pack(pady=5)

        ttk.Button(
            control_frame,
            text="Toggle Background",
            command=self.toggle_background_model,
        ).pack(pady=5)

        # Spacer frame to push the export button to the bottom
        ttk.Frame(control_frame).pack(expand=True, fill=tk.BOTH)

        ttk.Button(
            control_frame, text="Export Figure", command=self.export_figure
        ).pack(side=tk.BOTTOM, pady=10)

    def toggle_all_lines(self):
        if self.all_visible:
            for i, line in enumerate(self.lines):
                line.set_visible(False)
                self.legend_labels[i].config(fg="#888888")
            self.toggle_all_button.config(text="Show All")
        else:
            for i, line in enumerate(self.lines):
                line.set_visible(True)
                self.legend_labels[i].config(fg="#e0e0e0")
            self.toggle_all_button.config(text="Hide All")
        self.all_visible = not self.all_visible
        self.canvas.draw()

    def on_strain_rate_change(self, event=None):
        try:
            self.strain_rate = float(self.strain_entry.get())
            all_x = []

            # Update main lines
            for i, mat in enumerate(self.matlist):
                sigma_plot, z_plot = re.compute_dsigma(
                    mat, self.zs, self.T, self.strain_rate
                )
                self.lines[i].set_xdata(sigma_plot / 1e9)
                self.lines[i].set_ydata(z_plot / 1000)
                if self.lines[i].get_visible():
                    all_x.extend(sigma_plot / 1e9)

            # Update extra line if it exists
            if hasattr(self, "extra_line") and self.extra_line in self.ax.lines:
                extra_material = re.load_model()
                # Get the top boundary at the center for the extra model
                x, y = extra_material._xy_index_to_value()
                top_layer = extra_material.layers[0]
                top_z = top_layer.data[
                    (top_layer.data["x"] == x) & (top_layer.data["y"] == y)
                ]["z"].values[0]
                sigma_plot, z_plot = re.compute_dsigma(
                    extra_material, self.zs, self.T, self.strain_rate
                )
                z_plot_shifted = z_plot - top_z  # shift to actual model depth
                self.extra_line.set_xdata(sigma_plot / 1e9)
                self.extra_line.set_ydata(z_plot_shifted / 1000)
                all_x.extend(sigma_plot / 1e9)

            # Adjust x-limits
            if all_x:
                min_x, max_x = min(all_x), max(all_x)
                padding = 0.1 * (max_x - min_x) if max_x != min_x else 0.1
                self.ax.set_xlim(min_x - padding, max_x + padding)

            self.canvas.draw()

        except ValueError:
            print("Invalid strain rate input.")

    def set_grid(self, val):
        self.ax.grid(val)
        self.canvas.draw()

    def set_xlim(self, xmin=None, xmax=None):
        all_x = []

        # Collect x-data from visible main lines
        for line in self.lines:
            if line.get_visible():
                xdata = line.get_xdata()
                if len(xdata) > 0:
                    all_x.extend(xdata)

        # Include extra line if applicable
        if hasattr(self, "extra_line") and self.extra_line.get_visible():
            xdata = self.extra_line.get_xdata()
            if len(xdata) > 0:
                all_x.extend(xdata)

        if not all_x:
            return  # No data to determine limits

        # Compute new x limits
        data_min = min(all_x)
        data_max = max(all_x)
        actual_xmin = data_min if xmin is None else xmin
        actual_xmax = data_max if xmax is None else xmax

        if xmin is None:
            pad_min = 0.1 * (actual_xmax - actual_xmin)
            actual_xmin -= pad_min
        if xmax is None:
            pad_max = 0.1 * (actual_xmax - actual_xmin)
            actual_xmax += pad_max

        # Set limits
        self.ax.set_xlim(actual_xmin, actual_xmax)

        # Reposition background labels (to right edge)
        if hasattr(self, "bg_artists"):
            x_right = self.ax.get_xlim()[1]
            for artist in self.bg_artists:
                if isinstance(artist, plt.Text):
                    artist.set_x(x_right)

        self.canvas.draw()

    def toggle_background_model(self):
        if self.bg_model_visible:
            for artist in self.bg_artists:
                artist.remove()
            self.bg_artists.clear()
            if hasattr(self, "extra_line") and self.extra_line in self.ax.lines:
                self.extra_line.remove()
                del self.extra_line
            self.bg_model_visible = False
        else:
            # Only load and compute if not cached
            if self.model is None:
                try:
                    self.model = re.load_model()
                except Exception as e:
                    messagebox.showwarning(
                        "Model Load Error",
                        f"Failed to load default model:\n{e}\n\nPlease select a model.ini file.",
                    )
                    file_path = filedialog.askopenfilename(
                        title="Select model.ini file",
                        filetypes=[("INI files", "*.ini"), ("All files", "*.*")],
                    )
                    if file_path:
                        try:
                            self.model = re.load_model(file_path)
                        except Exception as e2:
                            messagebox.showerror(
                                "Model Load Error",
                                f"Failed to load selected model:\n{e2}",
                            )
                            return
                    else:
                        return  # User cancelled file dialog

                # Get the top boundary at the center
                x, y = self.model._xy_index_to_value()
                top_layer = self.model.layers[0]
                self.bg_top_z = top_layer.data[
                    (top_layer.data["x"] == x) & (top_layer.data["y"] == y)
                ]["z"].values[0]
                # Compute sigma and z arrays
                self.bg_sigma_plot, self.bg_z_plot = re.compute_dsigma(
                    self.model, self.zs, self.T, self.strain_rate
                )
                self.bg_z_plot_shifted = self.bg_z_plot - self.bg_top_z

            model = self.model
            sigma_plot = self.bg_sigma_plot
            z_plot_shifted = self.bg_z_plot_shifted

            layers = model.print_layers_at(tag=True)
            n_layers = len(layers)
            colormap = cm.lipari

            for i, (z_top, z_bottom, material) in enumerate(layers):
                color = colormap(i / max(1, n_layers - 1))
                top_km = -z_top / 1000
                bottom_km = -z_bottom / 1000
                span = self.ax.axhspan(
                    ymin=top_km, ymax=bottom_km, color=color, alpha=0.3, zorder=0
                )
                self.bg_artists.append(span)
                mid_km = (top_km + bottom_km) / 2
                label = self.ax.text(
                    self.ax.get_xlim()[1],
                    y=mid_km,
                    s=material.type if hasattr(material, "type") else str(material),
                    fontsize=8,
                    color="black",
                    ha="right",
                    va="center",
                    alpha=0.6,
                    zorder=1,
                )
                self.bg_artists.append(label)

            (self.extra_line,) = self.ax.plot(
                sigma_plot / 1e9,
                z_plot_shifted / 1000,  # Convert depth to km
                label=f"Custom Model",
                linestyle="--",
                linewidth=2,
                color="cornflowerblue",
                zorder=20,  # Ensure extra line is above all others
            )

            self.bg_model_visible = True

        self.canvas.draw()

    def draw_plot(self):
        # Placeholder for adding legend interactivity if needed
        self.canvas.draw()


if __name__ == "__main__":
    app = RheologyExplorerApp()
    app.mainloop()
