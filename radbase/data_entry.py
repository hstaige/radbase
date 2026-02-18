"""
A simple UI to assist with the compilation of previous nuclear radius information.
"""

import json
import re
import tkinter as tk
from collections import OrderedDict
from dataclasses import dataclass
from importlib.resources import files
from importlib.resources.abc import Traversable
from pathlib import Path
from tkinter import messagebox, scrolledtext, ttk
from typing import Any, Callable, Literal, NamedTuple, Optional, Protocol

import uncertainties

config = {
    'compilation_dir': files('radbase') / 'compilation',
    'reference_path': files('radbase') / 'compilation/references.json'
}


@dataclass
class Reference:
    title: str
    author: list[str]
    year: str
    month: Optional[str] = None
    journal: Optional[str] = None
    pages: Optional[int] = None
    volume: Optional[int] = None
    url: Optional[str] = None
    doi: Optional[str] = None


class ScrollableFrame(ttk.Frame):
    def __init__(self, container, *args, **kwargs):
        super().__init__(container, *args, **kwargs)

        canvas = tk.Canvas(self)
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)

        self.scrollable_frame = ttk.Frame(canvas)

        # Create window and store ID
        self.window_id = canvas.create_window(
            (0, 0),
            window=self.scrollable_frame,
            anchor="nw"
        )

        # Update scrollregion when content changes
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        # KEY FIX: make inner frame track canvas width
        canvas.bind(
            "<Configure>",
            lambda e: canvas.itemconfig(self.window_id, width=e.width)
        )

        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")


class ToolTip:
    def __init__(self, widget, text: str, delay: int = 500):
        self.widget = widget
        self.text = text
        self.delay = delay
        self._after_id = None
        self._tip = None

        widget.bind("<Enter>", self._schedule)
        widget.bind("<Leave>", self._hide)
        widget.bind("<ButtonPress>", self._hide)

    def _schedule(self, _):
        self._after_id = self.widget.after(self.delay, self._show)

    def _show(self):
        if self._tip or not self.text:
            return

        x = self.widget.winfo_rootx() + 20
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 5

        self._tip = tk.Toplevel(self.widget)
        self._tip.wm_overrideredirect(True)
        self._tip.wm_geometry(f"+{x}+{y}")

        label = ttk.Label(
            self._tip,
            text=self.text,
            background="#ffffe0",
            relief="solid",
            borderwidth=1,
            padding=6,
            wraplength=300
        )
        label.pack()

    def _hide(self, _=None):
        if self._after_id:
            self.widget.after_cancel(self._after_id)
            self._after_id = None
        if self._tip:
            self._tip.destroy()
            self._tip = None


class FieldProcessor(Protocol):
    def process(self, widget: tk.Widget) -> dict[str, Any]:
        pass


class CastProcessor:

    def __init__(self, field_type: type, key: str, allows_empty: bool = False):
        self.key = key
        self.field_type = field_type
        self.allows_empty = allows_empty

    def process(self, widget: tk.Widget):
        raw = widget.get()
        if raw == "" and not self.allows_empty:
            raise ValueError("Field is empty")
        return {self.key: self.field_type(raw)}


class GroupedProcessor:

    def __init__(self, processors: list[FieldProcessor]):
        self.processors = processors

    def process(self, widget: tk.Widget, allows_empty=False):
        result = {}

        for sub_widget, processor in zip(widget.widgets, self.processors):
            try:
                sub_result = processor.process(sub_widget)
            except ValueError:
                continue

            result = result | sub_result

        if result == {} and not allows_empty:
            raise ValueError(f'At least one sub_widget of {widget} must be non-empty for GroupedProcessor')

        return result


class XORProcessor:
    # Requires exactly one of the two fields to be filled, and then returns data from filled out field.
    def __init__(self, processors: list[FieldProcessor]):
        self.processors = processors

    def process(self, widget: tk.Widget):
        result = {}

        if len(widget.widgets) != 2:
            raise ValueError('XOR processor applied to widget that does not have two fields.')

        num_filled = 0
        for sub_widget, processor in zip(widget.widgets, self.processors):
            try:
                sub_result = processor.process(sub_widget)
            except ValueError:
                continue

            num_filled += 1
            result = result | sub_result

        if num_filled != 1:
            raise ValueError(f'Exactly one of the XOR fields must be filled, {num_filled} were filled.')

        return result


class NumberWithUncertaintyProcessor:

    def __init__(self, key: str):
        self.key = key

    def process(self, widget: tk.Widget):
        raw = widget.get().strip()
        try:
            value = float(raw)
            return {self.key: {'value': value, 'uncertainty': None}}
        except ValueError:
            try:
                uvar = uncertainties.ufloat_fromstr(raw)
                return {self.key: {'value': uvar.n, 'uncertainty': uvar.s}}
            except ValueError:
                raise ValueError(
                    f'Entered value of {raw} does not match the number with uncertainty pattern. Ex. 0.65(2) or 0.32')


class NuclearPolarizationProcessor:

    @staticmethod
    def process(widget: tk.Widget):
        selection = widget.selection
        fit_widget = widget.fit_widget
        prev_calced_widget = widget.prev_calced_widget

        calced_processor = PreviousDataProcessor(key='Calculated Nuclear Polarization')
        varied_processor = VariableNumberProcessor(
            GroupedProcessor([
                XORProcessor(
                    [TransitionProcessor(), CastProcessor(str, key='Level')]
                ),
                NumberWithUncertaintyProcessor('Energy [keV]')
            ]))

        mode = selection.get()
        result = {'Nuclear Polarization Method': mode}

        if mode == 'Calculated':
            return result | calced_processor.process(prev_calced_widget)
        elif mode == 'Vary':
            return result | varied_processor.process(fit_widget)
        elif mode == 'Mixed':
            return result | calced_processor.process(prev_calced_widget) | varied_processor.process(fit_widget)
        else:
            raise ValueError(f'{mode} is not in [Calculated, Vary, Mixed]')


class VariableNumberProcessor:

    def __init__(self, processor: FieldProcessor, suffixes: list[str] | None = None):
        self.processor = processor
        self.suffixes = suffixes if suffixes else [f'_{chr(ord('A') + i)}' for i in range(26)]

    def process(self, widget: tk.Widget) -> dict:
        data = {}
        for suffix, (label, widget) in zip(self.suffixes, widget.entries):
            sub_result = self.processor.process(widget)
            sub_result = {key + suffix: value for key, value in sub_result.items()}
            data = data | sub_result
        return data


class MultiTransitionProcessor:

    @staticmethod
    def process(widget: tk.Widget):
        data = {}
        for label, container in widget.transition_entries:
            data[label.cget("text")] = TransitionProcessor().process(container)
        return data


class NuclideProcessor:

    def __init__(self, key: str = 'Nuclide'):
        self.key = key

    def process(self, widget: tk.Widget):
        raw = widget.get()
        if re.match(r'[a-zA-Z]{1,2}\d{1,3}', raw) or re.match(r'[a-zA-Z]{1,2}nat', raw):
            return {self.key: raw}
        else:
            raise ValueError(
                'Entered Nuclide does not match pattern of two letters and one to three numbers (Ex. Pb208')


class PreviousDataProcessor:
    def __init__(self, key: str = 'Relies On'):
        self.key = key

    def process(self, widget: tk.Widget) -> dict[str, list[str]]:
        selected = [
            key
            for key, var in widget.selection_vars.items()
            if var.get()
        ]
        return {self.key: selected}


class ReferenceProcessor:
    def __init__(self, reference_path: Traversable, key: str = 'Reference'):
        self.key = key
        self.reference_path = reference_path

    @property
    def references(self):
        pass

    @references.getter
    def references(self):
        with open(self.reference_path, 'r') as f:
            return json.load(f)

    def process(self, widget: tk.Widget):
        value = widget.get().strip()

        if value not in sorted(self.references.keys()):
            raise ValueError(
                f"Reference '{value}' is not a known citation key"
            )

        return {self.key: value}


class SelectProcessor:
    def __init__(self, options: list[str], key: str):
        self.key = key
        self.options = set(options)

    def process(self, widget: tk.Widget) -> dict[str, str]:
        value = widget.get().strip()

        if not value:
            raise ValueError("No option selected")

        if value not in self.options:
            raise ValueError(f"'{value}' is not a valid selection")

        return {self.key: value}


class TransitionProcessor:
    SIMPLE_LEVEL_PATTERN = re.compile(
        r"""
        ^
        \d+                # principal quantum number
        [spdf]             # orbital
        (\d+/\d+)?         # optional j (e.g. 3/2)
        $
        """,
        re.VERBOSE
    )

    HYBRID_LEVEL_PATTERN = re.compile(
        r"""
        ^
        (\d+[+-],)?          # nuclear spin/parity (e.g. 2+, 0-)
        \d+                # principal quantum number
        [spdf]             # orbital
        \d/\d            # required j
        $
        """,
        re.VERBOSE
    )

    def __init__(self, key: str = 'Transition'):
        self.key = key

    def process(self, widget: tk.Widget) -> dict[str, dict[str, str]]:
        upper = widget.upper_entry.get().strip()
        lower = widget.lower_entry.get().strip()

        if not upper or not lower:
            raise ValueError("Both upper and lower transition levels must be provided")

        upper_simple = self.SIMPLE_LEVEL_PATTERN.match(upper)
        lower_simple = self.SIMPLE_LEVEL_PATTERN.match(lower)

        upper_hybrid = self.HYBRID_LEVEL_PATTERN.match(upper)
        lower_hybrid = self.HYBRID_LEVEL_PATTERN.match(lower)

        # both simple or both hybrid
        if upper_simple and lower_simple or upper_hybrid and lower_hybrid:
            return {self.key: {'Upper': upper, 'Lower': lower}}

        raise ValueError(
            "Invalid transition format.\n"
            "Both levels must be either:\n"
            "  - Simple atomic levels (e.g. 2p3/2)\n"
            "  - Hybrid nuclear+atomic levels (e.g. 2+,1s1/2)"
        )


class WidgetCreator(Protocol):
    def create_widget(self, frame, row) -> tk.Widget:
        pass


class DefaultWidgetCreator:
    @staticmethod
    def create_widget(frame, row) -> tk.Widget:
        entry = ttk.Entry(frame, width=40)
        entry.grid(row=row, column=1, pady=4, sticky='w')
        return entry


class GroupedWidgetCreator:

    def __init__(self, labels: list[str], widget_creators: list[WidgetCreator]):
        self.labels = labels
        self.widget_creators = widget_creators
        self.widgets = []

    def create_widget(self, frame, row) -> tk.Widget:
        grouped_frame = tk.Frame(frame)
        grouped_frame.grid(row=row, column=1, pady=4, sticky='w')
        grouped_frame.widgets = []

        for i, (label, widget_creator) in enumerate(zip(self.labels, self.widget_creators)):
            tk.Label(grouped_frame, text=label).grid(row=i, column=0, pady=4, sticky='w')
            grouped_frame.widgets.append(widget_creator.create_widget(grouped_frame, i))

        return grouped_frame


class NuclearPolarizationWidgetCreator:

    def __init__(self, compilation_path: Traversable):
        self.compilation_path = compilation_path

    def create_widget(self, frame, row) -> tk.Widget:
        container = ttk.Frame(frame)
        container.grid(row=row, column=1, sticky="w")

        # Selection combobox (Calculated / Vary / Mixed)
        selection = SelectWidgetCreator(
            ["Calculated", "Vary", "Mixed"]
        ).create_widget(container, row=0)
        selection.grid(row=0, column=0, sticky="w", pady=4)

        # Sub-widgets
        prev_calced_widget = PreviousDataWidgetCreator(
            self.compilation_path,
            filter_regex="muonic_nuclear"
        ).create_widget(container, row=1)

        pair = GroupedWidgetCreator(['', 'Energy [keV]'], [GroupedWidgetCreator(['Transition', 'Level'],
                                                                                [TransitionWidgetCreator(),
                                                                                 DefaultWidgetCreator()]),
                                                           DefaultWidgetCreator()])
        fit_widget = VariableNumberWidgetCreator(pair).create_widget(container, row=1)

        # Start hidden
        prev_calced_widget.grid_remove()
        fit_widget.grid_remove()

        def handle_selection(event=None):
            mode = selection.get()

            if mode == "Calculated":
                prev_calced_widget.grid(row=1, column=0)
                fit_widget.grid_remove()

            elif mode == "Vary":
                prev_calced_widget.grid_remove()
                fit_widget.grid(row=1, column=0)

            elif mode == "Mixed":
                prev_calced_widget.grid(row=1, column=0)
                fit_widget.grid(row=2, column=0)

        selection.bind("<<ComboboxSelected>>", handle_selection)
        handle_selection()

        container.selection = selection
        container.fit_widget = fit_widget
        container.prev_calced_widget = prev_calced_widget
        return container


class VariableNumberWidgetCreator:

    def __init__(self, widget_creator: WidgetCreator):
        self.widget_creator = widget_creator

    def create_widget(self, frame, row):

        class MultiContainer(ttk.Frame):  # Mainly for type hinting
            def __init__(self, entries, buttons, *args, **kwargs):
                super().__init__(frame, *args, **kwargs)
                self.entries: list[tuple[tk.Label, tk.Widget]] = entries
                self.buttons: list[tk.Button] = buttons

            def regrid(self):
                for i, (label, entry) in enumerate(self.entries):
                    label.grid(row=i), entry.grid(row=i)
                for button in self.buttons:
                    button.grid(row=len(self.entries))

        container = MultiContainer(entries=[], buttons=[])
        container.grid(row=row, column=1, pady=4, sticky="w")

        def add_new_entry(*_):
            number_of_rows = len(container.entries)
            new_letter = f'{chr(ord("A") + number_of_rows)}'

            label = tk.Label(container, text=f'{new_letter}')
            label.grid(row=number_of_rows, column=0)

            next_entry = self.widget_creator.create_widget(frame=container, row=number_of_rows)

            container.entries.append((label, next_entry))
            container.regrid()

        def delete_entry(*_):
            label, entry = container.entries.pop()
            label.grid_forget(), entry.grid_forget()
            container.regrid()

        add_button = ttk.Button(container, text="Add", command=add_new_entry, takefocus=0)
        add_button.grid(row=0, column=0, pady=10)

        delete_button = ttk.Button(container, text="Remove", command=delete_entry, takefocus=0)
        delete_button.grid(row=0, column=1, pady=10)
        container.buttons = [add_button, delete_button]
        add_new_entry(), add_new_entry()

        return container


class PreviousDataWidgetCreator:

    def __init__(self, compilation_path: Traversable, filter_regex: str = ""):
        self.filter_regex = filter_regex
        self.previous_measurements = self.load_previous_measurements_from_json(compilation_path)

    @staticmethod
    def load_previous_measurements_from_json(compilation_path: Traversable) -> list[str]:
        previous_measurements = []
        for reference_dir in compilation_path.iterdir():
            if not reference_dir.is_dir():
                continue
            for json_file in reference_dir.iterdir():
                if '.json' not in str(json_file):
                    continue
                with open(json_file, "r", encoding="utf-8") as f:
                    previous_measurements.extend(json.load(f).keys())

        return sorted(previous_measurements)

    def create_widget(self, frame, row) -> tk.Widget:
        container = ttk.Frame(frame)
        container.grid(row=row, column=1, pady=4, sticky="w")

        # --- filter entry ---
        filter_var = tk.StringVar()
        filter_var.set(self.filter_regex)
        filter_entry = ttk.Entry(container, width=80, textvariable=filter_var)
        filter_entry.grid(row=0, column=0, columnspan=2, sticky="w")

        # --- scrollable checklist ---
        canvas = tk.Canvas(container, height=150)
        scrollbar = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
        checklist_frame = ttk.Frame(canvas)

        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.grid(row=1, column=0, sticky="nsew")
        scrollbar.grid(row=1, column=1, sticky="ns")

        canvas.create_window((0, 0), window=checklist_frame, anchor="nw")

        # --- checkbox state ---
        vars_by_key: dict[str, tk.BooleanVar] = {}
        checkbox_widgets: dict[str, ttk.Checkbutton] = {}

        for data_key in self.previous_measurements:
            var = tk.BooleanVar(value=False)
            check = ttk.Checkbutton(
                checklist_frame,
                text=data_key,
                variable=var
            )
            check.pack(expand=True, fill='x')
            vars_by_key[data_key] = var
            checkbox_widgets[data_key] = check

        def update_scrollregion(*_):
            canvas.configure(scrollregion=canvas.bbox("all"))

        checklist_frame.bind("<Configure>", update_scrollregion)

        # --- regex filtering logic ---
        def filter_list(*_):
            pattern = filter_var.get()
            try:
                regex = re.compile(pattern, re.IGNORECASE)
            except re.error:
                # Invalid regex: hide everything
                for chk in checkbox_widgets.values():
                    chk.pack_forget()
                update_scrollregion()
                return

            for key, chk in checkbox_widgets.items():
                if regex.search(key):
                    chk.pack(anchor="w")
                else:
                    chk.pack_forget()

        filter_list()
        update_scrollregion()

        filter_var.trace_add("write", filter_list)

        # --- select / deselect buttons ---
        button_frame = ttk.Frame(container)
        button_frame.grid(row=2, column=0, columnspan=2, sticky="w", pady=(4, 0))

        def select_visible(state: bool):
            for key, chk in checkbox_widgets.items():
                if chk.winfo_manager():  # only visible (i.e., packed) widgets
                    vars_by_key[key].set(state)

        ttk.Button(
            button_frame,
            text="Select all",
            command=lambda: select_visible(True),
            width=12,
        ).pack(side="left", padx=(0, 4))

        ttk.Button(
            button_frame,
            text="Deselect all",
            command=lambda: select_visible(False),
            width=12,
        ).pack(side="left")

        container.selection_vars = vars_by_key
        return container


class ReferenceWidgetcreator:
    _reference_options: list[str]

    def __init__(self, references: dict | None = None, reference_path: Traversable | None = None):
        if references is None and reference_path is None:
            return
        if references is not None and reference_path is not None:
            raise ValueError('Both references and reference_file supplied to ReferenceWidgetcreator.')
        elif references is not None:
            self.references = references
            self.reference_path = None
        else:
            self.load_references_from_json(reference_path)

    # ---------- file helpers ----------
    def load_references_from_json(self, reference_path: Traversable):
        self.reference_path = reference_path
        with open(reference_path, "r", encoding="utf-8") as f:
            self.references = json.load(f)

    def _save_references(self, data: dict):
        with open(self.reference_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=4)

    @property
    def reference_options(self):
        return self._reference_options

    @reference_options.getter
    def reference_options(self):
        return sorted(self.references)

    # Main widget!
    def create_widget(self, frame, row) -> tk.Widget:
        container = ttk.Frame(frame)
        container.grid(row=row, column=1, pady=4, sticky="w")

        combo = ttk.Combobox(
            container,
            width=40,
            state="normal",
            values=self.reference_options
        )
        combo.grid(row=0, column=0, pady=4)

        def filter_references(_):
            typed = combo.get().lower()
            filtered = [
                key for key in self.reference_options
                if typed in key.lower()
            ]
            combo["values"] = filtered

        combo.bind("<KeyRelease>", filter_references)

        def auto_complete(_):
            typed = combo.get().lower()
            filtered = [name for name in self.reference_options if typed in name.lower()]
            if filtered:
                combo.set(filtered[0])

        combo.bind("<Tab>", auto_complete)

        add_new = ttk.Button(
            container,
            text="New Ref",
            command=lambda: self._open_new_reference_form(combo)
        )
        add_new.grid(row=0, column=1, padx=(6, 0))

        return combo

    def _open_new_reference_form(self, combo: ttk.Combobox):
        win = tk.Toplevel()
        win.title("Add New Reference")
        win.transient()
        win.grab_set()

        fields = [
            ("Title", "title"),
            ("Authors (one per line)", "author"),
            ("Year", "year"),
            ("Month", "month"),
            ("Journal", "journal"),
            ("Volume", "volume"),
            ("Pages", "pages"),
            ("URL", "url"),
            ("DOI", "doi"),
        ]

        entries = OrderedDict()

        long_entries = {'Title', 'Authors (one per line)'}

        for row, (label, key) in enumerate(fields):
            ttk.Label(win, text=label).grid(row=row, column=0, sticky="w", pady=4, padx=4)
            if label not in long_entries:
                entry = tk.Entry(win, width=50)
            else:
                entry = scrolledtext.ScrolledText(win, width=50, height=3)

            def focus_next_widget(event):
                event.widget.tk_focusNext().focus()
                return "break"

            entry.bind("<Tab>", focus_next_widget)

            entry.grid(row=row, column=1, pady=4, padx=4)
            entries[key] = entry

        # ---------- submit logic ----------
        def submit():
            try:
                authors = [
                    a.strip()
                    for a in entries["author"].get("1.0", 'end-1c').splitlines()
                    if a.strip()
                ]
                if not authors:
                    raise ValueError("At least one author is required")

                year = entries["year"].get().strip()
                if not year:
                    raise ValueError("Year is required")

                title = entries['title'].get("1.0", 'end-1c')
                citation_key = authors[0].split()[-1].lower() + '_' + title.lower().split()[0] + '_' + year

                references = self.references
                if citation_key in references:
                    raise ValueError(f"Reference '{citation_key}' already exists")

                references[citation_key] = {
                    "title": title,
                    "author": authors,
                    "year": year,
                    "month": entries["month"].get() or None,
                    "journal": entries["journal"].get() or None,
                    "volume": entries["volume"].get(),
                    "pages": entries["pages"].get(),
                    "url": entries["url"].get() or None,
                    "doi": entries["doi"].get() or None,
                }

                self._save_references(references)

                # update combobox
                combo["values"] = self.reference_options
                combo.set(citation_key)

                win.destroy()

            except ValueError as e:
                messagebox.showerror("Invalid Reference", str(e), parent=win)

        ttk.Button(win, text="Submit", command=submit).grid(
            row=len(fields), column=0, columnspan=2, pady=10
        )

    @staticmethod
    def _parse_int(value: str):
        value = value.strip()
        return int(value) if value else None


class SelectWidgetCreator:
    def __init__(self, options: list[str]):
        if not options:
            raise ValueError("SelectWidgetCreator requires a non-empty options list")
        self.options = options

    def create_widget(self, frame, row) -> tk.Widget:
        combo = ttk.Combobox(
            frame,
            values=self.options,
            state="readonly",
            width=40
        )
        combo.grid(row=row, column=1, pady=4, sticky="w")

        # auto-select first option
        combo.current(0)

        return combo


class TransitionWidgetCreator:

    @staticmethod
    def create_widget(frame, row) -> tk.Widget:
        container = ttk.Frame(frame)
        container.grid(row=row, column=1, pady=4, sticky="w")

        upper_entry = ttk.Entry(container, width=18)
        lower_entry = ttk.Entry(container, width=18)

        upper_entry.grid(row=0, column=0, padx=(0, 5))
        ttk.Label(container, text="→").grid(row=0, column=1)
        lower_entry.grid(row=0, column=2, padx=(5, 0))

        # Attach for processor access
        container.upper_entry = upper_entry
        container.lower_entry = lower_entry

        return container


class MultiTransitionWidgetCreator:

    @staticmethod
    def create_widget(frame, row):

        class MultiContainer(ttk.Frame):  # Mainly for type hinting
            def __init__(self, transition_entries, buttons, *args, **kwargs):
                super().__init__(frame, *args, **kwargs)
                self.transition_entries: list[tuple[tk.Label, tk.Frame]] = transition_entries
                self.buttons: list[tk.Button] = buttons

            def regrid(self):
                for i, (label, transition_frame) in enumerate(self.transition_entries):
                    label.grid(row=i), transition_frame.grid(row=i)
                for button in self.buttons:
                    button.grid(row=len(self.transition_entries))

        container = MultiContainer(transition_entries=[], buttons=[])
        container.grid(row=row, column=1, pady=4, sticky="w")

        def add_new_transition(*_):
            number_of_rows = len(container.transition_entries)
            new_letter = f'{chr(ord("A") + number_of_rows)}'

            label = tk.Label(container, text=f'{new_letter}')
            label.grid(row=number_of_rows, column=0)

            transition_frame = ttk.Frame(container)  # Contains label and transition_frame
            transition_frame.grid(row=number_of_rows, column=1)

            upper_entry = ttk.Entry(transition_frame, width=18)
            lower_entry = ttk.Entry(transition_frame, width=18)

            upper_entry.grid(row=0, column=0, padx=(0, 5))
            ttk.Label(transition_frame, text="→").grid(row=0, column=1)
            lower_entry.grid(row=0, column=2, padx=(5, 0))

            transition_frame.upper_entry = upper_entry
            transition_frame.lower_entry = lower_entry

            container.transition_entries.append((label, transition_frame))
            container.regrid()

        def delete_transition(*_):
            label, transition_frame = container.transition_entries.pop()
            label.grid_forget(), transition_frame.grid_forget()
            container.regrid()

        add_button = ttk.Button(container, text="Add", command=add_new_transition, takefocus=0)
        add_button.grid(row=0, column=0, pady=10)

        delete_button = ttk.Button(container, text="Remove", command=delete_transition, takefocus=0)
        delete_button.grid(row=0, column=1, pady=10)
        container.buttons = [add_button, delete_button]
        add_new_transition(), add_new_transition()

        return container


class FieldSpec(NamedTuple):
    label: str
    processor: FieldProcessor
    widget_creator: WidgetCreator = DefaultWidgetCreator()
    hovertext: str | None = None


class InputTemplate(NamedTuple):
    name: str
    fields: list[FieldSpec]
    data_key: Callable[[dict[str, Any]], str]


# Useful fields that get reused
reference_field = FieldSpec("Reference", ReferenceProcessor(reference_path=config['reference_path']),
                            ReferenceWidgetcreator(reference_path=config['reference_path']))

nuclide_field = FieldSpec("Nuclide", NuclideProcessor(), hovertext='Provide a nuclide in the format Pb208')

transition_field = FieldSpec("Transition", TransitionProcessor(), TransitionWidgetCreator(),
                             hovertext='Enter the upper and lower levels of the atomic/muonic transition')

transition_or_level_field = FieldSpec("Transition\nor Level",
                                      XORProcessor([TransitionProcessor(), CastProcessor(str, key='Level')]),
                                      GroupedWidgetCreator(['Transition', 'Level'],
                                                           [TransitionWidgetCreator(), DefaultWidgetCreator()]),
                                      hovertext='The transition/level for which the nuclear polarization calculation was done.')

notes_field = FieldSpec("Notes", CastProcessor(str, key='Notes', allows_empty=True))

# LIST TEMPLATES HERE. Each template must include a FieldSpec named "Reference".
muonic_transition_energy_template = InputTemplate(
    name="Muonic Transition Energy",
    fields=[
        reference_field,
        nuclide_field,
        transition_field,
        FieldSpec("Energy [keV]", NumberWithUncertaintyProcessor("Energy [keV]")),
        notes_field
    ],
    data_key=lambda values: '_'.join([values['Reference'], 'muonic', values['Nuclide'],
                                      values['Transition']['Upper'], values['Transition']['Lower']]))

muonic_transition_energy_difference_template = InputTemplate(
    name="Muonic Transition Energy Difference (different nuclides)",
    fields=[
        reference_field,
        FieldSpec("Nuclide A", NuclideProcessor(key="Nuclide_A")),
        FieldSpec("Nuclide B", NuclideProcessor(key="Nuclide_B")),
        transition_or_level_field,
        FieldSpec("Energy Difference [keV] (A-B)", NumberWithUncertaintyProcessor("Energy Difference [keV] (A-B)")),
        notes_field
    ],
    data_key=lambda values: '_'.join(
        [values['Reference'], 'muonic_difference', values['Nuclide_A'], values['Nuclide_B'],
         *([values['Transition']['Upper'], values['Transition']['Lower']] if 'Transition' in values else [
             values['Level'], ])])
)

muonic_transition_energy_difference_diff_transition_template = InputTemplate(
    name="Muonic Transition Energy Difference (different transitions)",
    fields=[
        reference_field,
        nuclide_field,
        FieldSpec("Transitions",
                  VariableNumberProcessor(XORProcessor([TransitionProcessor(), CastProcessor(str, key='Level')]), ),
                  VariableNumberWidgetCreator(GroupedWidgetCreator(['Transition', 'Level'],
                                                                   [TransitionWidgetCreator(),
                                                                    DefaultWidgetCreator()])),
                  hovertext='Enter the upper and lower levels of the atomic/muonic transition'),
        FieldSpec("Energy Difference [keV] (A-B)", NumberWithUncertaintyProcessor("Energy Difference [keV] (A-B)")),
        notes_field
    ],
    data_key=lambda values: '_'.join(
        [values['Reference'], 'muonic_difference', values['Nuclide'],
         *list(
             sum([(field_name, field['Upper'], field['Lower']) for field_name, field in values.items() if
                  'Transition' in field_name.lower()],
                 ()))])
)


def mu_nuc_pol_key(values):
    if 'Transition' in values:
        return '_'.join([values['Reference'], 'muonic_nuclear_polarization', values['Nuclide'],
                         values['Transition']['Upper'], values['Transition']['Lower']])
    else:
        return '_'.join([values['Reference'], 'muonic_nuclear_polarization', values['Nuclide'],
                         values['Level']])


muonic_nuclear_polarization_calculation_template = InputTemplate(
    name="Calculated Muonic Nuclear Polarization",
    fields=[
        reference_field,
        nuclide_field,
        transition_or_level_field,
        FieldSpec("Energy [keV]", NumberWithUncertaintyProcessor("Energy [keV]")),
        notes_field
    ],
    data_key=mu_nuc_pol_key)


def mu_nuc_pol_key(values):
    if 'Transition' in values:
        return '_'.join([values['Reference'], 'muonic_qed', values['Nuclide'],
                         values['Transition']['Upper'], values['Transition']['Lower']])
    else:
        return '_'.join([values['Reference'], 'muonic_qed', values['Nuclide'],
                         values['Level']])


muonic_qed_calculation_template = InputTemplate(
    name="Calculated QED Correction",
    fields=[
        reference_field,
        nuclide_field,
        transition_or_level_field,
        FieldSpec("Energy [keV]", NumberWithUncertaintyProcessor("Energy [keV]")),
        notes_field
    ],
    data_key=mu_nuc_pol_key)

muonic_barret_theory_template = InputTemplate(
    name="Muonic Barrett Moment",
    fields=[
        reference_field,
        FieldSpec("Data used as input", PreviousDataProcessor(key='Previous Muonic Measurements'),
                  PreviousDataWidgetCreator(compilation_path=config['compilation_dir'],
                                            filter_regex='muonic.*_')),
        nuclide_field,
        FieldSpec('Rka [fm]', NumberWithUncertaintyProcessor('Rka [fm]')),
        FieldSpec('k [-]', CastProcessor(str, key='k [-]')),
        FieldSpec('alpha [1/fm]', CastProcessor(str, key='alpha [1/fm]')),
        FieldSpec('Cz [fm/keV]', CastProcessor(str, key='Cz [fm/keV]')),
        FieldSpec('Nuclear Polarization Method',
                  NuclearPolarizationProcessor(),
                  NuclearPolarizationWidgetCreator(compilation_path=config['compilation_dir']),
                  hovertext='Were the NP corrections from theory ("Calculated") or varied as part of the optimization ("Fit")?'),
        notes_field
    ],
    data_key=lambda values: '_'.join([values['Reference'], 'barrett_moment', values['Nuclide'],
                                      *['-'.join(s.split('_')[-2:]) for s in
                                        values['Previous Muonic Measurements']]])
)

muonic_barret_shift_template = InputTemplate(
    name="Muonic Barrett Moment Difference",
    fields=[
        reference_field,
        FieldSpec("Data used as input", PreviousDataProcessor(),
                  PreviousDataWidgetCreator(compilation_path=config['compilation_dir'],
                                            filter_regex='muonic.*_')),
        FieldSpec("Nuclide A", NuclideProcessor(key="Nuclide_A")),
        FieldSpec("Nuclide B", NuclideProcessor(key="Nuclide_B")),
        FieldSpec('Rka (A-B) [fm]', NumberWithUncertaintyProcessor('Rka [fm]')),
        FieldSpec('k [-]', CastProcessor(str, key='k [-]')),
        FieldSpec('alpha [1/fm]', CastProcessor(str, key='alpha [1/fm]')),
        FieldSpec('Cz [fm/keV]', CastProcessor(str, key='Cz [fm/keV]')),
        FieldSpec('Nuclear Polarization Method',
                  NuclearPolarizationProcessor(),
                  NuclearPolarizationWidgetCreator(compilation_path=config['compilation_dir']),
                  hovertext='Were the NP corrections from theory ("Calculated") or varied as part of the optimization ("Fit")?'),
        FieldSpec('QED Calculations Used',
                  PreviousDataProcessor(key='Calculated QED'),
                  PreviousDataWidgetCreator(compilation_path=config['compilation_dir'])),
        notes_field
    ],
    data_key=lambda values: '_'.join(
        [values['Reference'], 'barrett_moment_difference', values['Nuclide_A'], values['Nuclide_B'],
         *['-'.join(s.split('_')[-2:]) for s in
           values['Previous Muonic Measurements']]])
)

muonic_radius_template = InputTemplate(
    name="Muonic Radius",
    fields=[
        reference_field,
        FieldSpec("Data used as input", PreviousDataProcessor(),
                  PreviousDataWidgetCreator(compilation_path=config['compilation_dir'],
                                            filter_regex='muonic.*_')),
        nuclide_field,
        FieldSpec('R [fm]', NumberWithUncertaintyProcessor(key='R [fm]')),
        FieldSpec('Nuclear Polarization Method',
                  NuclearPolarizationProcessor(),
                  NuclearPolarizationWidgetCreator(compilation_path=config['compilation_dir']),
                  hovertext='Were the NP corrections from theory ("Calculated") or varied as part of the optimization ("Fit")?'),
        FieldSpec('Reduced Chi-Squared', CastProcessor(float, key='Reduced Chi-Squared', allows_empty=True)),
        notes_field
    ],
    data_key=lambda values: '_'.join([values['Reference'], 'radius', values['Nuclide'],
                                      'NP', values['Nuclear Polarization Method']])
)

muonic_fermi_distribution_template = InputTemplate(
    name="Fermi Distribution",
    fields=[
        reference_field,
        FieldSpec("Data used as input", PreviousDataProcessor(),
                  PreviousDataWidgetCreator(compilation_path=config['compilation_dir'],
                                            filter_regex='muonic.*_')),
        nuclide_field,
        # TODO add fermi distribution processor and widget creator with different fermi distribution options.
        FieldSpec('c [fm]', NumberWithUncertaintyProcessor(key='c [fm]')),
        FieldSpec('a [fm]', NumberWithUncertaintyProcessor(key='a [fm]')),
        FieldSpec('beta2 [-]', NumberWithUncertaintyProcessor(key='beta2 [-]')),
        FieldSpec('beta4 [-]', NumberWithUncertaintyProcessor(key='beta4 [-]')),
        notes_field
    ],
    data_key=lambda values: '_'.join([values['Reference'], 'fermi', values['Nuclide']])
)

templates = [muonic_transition_energy_template,
             muonic_transition_energy_difference_template,
             muonic_transition_energy_difference_diff_transition_template,
             muonic_nuclear_polarization_calculation_template,
             muonic_qed_calculation_template,
             muonic_barret_theory_template,
             muonic_barret_shift_template,
             muonic_radius_template,
             muonic_fermi_distribution_template]


class DataEntryInterface:
    templates: list[InputTemplate]
    current_entries: OrderedDict[str, tuple[tk.Widget, FieldProcessor]]
    references: dict[str, Reference]
    reference_path: Traversable | None
    compilation_folder: Traversable

    def __init__(self,
                 references: dict[str, Reference] | Traversable | None = None,
                 compilation_folder: Optional[str | Traversable] = None,
                 start_interface: bool = True):

        if references is None:
            self._load_references_from_json(config['reference_path'])
        elif isinstance(references, Traversable):
            self._load_references_from_json(references)
        elif isinstance(references, dict):
            self.reference_path = None
            self.references = references

        compilation_folder = compilation_folder if compilation_folder is not None else config['compilation_dir']
        self.compilation_folder = Path(compilation_folder) if not isinstance(compilation_folder,
                                                                             Traversable) else compilation_folder

        # Register templates here
        self.templates = templates

        self.field_widgets = {}  # field_name -> widget
        self.keep_vars = {}  # field_name -> BooleanVar

        self.current_entries: OrderedDict[
            str, tuple[tk.Widget, FieldProcessor]] = OrderedDict()  # Stores tk entry objects for each field displayed.

        if start_interface:
            self.root = tk.Tk()
            self.root.minsize(650, 500)
            self.root.title("Data Entry Interface")

            def _on_mousewheel(event):
                widget = event.widget
                try:
                    widget.yview_scroll(int(-1 * (event.delta / 120)), "units")
                except AttributeError:
                    pass

            self.root.bind_all("<MouseWheel>", _on_mousewheel)

            self.last_submitted = tk.StringVar()

            self.show_template_selection()

    def _load_references_from_json(self, reference_path: Optional[str | Traversable] = None):
        with open(reference_path, "r", encoding="utf-8") as f:
            raw_refs = json.load(f)

        self.reference_path = reference_path
        self.references = {key: Reference(**value) for key, value in raw_refs.items()}

    def clear_window(self):
        for widget in self.root.winfo_children():
            widget.destroy()

    def show_template_selection(self):
        self.clear_window()

        ttk.Label(
            self.root,
            text="Select Data Entry Template",
            font=("Arial", 16, "bold")
        ).pack(pady=10)

        self.template_var = tk.StringVar()

        names = [t.name for t in self.templates]

        template_combo = ttk.Combobox(
            self.root,
            textvariable=self.template_var,
            values=names,
            font=("Arial", 12),
            state="readonly",
            width=60
        )
        template_combo.pack(pady=5)

        ttk.Button(
            self.root,
            text="Continue",
            command=self.load_selected_template
        ).pack(pady=10)

    def load_selected_template(self):
        selected_name = self.template_var.get()
        if not selected_name:
            messagebox.showerror("Error", "Please select a template.")
            return

        for template in self.templates:
            if template.name == selected_name:
                self.display_input(template)
                return

    def display_input(self, template: InputTemplate):
        self.clear_window()
        self.current_entries.clear()

        ttk.Label(
            self.root,
            text=template.name,
            font=("Arial", 16, "bold"),
            wraplength=300
        ).pack(pady=10)

        form_frame = ScrollableFrame(self.root)
        form_frame.pack(padx=10, pady=10, expand=True, fill='both')

        form_frame = form_frame.scrollable_frame
        form_frame.columnconfigure(1, weight=1)

        for row, field in enumerate(template.fields):
            label = ttk.Label(form_frame,
                              text=field.label,
                              font=("Arial", 11),
                              wraplength=120
                              )
            label.grid(row=row, column=0, sticky='w', pady=4)
            if field.hovertext:
                ToolTip(label, field.hovertext)

            entry = field.widget_creator.create_widget(form_frame, row)
            self.current_entries[field.label] = (entry, field.processor)

            keep_var = tk.BooleanVar(value=False)
            self.keep_vars[field.label] = keep_var

            keep_check = ttk.Checkbutton(
                form_frame,
                text="Keep",
                variable=keep_var,
                takefocus=0
            )
            keep_check.grid(row=row, column=2, padx=(6, 0))

        button_frame = ttk.Frame(self.root)
        button_frame.pack(pady=10)

        submit = ttk.Button(
            button_frame,
            text="Submit",
            command=lambda: self.submit_data(template)
        )
        submit.grid(row=0, column=0, padx=5)
        submit.bind("<Return>", lambda _: self.submit_data(template))

        ttk.Button(
            button_frame,
            text="Back",
            command=self.show_template_selection
        ).grid(row=0, column=1, padx=5)

        ttk.Button(
            button_frame,
            text="2D",
            command=lambda: self.select_for_twod_input(template)
        ).grid(row=0, column=2, padx=5)

        ttk.Label(
            button_frame,
            textvariable=self.last_submitted,
            font=("Arial", 10)
        ).grid(row=1, column=0, pady=3)

    def select_for_twod_input(self, template):
        dialog = tk.Toplevel(self.root)
        dialog.title("2D Input Selection")
        dialog.transient(self.root)
        dialog.grab_set()

        ttk.Label(
            dialog,
            text="Select the fields you would like to form the row index, column index, and the value you want to enter "
                 "for each row/column pair.",
            wraplength=400,
            justify="left"
        ).pack(padx=20, pady=15)

        field_names = [field.label for field in template.fields]
        btn_frame = ttk.Frame(dialog)
        btn_frame.pack(fill="both", expand=True, padx=20, pady=10)
        buttons = []

        for i, selection in enumerate(['Row Field', 'Column Field', 'Value Field']):
            tk.Label(btn_frame, text=selection).grid(row=i, column=0)
            button = SelectWidgetCreator(field_names).create_widget(btn_frame, row=i)
            buttons.append(button)

        def create_2d():
            field_answers = [btn.get() for btn in buttons]
            dialog.destroy()
            self.twod_input(template, *field_answers)

        ttk.Button(btn_frame, text="Create 2D",
                   command=create_2d).grid(row=3, column=0, padx=5)
        ttk.Button(btn_frame, text="Back", command=lambda: dialog.destroy()).grid(row=3, column=1, padx=5)

        self.root.wait_window(dialog)

    def twod_input(self, template, row_field, col_field, value_field):
        dialog = tk.Toplevel(self.root)
        dialog.title(f"2D Input Selection: Values={value_field}")
        dialog.transient(self.root)
        dialog.grab_set()

        frame = ttk.Frame(dialog)
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        ttk.Label(frame, text=row_field).grid(row=1, column=0, padx=5, pady=5)
        ttk.Label(frame, text=col_field).grid(row=0, column=1, padx=5, pady=5)

        grid_frame = ttk.Frame(frame)
        grid_frame.grid(row=1, column=1)

        row_btn_frame = ttk.Frame(frame)
        row_btn_frame.grid(row=2, column=0, sticky="n")
        ttk.Button(row_btn_frame, text="+", command=lambda: grow_grid(rows=1, columns=0)).pack()
        ttk.Button(row_btn_frame, text="-", command=lambda: shrink_grid(rows=1, columns=0)).pack()

        col_btn_frame = ttk.Frame(frame)
        col_btn_frame.grid(row=0, column=2, sticky="e")
        ttk.Button(col_btn_frame, text="+", command=lambda: grow_grid(rows=0, columns=1)).pack(side="left")
        ttk.Button(col_btn_frame, text="-", command=lambda: shrink_grid(rows=0, columns=1)).pack(side="left")

        main_btn_frame = ttk.Frame(frame)
        main_btn_frame.grid(row=0, column=0, sticky='w')
        tk.Button(main_btn_frame, text='Submit',
                  command=lambda: self.submit_2d(template, row_field, col_field, value_field, row_widgets, col_widgets,
                                                 cells)).pack()
        tk.Button(main_btn_frame, text='Back', command=lambda: dialog.destroy()).pack()

        row_widgets: list[tk.Widget] = []
        col_widgets: list[tk.Widget] = []
        cells: list[list[tk.Widget]] = []

        def make_widget(field_name):
            field, = [field for field in template.fields if field.label == field_name]
            widget = field.widget_creator.create_widget(grid_frame, row=0)
            widget.config(width=20)
            return widget

        def regrid():
            for w in grid_frame.winfo_children():
                w.grid_forget()

            # Column headers
            for j, cw in enumerate(col_widgets):
                cw.grid(row=0, column=j + 1, padx=5, pady=5, sticky="nsew")
                grid_frame.columnconfigure(j + 1, weight=1)

            # Row headers + cells
            for i, rw in enumerate(row_widgets):
                rw.grid(row=i + 1, column=0, padx=5, pady=5, sticky="nsew")
                grid_frame.rowconfigure(i + 1, weight=1)

                for j, cell in enumerate(cells[i]):
                    cell.grid(row=i + 1, column=j + 1, padx=5, pady=5, sticky="nsew")
                    cell.bind("<Tab>", lambda e, fi=i, fj=j: (cells[fi][(fj + 1) % len(cells[fi])].focus(), "break")[1])
                    cell.bind("<Return>", lambda e, fi=i, fj=j: (cells[(fi + 1) % len(cells)][fj].focus(), "break")[1])

        def grow_grid(rows, columns):
            if rows:
                rw = make_widget(row_field)
                rw.config(takefocus=0)
                row_widgets.append(rw)

                new_row = [make_widget(value_field) for _ in col_widgets]
                cells.append(new_row)

            if columns:
                cw = make_widget(col_field)
                cw.config(takefocus=0)
                col_widgets.append(cw)

                for row in cells:
                    row.append(make_widget(value_field))

            regrid()

        def shrink_grid(rows, columns):
            if rows and row_widgets:
                row_widgets.pop().destroy()
                for cell in cells.pop():
                    cell.destroy()

            if columns and col_widgets:
                col_widgets.pop().destroy()
                for row in cells:
                    row.pop().destroy()

            regrid()

        grow_grid(rows=2, columns=2)

    def submit_2d(self, template, row_field, col_field, value_field, row_widgets, col_widgets, cells):
        # We take each (row_widget, col_widget, cell_widget) triplet, fill them in to the main data entry tab,
        # and trigger the data submission.

        original_entries = {field_name: self.current_entries[field_name] for field_name in
                            (row_field, col_field, value_field)}

        for row_widget, cell_row in zip(row_widgets, cells):
            for col_widget, cell in zip(col_widgets, cell_row):

                try:
                    if not cell.get():
                        continue
                except AttributeError:
                    continue

                for field_name, widget in [(row_field, row_widget), (col_field, col_widget), (value_field, cell)]:
                    prev_entry, processor = self.current_entries[field_name]
                    self.current_entries[field_name] = (widget, processor)
                self.submit_data(template,
                                 clear_fields=False)  # We don't want to erase our row/col entries as we submit

        for field_name in (row_field, col_field,
                           value_field):  # make sure we restore the original widgets once we are done submitting.
            self.current_entries[field_name] = original_entries[field_name]

        pass

    def submit_data(self, template: InputTemplate, clear_fields=True):

        def clear_widget(widg: tk.Widget):
            if isinstance(widg, ttk.Entry):
                widg.delete(0, tk.END)
                return

            # Composite widgets (e.g. Transition frame)
            for child in widget.winfo_children():
                if isinstance(child, ttk.Entry):
                    child.delete(0, tk.END)

        data = OrderedDict()

        try:
            for field_name, (widget, processor) in self.current_entries.items():
                result = processor.process(widget)
                if overlap := set(data.keys()).intersection(set(result.keys())):  # Key collision
                    raise ValueError(
                        f'Tried to add {field_name} data, but {overlap} field(s) are already present from another widget')
                data = data | result

        except ValueError as e:
            messagebox.showerror("Invalid Input", str(e))
            return

        self.save_data(template, data)

        focus_widget = None
        for field_name, (widget, _) in self.current_entries.items():
            if not self.keep_vars[field_name].get() and clear_fields:
                clear_widget(widget)
                if focus_widget is None:
                    widget.focus()
                    focus_widget = widget

        self.last_submitted.set(f'Last Submitted: {template.data_key(data)}')

    def save_data(
            self,
            template: InputTemplate,
            data: dict,
            replacement_strategy: Literal["Ask", "AlwaysReplace", "NeverReplace"] = "Ask",
    ) -> None:

        data_key = template.data_key(data)
        file_dir = self.compilation_folder / data["Reference"]
        file_dir.mkdir(parents=True, exist_ok=True)

        file = file_dir / f"{template.name}.json"

        new_entry = {data_key: data}

        if file.exists():
            with file.open("r", encoding="utf-8") as f:
                existing_data = json.load(f)

            if data_key in existing_data:
                match replacement_strategy:
                    case "AlwaysReplace":
                        existing_data[data_key] = data

                    case "NeverReplace":
                        return

                    case "Ask":
                        action = self._ask_collision_action(self.root, data_key, file)

                        if action == "replace":
                            existing_data[data_key] = data

                        elif action == "suffix":
                            new_key = self._next_available_key(data_key, existing_data)
                            existing_data[new_key] = data

                        else:  # do nothing
                            return
            else:
                existing_data[data_key] = data

            data_to_write = existing_data

        else:
            data_to_write = new_entry

        with file.open("w", encoding="utf-8") as f:
            json.dump(data_to_write, f, indent=4)

        # print(f"Data successfully saved to {file}")

    @staticmethod
    def _next_available_key(base_key: str, existing: dict) -> str:
        suffix = ord("a")
        while True:
            candidate = f"{base_key}_{chr(suffix)}"
            if candidate not in existing:
                return candidate
            suffix += 1

    @staticmethod
    def _ask_collision_action(parent, key: str, file: Path) -> str | None:
        """
        Returns: 'replace', 'suffix', or None (do nothing)
        """
        result = None
        dialog = tk.Toplevel(parent)
        dialog.title("Duplicate Entry")
        dialog.transient(parent)
        dialog.grab_set()

        ttk.Label(
            dialog,
            text=f"'{key}' already exists in:\n{file}\n\nWhat would you like to do?",
            wraplength=400,
            justify="left"
        ).pack(padx=20, pady=15)

        def choose(action):
            nonlocal result
            result = action
            dialog.destroy()

        btns = ttk.Frame(dialog)
        btns.pack(pady=10)

        ttk.Button(btns, text="Replace", command=lambda: choose("replace")).grid(row=0, column=0, padx=5)
        ttk.Button(btns, text="Add Suffix", command=lambda: choose("suffix")).grid(row=0, column=1, padx=5)
        ttk.Button(btns, text="Do Nothing", command=lambda: choose(None)).grid(row=0, column=2, padx=5)

        parent.wait_window(dialog)
        return result


if __name__ == "__main__":
    interface = DataEntryInterface(references=config['reference_path'])
    interface.root.mainloop()
