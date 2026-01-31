import json
# Allows for
import os
import pathlib
import tkinter as tk
from tkinter import ttk

import pytest

from radbase.data_entry import (CastProcessor, DataEntryInterface, FieldSpec,
                                InputTemplate, NuclideProcessor,
                                NumberWithUncertaintyProcessor,
                                ReferenceProcessor, templates)

if os.environ.get('DISPLAY','') == '':
    print('no display found. Using :0.0')
    os.environ.__setitem__('DISPLAY', ':0.0')

example_references = {
    "Smith_nuclear_1978": {
        "title": "Nuclear Charge Radii from Muonic Atoms",
        "author": ["John Smith", "Alice Jones"],
        "year": "1978",
    },
    "Brown_isotope_1992": {
        "title": "Isotope Shifts in Heavy Nuclei",
        "author": ["Michael Brown"],
        "year": "1992"
    }
}

# Maps field names in data_entry.templates to example outputs.
example_outputs = {
    'Reference': list(example_references.keys())[0],
    'Previous Muonic Measurements': list(example_references.keys())[0] + '_muonic_2p1/2_1s1/2',
    'Nuclide': 'Pb208'
    ''


}


def make_dummy_references(tmp_path):
    temp_references_file = tmp_path / 'references.json'
    temp_references_file.write_text(json.dumps(example_references))
    return temp_references_file


class DummyProcessor:
    @staticmethod
    def process(widget):
        return 'nothing'


dummy_template = InputTemplate(name='test',
                               fields=[FieldSpec('Reference', CastProcessor(str)),
                                       FieldSpec('testfield', DummyProcessor())],
                               data_key=lambda x: 'testkey')


class DummyWidget(ttk.Widget):
    def __init__(self, value):
        self.value = value

    def get(self):
        return self.value


def test_cast_processor():
    widget = DummyWidget('30.2')
    blank_widget = DummyWidget('')

    assert CastProcessor(str).process(widget) == '30.2'
    assert CastProcessor(float).process(widget) == 30.2

    with pytest.raises(ValueError) as e_info:  # should raise error if it can not cast
        CastProcessor(int).process(widget)
    with pytest.raises(ValueError) as e_info:  # should raise error if field was left blank.
        CastProcessor(str).process(blank_widget)


def test_nuclide_processor():
    widget = DummyWidget('Pb208')
    natwidget = DummyWidget('Pbnat')
    bad_widget = DummyWidget('082203')

    assert NuclideProcessor().process(widget) == 'Pb208'
    assert NuclideProcessor().process(natwidget) == 'Pbnat'
    with pytest.raises(ValueError) as e_info:
        NuclideProcessor().process(bad_widget)


def test_reference_processor(tmp_path):

    global example_references

    good_widget = DummyWidget('Smith_nuclear_1978')
    bad_widget = DummyWidget('Noone_title_1340')

    tmp_ref_path = make_dummy_references(tmp_path)
    ref_proc = ReferenceProcessor(tmp_ref_path)

    assert ref_proc.process(good_widget) == 'Smith_nuclear_1978'
    with pytest.raises(ValueError) as e_info:
        ref_proc.process(bad_widget)


def test_number_with_uncertainty_processor():

    float_widget = DummyWidget('100.1')
    uncertainty_widget = DummyWidget('100.1(3)')
    bad_widget = DummyWidget('bad')

    ref_proc = NumberWithUncertaintyProcessor()

    assert ref_proc.process(float_widget) == {'value': 100.1, 'uncertainty': None}
    assert ref_proc.process(uncertainty_widget) == {'value': 100.1, 'uncertainty': 0.3}
    with pytest.raises(ValueError) as e_info:
        ref_proc.process(bad_widget)


def test_init_data_entry_interface(tmp_path):

    global example_references

    temp_references_file = tmp_path / 'references.json'
    temp_references_file.write_text(json.dumps(example_references))

    DataEntryInterface()


def test_save_data(tmp_path):

    ref_key = list(example_references.keys())[0]
    test_data = {'Reference': ref_key, 'testfield': 3}
    tmp_path.mkdir(exist_ok=True)
    tmp_data_loc = tmp_path / ref_key / f'{dummy_template.name}.json'

    interface = DataEntryInterface(compilation_folder=tmp_path)

    interface.save_data(dummy_template, test_data, replacement_strategy='AlwaysReplace')
    assert json.loads(tmp_data_loc.read_text()) == {'testkey': test_data}

    new_test_data = test_data.copy()
    new_test_data['testfield'] = 5
    interface.save_data(dummy_template, new_test_data, replacement_strategy='AlwaysReplace')
    assert json.loads(tmp_data_loc.read_text()) == {'testkey': new_test_data}

    interface.save_data(dummy_template, test_data, replacement_strategy='NeverReplace')
    assert json.loads(tmp_data_loc.read_text()) == {'testkey': new_test_data}


# def test_previous_measurement_processor(tmp_path):
#
#     previous_measurements = {'Authors'}
#     proc = PreviousMeasurementProcessor(previous_measurements = {})
