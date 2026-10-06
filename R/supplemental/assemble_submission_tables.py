#!/usr/bin/env python3
"""Assemble final Tables 1–12 from frozen Tables 1–11 and marker annotations.

No models are fitted and no numerical results are recalculated. The input
workbook must contain the reviewed Tables 1–11 and may already contain the two
Table 12 annotation sheets. Existing files are never overwritten.
"""
import argparse
import importlib.util
from pathlib import Path
import sys
import tempfile
import zipfile
import xml.etree.ElementTree as ET

CORE = ['ST1_Client_inventory', 'ST2_Pathway_coverage', 'ST3_Fig5_null_input',
        'ST4_Fig5_null_summary', 'ST5_Fig6_client_priority', 'ST6_Fig6_layer_summary',
        'ST7_Sensitivity_checks', 'ST8_Demographics_APOE', 'ST9_Conventional_DE_DA',
        'ST10_Pathology_models', 'ST11_Network_cognition']
ANNOTATIONS = ['S12_Marker_annotations', 'S12_GO_annotations']
INTERNAL = {'Source_notes', 'S12_Definitions'}
M = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'


def sheets(path):
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read('xl/workbook.xml'))
    return [n.get('name') for n in root.find('{'+M+'}sheets')]


def invoke(name, arguments, remove=None):
    path = Path(__file__).with_name(name + '.py')
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if remove is not None:
        module.REMOVE = remove
    original = sys.argv
    try:
        sys.argv = [str(path)] + [str(x) for x in arguments]
        module.main()
    finally:
        sys.argv = original


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--existing', required=True, type=Path)
    parser.add_argument('--table12', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists() or args.output.resolve() in {args.existing.resolve(), args.table12.resolve()}:
        parser.error('Choose a new output filename, distinct from both inputs.')
    names = sheets(args.existing)
    if [n for n in names if n not in set(ANNOTATIONS) | INTERNAL] != CORE:
        parser.error('Input must contain the reviewed Tables 1–11 in order; historical APOE appendix workbooks are not accepted.')
    with tempfile.TemporaryDirectory() as temporary:
        directory = Path(temporary)
        base = args.existing
        remove = set(names) & (set(ANNOTATIONS) | INTERNAL)
        if remove:
            base = directory / 'tables_1_11.xlsx'
            invoke('clean_submission_workbook', ['--input', args.existing, '--output', base], remove)
        appended = directory / 'tables_with_annotations.xlsx'
        invoke('add_supplementary_table_12', ['--existing', base, '--table12', args.table12, '--output', appended])
        invoke('clean_submission_workbook', ['--input', appended, '--output', args.output])
    if sheets(args.output) != CORE + ANNOTATIONS:
        raise ValueError('Unexpected submission workbook sheet order.')
    print('Verified: final workbook contains Tables 1–12 and no internal-note sheets.')


if __name__ == '__main__':
    main()
