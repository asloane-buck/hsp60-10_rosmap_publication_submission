#!/usr/bin/env python3
"""Append annotation worksheets without reserializing existing worksheet XML.

Uses only the Python standard library. Writes a new workbook; never overwrites
the input. Existing worksheet bytes and scientific values are preserved.
"""
import argparse
import copy
from pathlib import Path
import posixpath
import zipfile
import xml.etree.ElementTree as ET

MAIN = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
REL = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
PKG = 'http://schemas.openxmlformats.org/package/2006/relationships'
CT = 'http://schemas.openxmlformats.org/package/2006/content-types'
ET.register_namespace('', MAIN)
ET.register_namespace('r', REL)

def xmlbytes(root):
    return ET.tostring(root, encoding='utf-8', xml_declaration=True)

def target_path(target):
    return target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/' + target)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--existing', required=True, type=Path)
    parser.add_argument('--table12', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.output.resolve() in {args.existing.resolve(), args.table12.resolve()}:
        parser.error('Output must be a new file, different from both inputs.')
    if args.output.exists():
        parser.error('Output already exists. Choose another filename or remove that output first.')
    with zipfile.ZipFile(args.existing) as z:
        original = {n: z.read(n) for n in z.namelist()}
    with zipfile.ZipFile(args.table12) as z:
        source = {n: z.read(n) for n in z.namelist()}
    wb = ET.fromstring(original['xl/workbook.xml'])
    rels = ET.fromstring(original['xl/_rels/workbook.xml.rels'])
    content = ET.fromstring(original['[Content_Types].xml'])
    sheets = wb.find(f'{{{MAIN}}}sheets')
    names = {s.get('name') for s in sheets}
    additions = [('Marker panel', 'S12_Marker_annotations'),
                 ('GO annotations', 'S12_GO_annotations'),
                 ('Definitions', 'S12_Definitions')]
    if names.intersection(n for _, n in additions):
        parser.error('Table 12 worksheets already exist in the destination workbook.')
    src_wb = ET.fromstring(source['xl/workbook.xml'])
    src_rels = {r.get('Id'): target_path(r.get('Target'))
                for r in ET.fromstring(source['xl/_rels/workbook.xml.rels'])}
    src_sheets = {s.get('name'): s for s in src_wb.find(f'{{{MAIN}}}sheets')}
    strings = []
    if 'xl/sharedStrings.xml' in source:
        strings = list(ET.fromstring(source['xl/sharedStrings.xml']))
    result = dict(original)
    next_id = max(int(s.get('sheetId')) for s in sheets) + 1
    used_rel_ids = {r.get('Id') for r in rels}
    counts = []
    for oldname, newname in additions:
        if oldname not in src_sheets:
            raise ValueError('Missing source sheet: ' + oldname)
        sheet = src_sheets[oldname]
        tree = ET.fromstring(source[src_rels[sheet.get(f'{{{REL}}}id')]])
        # New sheets use default formatting; existing workbook styles are untouched.
        # Resolve source shared strings to inline strings to avoid altering the
        # destination shared-string table or existing string indices.
        for cell in tree.iter(f'{{{MAIN}}}c'):
            cell.attrib.pop('s', None)
            if cell.get('t') == 's':
                value = cell.find(f'{{{MAIN}}}v')
                item = strings[int(value.text)]
                cell.remove(value)
                cell.set('t', 'inlineStr')
                inline = ET.SubElement(cell, f'{{{MAIN}}}is')
                for child in item:
                    if child.tag in {f'{{{MAIN}}}t', f'{{{MAIN}}}r'}:
                        inline.append(copy.deepcopy(child))
        for tag in ['tableParts', 'drawing', 'legacyDrawing', 'hyperlinks',
                    'conditionalFormatting', 'extLst']:
            for node in tree.findall(f'{{{MAIN}}}{tag}'):
                tree.remove(node)
        # Column and row widths, wrapping and freeze panes are kept, but source
        # style indices must not refer to the other workbook's style table.
        for node in tree.iter():
            if node.tag in {f'{{{MAIN}}}row', f'{{{MAIN}}}col'}:
                node.attrib.pop('s', None)
                node.attrib.pop('style', None)
                node.attrib.pop('customFormat', None)
        for view in tree.iter(f'{{{MAIN}}}sheetView'):
            view.attrib.pop('tabSelected', None)
        path = f'xl/worksheets/table12_{next_id}.xml'
        rid = 'rIdTable12_' + str(next_id)
        while rid in used_rel_ids:
            rid += '_'
        used_rel_ids.add(rid)
        ET.SubElement(sheets, f'{{{MAIN}}}sheet',
                      {'name': newname, 'sheetId': str(next_id), f'{{{REL}}}id': rid})
        ET.SubElement(rels, f'{{{PKG}}}Relationship',
                      {'Id': rid, 'Type': REL + '/worksheet',
                       'Target': 'worksheets/' + path.rsplit('/', 1)[1]})
        ET.SubElement(content, f'{{{CT}}}Override',
                      {'PartName': '/' + path,
                       'ContentType': 'application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml'})
        result[path] = xmlbytes(tree)
        counts.append((newname, len(tree.findall(f'{{{MAIN}}}sheetData/{{{MAIN}}}row')) - 1))
        next_id += 1
    assert counts[0][1] == 44, counts
    assert counts[1][1] == 1511, counts
    result['xl/workbook.xml'] = xmlbytes(wb)
    result['xl/_rels/workbook.xml.rels'] = xmlbytes(rels)
    result['[Content_Types].xml'] = xmlbytes(content)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(args.output, 'w', zipfile.ZIP_DEFLATED) as z:
        for name, payload in result.items():
            z.writestr(name, payload)
    changed = {'xl/workbook.xml', 'xl/_rels/workbook.xml.rels', '[Content_Types].xml'}
    with zipfile.ZipFile(args.output) as z:
        assert all(z.read(n) == data for n, data in original.items() if n not in changed)
    print('Saved:', args.output)
    print('Verified: all original worksheet XML, styles, shared strings and other original parts unchanged.')
    for name, count in counts:
        print(f'Added {name}: {count} records')
    print('The three added sheets use plain formatting. Their first rows contain filterable column headings.')

if __name__ == '__main__':
    main()
