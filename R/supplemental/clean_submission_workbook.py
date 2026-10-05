#!/usr/bin/env python3
"""Remove internal audit-note worksheets without changing scientific data.

Python standard library only. Creates a new workbook and preserves all retained
worksheet XML byte-for-byte. Sources and GO evidence in Table 12 are retained.
"""
import argparse
from pathlib import Path
import posixpath
import zipfile
import xml.etree.ElementTree as ET

M = 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'
R = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships'
P = 'http://schemas.openxmlformats.org/package/2006/relationships'
ET.register_namespace('', M)
ET.register_namespace('r', R)
REMOVE = {'Source_notes', 'S12_Definitions'}

def resolve(target):
    return target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/' + target)

def encode(tree):
    return ET.tostring(tree, encoding='utf-8', xml_declaration=True)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    if a.input.resolve() == a.output.resolve() or a.output.exists():
        ap.error('Choose a new output filename. Input and existing files will not be overwritten.')
    with zipfile.ZipFile(a.input) as z:
        data = {n: z.read(n) for n in z.namelist()}
    wb = ET.fromstring(data['xl/workbook.xml'])
    sheets = wb.find(f'{{{M}}}sheets')
    rels = ET.fromstring(data['xl/_rels/workbook.xml.rels'])
    rel_map = {r.get('Id'): r for r in rels}
    removed = []
    drop = set()
    retained = []
    old_to_new = {}
    for oldindex, s in enumerate(list(sheets)):
        rid = s.get(f'{{{R}}}id')
        path = resolve(rel_map[rid].get('Target'))
        if s.get('name') in REMOVE:
            removed.append(s.get('name'))
            sheets.remove(s)
            rels.remove(rel_map[rid])
            drop.add(path)
            rel_path = posixpath.dirname(path) + '/_rels/' + posixpath.basename(path) + '.rels'
            drop.add(rel_path)
        else:
            old_to_new[oldindex] = len(retained)
            retained.append((s.get('name'), path))
    if not removed:
        ap.error('Neither internal-notes sheet is present. No changes made.')
    # Do not silently break a retained scientific formula that refers to notes.
    for name, path in retained:
        for f in ET.fromstring(data[path]).iter(f'{{{M}}}f'):
            if any(n in (f.text or '') for n in removed):
                raise ValueError('Retained formula refers to a removed sheet: ' + name)
    defined = wb.find(f'{{{M}}}definedNames')
    if defined is not None:
        for n in list(defined):
            local = n.get('localSheetId')
            if any(t in (n.text or '') for t in removed) or (local is not None and int(local) not in old_to_new):
                defined.remove(n)
            elif local is not None:
                n.set('localSheetId', str(old_to_new[int(local)]))
    for v in wb.iter(f'{{{M}}}workbookView'):
        if v.get('activeTab'):
            v.set('activeTab', str(old_to_new.get(int(v.get('activeTab')), 0)))
        if v.get('firstSheet'):
            v.set('firstSheet', str(old_to_new.get(int(v.get('firstSheet')), 0)))
    ct = ET.fromstring(data['[Content_Types].xml'])
    for n in list(ct):
        if n.get('PartName', '').lstrip('/') in drop:
            ct.remove(n)
    result = {n: b for n, b in data.items() if n not in drop}
    result['xl/workbook.xml'] = encode(wb)
    result['xl/_rels/workbook.xml.rels'] = encode(rels)
    result['[Content_Types].xml'] = encode(ct)
    # Remove stale worksheet-name lists in extended file properties.
    if 'docProps/app.xml' in result:
        app = ET.fromstring(result['docProps/app.xml'])
        ns = 'http://schemas.openxmlformats.org/officeDocument/2006/extended-properties'
        vt = 'http://schemas.openxmlformats.org/officeDocument/2006/docPropsVTypes'
        titles = app.find(f'{{{ns}}}TitlesOfParts/{{{vt}}}vector')
        if titles is not None:
            for child in list(titles):
                if child.text in removed:
                    titles.remove(child)
            titles.set('size', str(len(titles)))
        hp = app.find(f'{{{ns}}}HeadingPairs/{{{vt}}}vector')
        if hp is not None:
            variants = list(hp)
            for i in range(0, len(variants)-1, 2):
                label = variants[i].find(f'{{{vt}}}lpstr')
                count = variants[i+1].find(f'{{{vt}}}i4')
                if label is not None and label.text == 'Worksheets' and count is not None:
                    count.text = str(len(retained))
        result['docProps/app.xml'] = encode(app)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(a.output, 'w', zipfile.ZIP_DEFLATED) as z:
        for n, b in result.items():
            z.writestr(n, b)
    with zipfile.ZipFile(a.output) as z:
        assert all(z.read(path) == data[path] for _, path in retained)
        assert z.read('xl/styles.xml') == data['xl/styles.xml']
        if 'xl/sharedStrings.xml' in data:
            assert z.read('xl/sharedStrings.xml') == data['xl/sharedStrings.xml']
    print('Saved:', a.output)
    print('Removed:', ', '.join(removed))
    print('Verified: every retained worksheet and its scientific values are unchanged.')
    print('Protein-source links and GO evidence columns are retained.')

if __name__ == '__main__':
    main()
