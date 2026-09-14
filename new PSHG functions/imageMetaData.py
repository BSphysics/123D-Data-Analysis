import os
import re
from ScanImageTiffReader import ScanImageTiffReader


def imageMetaData(data_path, files_tiff, out_name='Image meta data.txt'):
    """Read the ScanImage header from one representative TIFF and dump every
    metadata field (including empty ones) to a text file.

    files_tiff may be a list/tuple (first entry is used, since the whole
    acquisition shares one non-varying header) or a single filename string.
    A full path in place of a bare filename also works.
    """
    # --- resolve to a single file: works for a group OR a lone image ------
    first = files_tiff[0] if isinstance(files_tiff, (list, tuple)) else files_tiff
    if not first:
        raise ValueError('No file supplied to imageMetaData')
    filename = first if os.path.isabs(first) else os.path.join(data_path, first)

    # --- read metadata once -----------------------------------------------
    with ScanImageTiffReader(filename) as reader:
        md = reader.metadata()                 # non-varying SI header (+RoiGroups)
        try:
            desc0 = reader.description(0)       # per-frame tag of frame 0
        except Exception:
            desc0 = ''

    # --- parse "SI.* = value" lines into a dict (keeps empty RHS) ----------
    meta = {}
    for line in md.splitlines():
        if '=' in line:
            key, _, val = line.partition('=')
            meta[key.strip()] = val.strip()

    # --- write EVERYTHING to disk -----------------------------------------
    out_path = os.path.join(data_path, out_name)
    with open(out_path, 'w') as f:
        f.write('# ScanImage metadata for: %s\n\n' % filename)
        f.write('## Parsed header (all SI.* fields, empties included)\n')
        for key in sorted(meta):
            f.write('%s = %s\n' % (key, meta[key]))
        f.write('\n## Frame 0 per-frame description\n%s\n' % desc0)
        f.write('\n## Raw metadata string (verbatim)\n%s\n' % md)

    # --- robust lookups: match exact key, else the trailing token ----------
    num = re.compile(r'[-+]?(?:\d*\.\d+|\d+\.?)(?:[Ee][+-]?\d+)?')

    def find_key(short):
        if short in meta:
            return short
        cands = [k for k in meta if k.split('.')[-1] == short]
        if len(cands) > 1:
            print('WARNING: %r is ambiguous: %s' % (short, cands))
        return cands[0] if cands else None

    def get_number(short):
        key = find_key(short)
        if key is None:
            print('WARNING: %r not found in metadata' % short)
            return None
        m = num.search(meta[key])
        return float(m.group()) if m else None

    zoomFactor = get_number('scanZoomFactor')
    frames     = get_number('framesPerSlice')
    zStackStep = get_number('stackZStepSize')

    print(' Zoom factor = %s' % zoomFactor)
    print(' Number of frames = %s' % frames)

    return zoomFactor, frames, zStackStep, meta
