#!/usr/bin/env python3

# %%
import argparse
import sys
import os
import re
from datetime import datetime, timezone
from astropy.io import fits
from xisf import XISF
from glob import iglob 

DESC="""
Batch rename FITS files in the specified path, using a template based on custom text and FITS keyword values.

Example: 
  batch_fits_rename.py e:\Astrofoto\_Flats\**\*.fits '{FRAME}_{FILTER}_{DATEOBS}_{SEQ}{EXT}' -r
"""

# Debug
#sys.argv = ['this.py', 'E:\\Astrofoto\\_Flats\\**\\*.xisf', 'masterFlat_{DATEOBS}_BIN-1_{FILTER}_432mm_F4.8{EXT}', '-r', '-d']


RE_KEYWORD = r"\{(\w+)\}"


parser = argparse.ArgumentParser(description=DESC, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("path_to_fits_or_xisf_files", 
                    help="Path to the FITS or XISF files, using wildcards. With -r, you may also specify '**' " +
                    "to expand for any number of directory levels.")
parser.add_argument("filename_template", 
                    help="String that specifies the target filenames, combining fixed text and FITS keywords, " + 
                    "between curly braces. Also available as keywords are " + 
                    "DATEOBS (local datetime), SEQ (sequence number, if the original filename had it) " + 
                    "and EXT (original file extension, including the '.')")
parser.add_argument("-d", "--dry-run", action='store_true',
                    help="Dry run")
parser.add_argument("-r", "--recursive", action='store_true',
                    help="Recurse directories")
args = parser.parse_args()


def get_fits_keyword_values(header, keys):
    return [str(header[k]) if k in hdr else 'MISSING' for k in key]


# FITS path
files_path = args.path_to_fits_or_xisf_files
recursive = args.recursive


def get_file_keywords(fname):
    hdr = None

    try:
        with fits.open(fname) as hdul:
            hdr = dict(hdul[0].header)
    except:
        pass

    try:
        xisf = XISF(fname)
        hdr = xisf.get_images_metadata()[0]
        fits_keyw = hdr.pop('FITSKeywords')
        fits_keyw = {k:fits_keyw[k][0]['value'] for k in fits_keyw.keys()} # Only the first value
        xisf_keyw = hdr.pop('XISFProperties')
        hdr = {**hdr, **fits_keyw, **xisf_keyw} # TODO: possible collisions
    except:
        pass

    if hdr is None:
        print(f"Error: {fname} not recognized as valid FITS or XISF file!", file=sys.stderr)
        sys.exit(1)

    return hdr


def get_keywords_in_tpl(fname_tpl):
    return re.findall(RE_KEYWORD, fname_tpl)


# def format_fname_tpl(fname_tpl):
#     return re.sub(RE_KEYWORD, r"{hdr['\1']}", fname_tpl)


def augment_keywords(hdr, fname):
    # Add "keywords" not present in FITS header, based on current filename
    hdr['NAME'], hdr['EXT'] = os.path.splitext(os.path.basename(fname))
    try:
        hdr['SEQ'] = re.search(r'\d+$', hdr['NAME']).group()
    except:
        hdr['SEQ'] = "NULLSEQ"

    # Add DATEOBS as local time
    hdr['DATEOBS'] = datetime.fromisoformat(hdr['DATE-OBS']).replace(tzinfo=timezone.utc).astimezone(tz=None).isoformat()[:19]

    # Trim DATE-OBS (UTC) to second accuracy
    hdr['DATE-OBS'] = hdr['DATE-OBS'][:19]

    # If FILTER is missing (KStars/Ekos bug), use the current directory
    if 'FILTER' not in hdr:
        hdr['FILTER'] = os.path.basename(os.path.dirname(fname))
        print(f"  WARN in {hdr['NAME']}: FILTER not present, filling with current directory name ({hdr['FILTER']})")


def to_safe_fname(s):
    s = s.replace(':', '-')
    return "".join( c for c in s if (c.isalnum() or c in "._- "))


def to_filename(tpl, hdr):
    return to_safe_fname(tpl.format(**hdr))

    
keys = get_keywords_in_tpl(args.filename_template)
prev_dir = ''
fname_gen = ( fname for fname in iglob(files_path, recursive=recursive) if os.path.isfile(fname) )
for fname in fname_gen:
    hdr = get_file_keywords(fname)
    augment_keywords(hdr, fname)

    keys_not_present = set(keys) - set(hdr.keys())
    if keys_not_present:
        raise ValueError("Keywords not present! ", ", ".join(keys_not_present))

    cur_dir = os.path.dirname(fname)
    if cur_dir != prev_dir:
        print('\n\n__/', cur_dir, '\__________')
        prev_dir = cur_dir

    new_fname = to_filename(args.filename_template, hdr)        
    print(os.path.basename(fname), "->", new_fname, end=' ... ')
    new_fname = os.path.join(cur_dir, new_fname)

    if not args.dry_run:
        os.rename(fname, new_fname)
        print("done.")
    else:
        print("skipped.")



# %%
