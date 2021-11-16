#!/usr/bin/env python3

# %%
import argparse
import sys
import os
import re
from datetime import datetime, timezone, time, timedelta
from astropy.io import fits
from xisf import XISF
from glob import iglob 

DESC="""
Batch rename FITS files in the specified path, using a template based on custom text and FITS keyword values.

Examples: 
- Basic usage:
  ```
  batch_fits_rename.py E:\Astro\Flats\**\*.fits '{FRAME}_{FILTER}_{DATEOBS}_{SEQ}{EXT}' -r`
  ```

- Renaming and organizing in subdirectories in one go:
  ```
  batch_fits_rename.py Light\**\*.fits "SESSION_{SESSION}\LionNeb_Light_{FILTER}_600secs_{LOC_ISO_DATETIME}_{SEQ}{EXT}" -r 
  ```
"""

# Debug
#sys.argv = ['this.py', 'E:\\Astrofoto\\_Flats\\**\\*.xisf', 'masterFlat_{DATEOBS}_BIN-1_{FILTER}_432mm_F4.8{EXT}', '-r', '-d']


RE_KEYWORD = r"\{(\w+)\}"


parser = argparse.ArgumentParser(description=DESC, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("path_to_fits_or_xisf_files", 
                    help="Path to the FITS or XISF files, using wildcards. With -r, you may also specify '**' " +
                    "to expand for any number of directory levels.")
parser.add_argument("filename_template", 
                    help="String that specifies the target filenames, combining fixed text and " + 
                    "FITS keywords between curly braces. Also available as keywords are " + 
                    "LOC_ISO_DATETIME (local datetime), LOC_ISO_DATE (local date), " +
                    "SESSION (date 'rounded' to the nearest midnight), " + 
                    "SEQ (sequence number, if the original filename had it) " + 
                    "and EXT (original file extension, including the '.')." +
                    "This template could include pathname separators: in this case, the directories " + 
                    "will be automatically created. See the examples.")
parser.add_argument("-d", "--dry-run", action='store_true',
                    help="Dry run")
parser.add_argument("-r", "--recursive", action='store_true',
                    help="Recurse directories")
args = parser.parse_args()

# FITS path
files_path = args.path_to_fits_or_xisf_files
recursive = args.recursive


def get_fits_keyword_values(header, keys):
    return [str(header[k]) if k in hdr else 'MISSING' for k in key]


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

    # Add local time
    dateobs_dt = datetime.fromisoformat(hdr['DATE-OBS']).replace(tzinfo=timezone.utc).astimezone(tz=None)
    dateobs_iso = dateobs_dt.isoformat()
    hdr['LOC_ISO_DATETIME'] = dateobs_iso[:19]
    hdr['LOC_ISO_DATE'] = dateobs_iso[:10]

    # Trim DATE-OBS (UTC) to second accuracy
    hdr['DATE-OBS'] = hdr['DATE-OBS'][:19]

    # Derive SESSION identifier from DATEOBS, "rounding" to the nearest midnight
    session_date = dateobs_dt if dateobs_dt.time() < time(12, 00) else dateobs_dt + timedelta(days=1)
    hdr['SESSION'] =  session_date.strftime("%Y%m%d")

    # If FILTER is missing (KStars/Ekos bug), use the current directory
    if 'FILTER' not in hdr:
        hdr['FILTER'] = os.path.basename(os.path.dirname(fname))
        print(f"  WARN in {hdr['NAME']}: FILTER not present, filling with current directory name ({hdr['FILTER']})")


def to_safe_fname(s):
    s = s.replace(':', '-')
    return "".join( c for c in s if (c.isalnum() or c in "._- \\"))


def to_filename(tpl, hdr):
    return to_safe_fname(tpl.format(**hdr))

    
keys = get_keywords_in_tpl(args.filename_template)
prev_dir = ''
# We snapshot the fnames here to avoid recursion (renaming the renamed files!)
fname_list = [fname for fname in iglob(files_path, recursive=recursive) if os.path.isfile(fname)]
for fname in fname_list:
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
        os.makedirs(os.path.dirname(new_fname), exist_ok=True)
        os.rename(fname, new_fname)
        print("done.")
    else:
        print("skipped.")



# %%
