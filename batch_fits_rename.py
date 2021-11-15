#!/usr/bin/env python3

# %%
import argparse
import os
import re
from datetime import datetime, timezone
from astropy.io import fits
from glob import glob 

DESC="""
Batch rename FITS files in the specified path, using a template based on custom text and FITS keyword values.

Example: 
  batch_fits_rename.py e:\Astrofoto\_Flats\ '{FRAME}_{FILTER}_{DATEOBS}_{SEQ}{EXT}' -r
"""

RE_KEYWORD = r"\{(\w+)\}"

DEF_FITS_EXT = "fits"

# Debug
#sys.argv = ['this.py', 'E:\\Astrofoto\\_Flats', '{FRAME}_{FILTER}_{DATEOBS}_{SEQ}{EXT}', '-r', '-d']

parser = argparse.ArgumentParser(description=DESC, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("path_to_fits_files", 
                    help="Path to the FITS files")
parser.add_argument("filename_template", 
                    help="String that specifies the target filenames, combining fixed text and FITS keywords, " + 
                    "between curly braces. Also available as keywords are " + 
                    "DATEOBS (local datetime), SEQ (sequence number, if the original filename had it) " + 
                    "and EXT (original file extension, including the '.')")
parser.add_argument("-d", "--dry-run", action='store_true',
                    help="Dry run")
parser.add_argument("-r", "--recursive", action='store_true',
                    help="Recurse directories")
parser.add_argument("-e", "--extension", default=DEF_FITS_EXT, 
                    help="FITS extension, '%s' by default" % (DEF_FITS_EXT,))
args = parser.parse_args()


def get_fits_keyword_values(header, keys):
    return [str(header[k]) if k in hdr else 'MISSING' for k in key]


# FITS filename template (wildcard) from extension
fits_template = "*.%s" % (args.extension,)
recursive = args.recursive
if recursive:
    fits_template = os.path.join('**', fits_template)

# FITS path
fits_directory = args.path_to_fits_files
if not os.path.isdir(fits_directory):
    print("Path doesn't exist or is not valid!")
    exit(1)


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
        raise ValueError("No sequence value detected on filename!")

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



# %%
keys = get_keywords_in_tpl(args.filename_template)

fits_template_full_path = os.path.join(fits_directory, fits_template)
print(fits_template_full_path)

prev_dir = ''
for fname in glob(fits_template_full_path, recursive=recursive):
    with fits.open(fname) as hdul:
        hdr = dict(hdul[0].header)
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
