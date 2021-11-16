#!/usr/bin/env python3

# %%
import argparse
import os
from astropy.io import fits
from glob import iglob 

DESC="""
Batch-process FITS files in the specified path.

Usage:
  -k <keyword> -v <value> [-c "Comment ..."] <path> : Adds or updates the value (and optionally, the comment) of the given keyword. 
  -k <keyword_1>,...,<keyword_N> <path>             : Prints the keywords values in CSV format

With -r, you may specify '**' in the <path> to expand for any number of directory levels.

Examples:
* Get the value of the FILTER keyword for all files named like '*Ha*.fits' under the given directory, recursively:
```
  batch_fits_keyword.py e:\astrofoto\_Flats\**\*Ha*.fits -r -k FILTER
```  

* Same as above, but setting its value to 'Ha':
```
  batch_fits_keyword.py e:\astrofoto\_Flats\**\*Ha*.fits -r -k FILTER -v 'Ha'  
```

"""

# Debug
#sys.argv = ['this.py', 'E:\\Astrofoto\\_Flats\\**\\*Ha*.fits', '-r', '-k FILTER']


KEY_LIST_SEP = ","

parser = argparse.ArgumentParser(description=DESC, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("path_to_fits_files", 
                    help="Path to the FITS files")
parser.add_argument("-r", "--recursive", action='store_true',
                    help="Recurse directories")
parser.add_argument("-k", "--key", required=True,
                    help="Name of the keyword(s) to print, or the keyword (only one) to add or modify (required)")
parser.add_argument("-v", "--value", 
                    help="Value for the keyword. Use single quotes for string literals (optional)")
parser.add_argument("-c", "--comment", default='',
                    help="Comment for the keyword/value (optional)")
args = parser.parse_args()


# FITS path
files_path = args.path_to_fits_files
recursive = args.recursive

# FITS keyword
key = args.key.split(KEY_LIST_SEP)

# FITS value, FITS comment
value, comment = args.value, args.comment
if value:
    if len(key) > 1:
        print("Multiple keywords given for update!")
        exit(1)

    # Guess value type
    if args.value[0] == "'":
        value = args.value[1:-1]
    else:
        try:
            value = int(args.value)
        except ValueError:
            value = float(args.value)


def get_keyword(fits_image_filename, key):
    with fits.open(fits_image_filename) as hdul:
        hdr = hdul[0].header
        return [str(hdr[k]) if k in hdr else 'MISSING!' for k in key]


def add_replace_keyword(fits_image_filename, key, value, comment):
    with fits.open(fits_image_filename, mode='update') as hdul:
        hdr = hdul[0].header
        if not comment and key in hdr.comments:
            comment = hdr.comments[key]
        hdr[key] = (value, comment)
        hdul.flush()


if not value:
    print(KEY_LIST_SEP.join(["Filename"] + key))
else:
    print(f"Adding or replacing keyword: {key} = {value} // {comment}")

fname_gen = ( fname for fname in iglob(files_path, recursive=recursive) if os.path.isfile(fname) )
for fname in fname_gen:
    current_value = get_keyword(fname, key)
    if not value:
        print(KEY_LIST_SEP.join([fname] + current_value))
    else:
        add_replace_keyword(fname, key[0], value, comment)
        print(f"{fname},{current_value[0]}->{value} ({comment})")



# %%
