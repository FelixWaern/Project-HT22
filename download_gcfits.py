# Downloading GC fits from SkewDB.
import ssl
from urllib.request import urlopen
from urllib import request
import zipfile

try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    # Legacy Python that doesn't verify HTTPS certificates by default
    pass
else:
    # Handle target environment that doesn't support HTTPS verification
    ssl._create_default_https_context = _create_unverified_https_context

def run_download_gcfits(file_path):
    # Downloading the gcfit file from SkewDB and unzipping it in the same location as the chromosome data. 

    gc_fits_zip_path = file_path+'gcfits.tar.bz2'
    zipurl = "https://berthub.eu/antonie/gcfits.tar.bz2"
    request.urlretrieve(zipurl, gc_fits_zip_path)
    with zipfile.ZipFile(gc_fits_zip_path, 'r') as zip_ref:
        zip_ref.extractall(file_path)

# Testing the script
csv_path = "C:/Users/Felix/Documents/FilteredDataFile.csv"
file_path = csv_path.rstrip("FilteredDataFile.csv")
run_download_gcfits(file_path)



