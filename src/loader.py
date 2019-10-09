"""
Basic module for loading 1000 genomes data
"""
from contextlib import contextmanager, closing
import shutil
import urllib.request as request

from ftplib import FTP

from src.parser import BASE_PATH


class Loader:

    def __init__(self, base_url=None, base_path=None):
        self.base_url = base_url or 'ftp.1000genomes.ebi.ac.uk/vol1/ftp/'
        self.base_path = base_path or BASE_PATH

    def local_path(self, filepath):
        filename = filepath.split('/')[-1]
        return str(self.base_path / 'data' / filename)

    def load(self, filepath):
        with closing(request.urlopen(self.base_url + filepath)) as remote:
            with open(self.local_path(filepath), 'wb') as f:
                shutil.copyfileobj(remote, f)

    @contextmanager
    def ftp_keeper(self, url):
        ftp = FTP(url)
        yield ftp
        ftp.quit()

    def ftp_load(self, filepath):
        with self.ftp_keeper(self.base_url) as ftp, open(self.local_path(filepath), 'wb') as f:
            ftp.retrbinary('RETR ' + filepath, f.write)

