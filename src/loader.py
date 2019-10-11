"""
Basic module for loading 1000 genomes data
"""
import shutil
import pathlib
import subprocess
import urllib.request as request

from contextlib import contextmanager, closing
from ftplib import FTP


BASE_PATH = pathlib.Path(__file__).parents[1]


class Loader:

    def __init__(self, base_url=None, base_path=None):
        self.base_url = base_url or 'ftp.1000genomes.ebi.ac.uk/vol1/ftp/'
        self.base_path = base_path or BASE_PATH

        self._set_credentials()
        self._configure()

    def _configure(self,  collection_path=None, file_pattern=None):
        self.collection_path = self.base_url + (collection_path or 'release/20130502/')

        default_pattern = 'ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
        self.file_pattern = file_pattern or default_pattern

    def local_path(self, filepath):
        filename = filepath.split('/')[-1]
        return str(self.base_path / 'data' / filename)

    def load(self, filepath):
        with closing(request.urlopen(self.base_url + filepath)) as remote:
            with open(self.local_path(filepath), 'wb') as f:
                shutil.copyfileobj(remote, f)

    def _set_credentials(self, login=None, password=None):
        self.login = login
        self.password = password

    @contextmanager
    def ftp_keeper(self, url):
        ftp = FTP(url)
        ftp.login(user=self.login or '', passwd=self.password or '')
        yield ftp
        ftp.quit()

    def ftp_load(self, filepath):
        with self.ftp_keeper(self.base_url) as ftp, open(self.local_path(filepath), 'wb') as f:
            ftp.retrbinary('RETR ' + filepath, f.write)

    def slice(self, chrom, start, end):
        file = f'ftp://{self.collection_path}{self.file_pattern.format(chrom=chrom)}'
        query = f'{chrom}:{start}-{end}'
        local_filename = f'chr{chrom}.{start}_{end}.vcf.gz'

        load_command = ' '.join(['tabix', '-h', file, query, '|', 'bgzip', '-c', '>', local_filename])
        subprocess.call(load_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        index_command = ' '.join(['tabix', '-p', 'vcf', local_filename])
        subprocess.call(index_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)


