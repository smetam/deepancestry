"""
Basic module for loading 1000 genomes data
"""
import ftplib
import shutil
import pathlib
import subprocess
import urllib.request as request

from contextlib import closing


BASE_PATH = pathlib.Path(__file__).parents[1]


# TODO: add logging
class Loader:
    def __init__(self, chromosome, start, end, out_dir):
        self.start = start
        self.end = end
        self.chromosome = chromosome
        self._ensure_dir(out_dir)
        self.server = 'ftp.1000genomes.ebi.ac.uk'
        self.collection = 'vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/'
        self.chr_file = f'ALL.chr{chromosome}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz'

    def _ensure_dir(self, out_dir):
        self.out_dir = out_dir or f'chr{self.chromosome}_data'

        path = pathlib.Path(BASE_PATH / self.out_dir)
        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)

    @property
    def slice_mode(self):
        return self.start and self.end

    @property
    def query(self):
        return f'{self.chromosome}:{self.start}-{self.end}'

    @property
    def ftp_path(self):
        return f'ftp://{self.server}/{self.collection}{self.chr_file}'

    @property
    def local_filename(self):
        if self.slice_mode:
            return f'chr{self.chromosome}.{self.start}_{self.end}.vcf.gz'
        return f'chr{self.chromosome}_full.vcf.gz'

    @property
    def local_filepath(self):

        return str(BASE_PATH / self.out_dir / self.local_filename)

    def _fetch(self):
        ftp = ftplib.FTP(self.server)
        ftp.login()

        ftp.cwd(self.collection)
        ftp.retrlines('LIST')

        ftp.retrbinary("RETR {}".format(self.chr_file), open(self.local_filepath, 'wb').write)
        ftp.quit()

    def _slice(self):
        load_command = ' '.join([f'mkdir {self.out_dir} && cd {self.out_dir} &&', 'tabix', '-h',
                                 self.ftp_path, self.query, '|', 'bgzip', '-c', '>', self.local_filename])
        subprocess.call(load_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        index_command = ' '.join([f'cd {self.out_dir} &&', 'tabix', '-p', 'vcf', self.local_filename])
        subprocess.call(index_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    def _load(self):
        with closing(request.urlopen(self.ftp_path)) as remote:
            with open(self.local_filepath, 'wb') as f:
                shutil.copyfileobj(remote, f)

    def fetch(self):
        if self.slice_mode:
            self._slice()
        else:
            self._fetch()


