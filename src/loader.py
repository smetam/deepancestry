from ftplib import FTP


ftp = FTP('ftp.1000genomes.ebi.ac.uk')
ftp.login()

data = ftp.retrlines('LIST')
print(data)
ftp.cwd('vol1/ftp')
data = ftp.retrlines('LIST')


