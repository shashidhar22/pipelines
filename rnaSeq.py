import os
import re
import sys
import csv
import gzip
import time
import logging
import argparse
import subprocess
from collections import namedtuple

class QualCheck:

    def __init__(self, bbduk_path, adp_path, out_path):
        self.bbduk_path = os.path.abspath(bbduk_path)
        self.adp_path = os.path.abspath(adp_path)
        self.out_path = '{0}/CleanedFastq'.format(os.path.abspath(out_path))
        self.logger = logging.getLogger('RNASeq.BBDuk')
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)

    def bbduk(self, rone_path, rtwo_path):
        #java -ea -Xmx6000m -cp .../current jgi.BBDukF ktrim=r k=27 hdist=1
        #edist=0 mink=4 ref=adapters.fa qtrim=rl trimq=30 minlength=50 qin=33
        #in=input1.fastq in2=input2.fastq out=output1.fastq out2=output2.fastq
        #Setup output paths
        brone = os.path.splitext(os.path.basename(rone_path))[0]
        brtwo = os.path.splitext(os.path.basename(rtwo_path))[0]
        orone_path = '{0}/{1}_cleaned.fq'.format(self.out_path, brone)
        ortwo_path = '{0}/{1}_cleaned.fq'.format(self.out_path, brtwo)
        stats_path = '{0}/{1}_stats.txt'.format(self.out_path, brone)
        #Set up the command
        bbcmd = [self.bbduk_path, '-Xmx1g', 'k=27', 'hdist=1', 'edist=0', 'ktrim=l',
                'mink=4', 'ref={0}'.format(self.adp_path), 'qtrim=rl',
                'trimq=30', 'minlength=50', 'qin=33', 'overwrite=true',
                'in={0}'.format(rone_path), 'in2={0}'.format(rtwo_path),
                'out={0}'.format(orone_path), 'out2={0}'.format(ortwo_path),
                'stats={0}'.format(stats_path)]

        #Run bbduk
        bbrun = subprocess.Popen(bbcmd, shell=False)
        bbrun.wait()
        if bbrun.returncode != 0:
            self.logger.error('BBDuk failed running the following command : {0}'.format(' '.join(bbcmd)))
        return(orone_path, ortwo_path, bbrun.returncode)

class SplicedAligner:

    def __init__(self, tophat_path, ref_path, gtf_path, out_path):
        self.tophat_path = os.path.abspath(tophat_path)
        self.ref_path = os.path.abspath(ref_path)
        self.gtf_path = os.path.abspath(gtf_path)
        self.out_path = '{0}/MappedReads'.format(os.path.abspath(out_path))
        self.logger = logging.getLogger('RNASeq.Tophat2')
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)

    def tophat2(self, fastq_dict):
        #tophat2 --no-coverage-search --b2-very-sensitive --no-novel-juncs -p 4
        #-G <genome_gtf> -o <out_path> <bowtie_index> <input_path>
        tprone = ','.join(fastq_dict['R1'])
        tprtwo = ','.join(fastq_dict['R2'])
        #Setup the command
        tpcmd = [self.tophat_path, '--no-coverage-search', '--b2-very-sensitive',
                '--no-novel-juncs', '-p', '4', '-G', self.gtf_path, '-o',
                self.out_path, self.ref_path, tprone, tprtwo]

        #Run Tophat2
        tprun = subprocess.Popen(tpcmd, shell=False)
        tprun.wait()
        if tprun.returncode != 0:
            bam_path = None
            self.logger.error('Tophat2 failed running the following command : {0}'.format(' '.join(tpcmd)))
        else:
            bam_path = '{0}/accepted_hits.bam'.format(self.out_path)
            self.logger.info('Tophat2 completed successfully')
        return(bam_path, tprun.returncode)

class GeneQuant:

    def __init__(self, cuff_path, gtf_path, bam_path, out_path):
        self.cuff_path = os.path.abspath(cuff_path)
        self.gtf_path = os.path.abspath(gtf_path)
        self.bam_path = os.path.abspath(bam_path)
        self.out_path = '{0}/GeneExpression'.format(os.path.abspath(out_path))
        self.logger = logging.getLogger('RNASeq.Cufflinks')
        if not os.path.exists(self.out_path):
            os.mkdir(self.out_path)

    def cuffLinks(self):
        #cufflinks -p 4 -G <gtf_path> -o <out_path> <bam_file>
        #Setup the command
        clcmd = [self.cuff_path, '-p', '4', '-G', gtf_path, '-o', self.out_path,
                self.bam_path]

        #Run Cufflinks
        clrun = subprocess.Popen(clcmd, shell=False)
        clrun.wait()
        if clrun.returncode != 0:
            self.logger.error('Cufflinks failed running the following command : {0}'.format(' '.join(clcmd)))
        else:
            self.logger.info('Cufflinks completed successfully')
        return(clrun.returncode)

class Fastq:

    def __init__(self, fastqfile, outdir, phred):
        self.fastq = fastqfile
        self.outdir = outdir
        self.phred = phred
        FORMAT = '%(asctime)-15s : %(levelname)-8s :  %(message)s'
        logging.basicConfig(format=FORMAT)
        self.phreddict = self.phredmap()


    def phredmap(self):
        phreddict = dict()
        if self.phred == "phred33":
            for asciis, quals in zip(range(33,126), range(0,92)):
                phreddict[asciis] = quals
        return(phreddict)

    def formatChecker(self, header, seq, sheader, quals, line, phred):
        if header[0] != "@":
            return(1)
        if sheader == "" or sheader[0] != "+":
            return(2)
        if len(seq) != len(quals):
            return(3)
        if seq == '':
            return(4)

    def qualmasker(self, seq, quals):
        maskedseq = list()
        for base, qual in zip(seq, quals):
            if qual < 20:
                maskedseq.append('N')
            else:
                maskedseq.append(base)
        return(maskedseq, quals)

    def read(self):
        if '.gz' in self.fastq or '.fastqgz' in self.fastq :
            fqhandle = gzip.open(self.fastq, 'rb')
        else:
            fqhandle = open(self.fastq, 'r')
        Record = namedtuple('Record',['header', 'sheader', 'seq', 'quals'])
        lineno = 0
        while True:
            try:
                #Get fastq record
                if '.gz' in self.fastq or '.fastqgz' in self.fastq:
                    header = str(next(fqhandle),'utf-8').strip()
                    seq = [base for base in str(next(fqhandle),'utf-8').strip()]
                    sheader = str(next(fqhandle), 'utf-8').strip()
                    quals = str(next(fqhandle), 'utf-8').strip()
                else:
                    header = next(fqhandle).strip()
                    seq = [base for base in next(fqhandle).strip()]
                    sheader = next(fqhandle).strip()
                    quals = next(fqhandle).strip()
                #quals = [self.phreddict[int(ord(str(qual)))] for qual in next(fqhandle).strip()]
                lineno += 1

                #Check if record adheres to fastq format
                check = self.formatChecker(header, seq, sheader, quals, lineno, self.phred)
                if check == 1:
                    logging.error('Invalid header in fastq read ; record number : {0}'.format(lineno))
                    raise NameError('Invalid header in fastq read; record number : {0}'.format(lineno))
                if check == 2:
                    logging.error('Invalid secondary header in fastq read ; record number : {1}'.format(lineno))
                    raise NameError('Invalid secondary header in fastq read; record number : {0}'.format(lineno))
                if check == 3:
                    logging.error('Sequence and quality strings of unequl length in fastq read; record number : {0}'.format(lineno))
                    raise NameError('Sequence and quality strings of unequal length in fastq read; record number : {0}'.format(lineno))
                if check == 4:
                    logging.error('Sequence data missing; record number : {0}'.format(lineno))
                    raise NameError('Sequence data missing; record number : {0}'.format(lineno))

                #Optionally mask low quality reads
                #seq, quals = self.qualmasker(seq, quals)

                #Return record
                record = Record(header, sheader, seq, quals)
                yield(record)
            except StopIteration:
                logging.info('End of file')
                break
            except NameError:
                break

        return

class Metrics:

    def __init__(self, fastq):
        self.fastq = fastq

    def avgReadLen(self):
        fastq_reader = Fastq(self.fastq, './', 'phred33')
        total_length = 0
        total_reads = 0
        for lines in fastq_reader.read():
            total_length += len(lines.seq)
            total_reads += 1
            if total_reads >= 100:
                break

        avg_length = total_length/float(total_reads)
        return(avg_length)

class Prepper:

    def __init__(self, input_path):
        self.input_path = os.path.abspath(input_path)
        self.prep_logger = logging.getLogger('MaRS.Prepper')

    def getFastqPaths(self):
        filenames = list()
        for subdir, dirname, files in os.walk(self.input_path):
            for filename in files:
                if '.fastq' in filename or '.fastq.gz' in filename:
                    filepath = subdir + os.sep + filename
                    filenames.append(filepath)
        self.prep_logger.debug('Found {0} fastq files in {1}'.format(len(filenames), self.input_path))
        return(filenames)

    def getReadNumbers(self, file_name):
        reader = Fastq(file_name, None, None)
        read_number = 0
        for rec in reader.read():
            read_number += 1
        return(read_number)

    def prepInputs(self):
        files = self.getFastqPaths()
        experiment = dict()
        for fastq in files:
            reader = Fastq(fastq, './', 'phred33')
            Sample = namedtuple('Sample', ['sample', 'libname', 'library', 'files', 'prep', 'paired', 'numreads'])
            rec = next(reader.read())
            identifier = Identifier(rec)
            metric = Metrics(fastq)
            isIllOld =  identifier.isIlluminaOld()
            isIllNew =  identifier.isIlluminaNew()
            isSraOld = identifier.isSraOld()
            isSraNew = identifier.isSraNew()
            isPac = identifier.isPacbio()
            seqType = ''
            libType = ''
            sample_regex = re.compile('_?[r,R][1,2]|_?[l,L]00[1,2,3,4]')
            sample = sample_regex.split(os.path.basename(fastq))[0]
            if isIllOld:
                paired_regex = re.compile('@\w+-?\w+:\d+:\d+:\d+:\d+#\d')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'
            elif isIllNew:
                paired_regex = re.compile('@\w+-?\w+:\d+:\w+-?\w+:\d+:\d+:\d+:\d+\s')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'
            elif isSraOld:
                paired_regex = re.compile('@\w+\.?\w? \w+-\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'
            elif isSraNew:
                paired_regex = re.compile('@\w+\.?\w? \w+-\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
                lib = re.findall(paired_regex, rec.header)[0]
                paired = False
                seqType = 'Illumina'
                if metric.avgReadLen():
                    libType = 'Short'
            elif isPac:
                lib_regex = re.compile('@\w+_\d+_\d+_\w+')
                lib = re.findall(lib_regex, rec.header)[0]
                paired = False
                seqType = 'Pacbio'
                if metric.avgReadLen():
                    libType = 'Long'
            else:

                self.prep_logger.warning('Read from {0} with header : {1} does not follow any defined fastq header format.Please correct it'.format(fastq, rec_header))
            try:
                paired = True
                experiment[sample] = Sample(sample, lib, seqType, experiment[sample].files + [fastq], libType, paired, 0)
            except (KeyError, AttributeError):
                experiment[sample] = Sample(sample, lib, seqType, [fastq], libType, paired, 0)
        self.prep_logger.info('A total of {0} libraries were identified from the given folder {1}'.format(len(experiment), self.input_path))
        self.prep_logger.debug('The following libraries were detected in the given folder : {0}'.format(self.input_path))
        for sample, values in experiment.items():
            self.prep_logger.debug('Sample : {0}; Library: {1} ; Sequence type: {2} ; Files: {3} ; Library type: {4} ; Paired: {5} ; Total number of reads: {6}'.format(
                    values.sample, values.libname, values.library, ''.join(values.files), values.prep, values.paired, values.numreads))
        return(experiment)

class Identifier:

    def __init__(self, record):
        self.rec = record


    def isIlluminaOld(self):
        #@HWUSI-EAS100R:6:73:941:1973#0/1
        header_regex = re.compile('@\w+-?\w+:\d+:\d+:\d+:\d+#\d*')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)


    def isIlluminaNew(self):
        #@D00468:24:H8ELMADXX:1:1101:1470:2237 1:N:0:2
        header_regex = re.compile('@\w+-?\w+:\d+:\w+-?\w+:\d+:\d+:\d+:\d+\s\d:\w+:\w+:\w*')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)


    def isSraOld(self):
        #@SRR037455.1 HWI-E4_6_30ACL:4:1:0:29 length=35
        #@SRR902931.1 HWI-ST1384:61:D1DJ4ACXX:8:1101:1240:2015 length=50
        header_regex = re.compile('@\w+\.?\w? \w+-\w+:\d+:\d+:\d+:\d+ length=\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)

    def isSraNew(self):
        header_regex = re.compile('@\w+\.?\w? \w+-\w+:\d+:\w+:\d+:\d+:\d+:\d+ length=\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)

    def isPacbio(self):
        #@m160113_152755_42135_c100906712550000001823199104291667_s1_p0/15/7044_26271
        header_regex = re.compile('@\w+/\d+/\d+_\d+')
        match = re.fullmatch(header_regex, self.rec.header)
        if match != None:
            return(True)
        else:
            return(False)


if __name__ == '__main__':
    #Define deffault paths and aligner informations
    def_path = "{0}/lib".format(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
    bbduk_def = "{0}/bbmap/bbduk.sh".format(def_path)
    adp_def = "{0}/bbmap/resources/adapters.fa".format(def_path)
    tph_def = "{0}/tophat2/tophat2".format(def_path)
    cuf_def = "{0}/cufflinks/cufflinks".format(def_path)
    out_def = "{0}/local".format(os.path.abspath(os.path.dirname(os.path.realpath(__file__))))
    #Get arguments
    parser = argparse.ArgumentParser(prog='Nestv2')
    parser.add_argument('-i', '--inp_path', type=str,
                        help='Path to input directory')
    parser.add_argument('-r', '--ref', dest='ref_path', type=str,
                        help='Path to Reference bowtie index', required=True)
    parser.add_argument('-g', '--gtf', dest='gtf_path', type=str,
                        help='Path to Reference GTF file', required=True)
    parser.add_argument('-o', '--outpath', dest='out_path', type=str,
                        help='Path where all outputs will be stored', required=True)
    parser.add_argument('-a', '--adapter', dest='adp_path', type=str,
                        help='Path to Adpater fasta file', default=adp_def)
    parser.add_argument('-t', '--tophat', dest='tph_path', type=str,
                        help='Path to Tophat2 executable', default=tph_def)
    parser.add_argument('-c', '--cufflinks', dest='cuf_path', type=str,
                        help='Path to Cufflinks executable', default=cuf_def)
    parser.add_argument('-b', '--bbduk', dest='bbd_path', type=str,
                        help='Path to BBDuk executable', default=bbduk_def)
    args = parser.parse_args()

    #Validate parsed arguments
    if os.path.exists(os.path.abspath('{0}.fa'.format(args.ref_path))):
        ref_path = os.path.abspath(args.ref_path)
    else:
        logging.error('Path does not exists : {0}'.format(args.ref_path))
    if os.path.exists(os.path.abspath(args.gtf_path)):
        gtf_path = os.path.abspath(args.gtf_path)
    else:
        logging.error('Path does not exists : {0}'.format(args.gtf_path))
    if os.path.exists(os.path.abspath(args.out_path)):
        out_path = os.path.abspath(args.out_path)
    else:
        out_path = os.path.abspath(args.out_path)
        os.mkdir(out_path)
    if os.path.exists(os.path.abspath(args.adp_path)):
        adp_path = os.path.abspath(args.adp_path)
    else:
        logging.error('Path does not exists : {0}'.format(args.adp_path))
    if os.path.exists(os.path.abspath(args.tph_path)):
        tph_path = os.path.abspath(args.tph_path)
    else:
        logging.error('Path does not exists : {0}'.format(args.tph_path))
    if os.path.exists(os.path.abspath(args.cuf_path)):
        cuf_path = os.path.abspath(args.cuf_path)
    else:
        logging.error('Path does not exists : {0}'.format(args.cuf_path))
    if os.path.exists(os.path.abspath(args.bbd_path)):
        bbd_path = os.path.abspath(args.bbd_path)
    else:
        logging.error('Path does not exists : {0}'.format(args.bbd_path))

    #Prepare inputs
    prepper = Prepper(args.inp_path)
    experiment = prepper.prepInputs()

    #Iterate through experiments and run RNA Seq pipeline
    for study, details in experiment.items():
        rones = list()
        rtwos = list()
        sout_path = '{0}/{1}'.format(out_path, study)
        if not os.path.exists(sout_path):
            os.mkdir(sout_path)
        for files in details.files:
            if re.findall('.*_R1.*', files):
                rones.append(files)
            else:
                rtwos.append(files)
        #Preprocess reads
        crones = list()
        crtwos = list()
        cleaner = QualCheck(bbd_path, adp_path, sout_path)
        for rone, rtwo in zip(rones, rtwos):
            crone, crtwo, ret = cleaner.bbduk(rone, rtwo)
            if ret != 0 :
                logging.error('BBDuk failed to complete')
                sys.exit()
            crones.append(crone)
            crtwos.append(crtwo)
        #Align to Reference
        fastq_dict = {'R1' : crones, 'R2': crtwos}
        aligner = SplicedAligner(tph_path, ref_path, gtf_path, sout_path)
        bam_path, ret = aligner.tophat2(fastq_dict)
        #Run Cufflinks
        counter = GeneQuant(cuf_path, gtf_path, bam_path, sout_path)
        ret = counter.cuffLinks()
