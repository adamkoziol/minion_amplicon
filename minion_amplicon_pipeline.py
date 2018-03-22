#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import make_path, printtime
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Application import ApplicationError
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
from argparse import ArgumentParser
from subprocess import call
from csv import DictReader
import multiprocessing
from time import time
import shutil
import pysam
import os

__author__ = 'adamkoziol'


class MinionPipeline(object):

    def main(self):
        """
        Run the methods in the correct order
        """
        if self.trim:
            self.combine_fastq()
            self.trim_fastq()
            self.clump_reads()
        self.bait()
        self.fastq_to_fasta()
        self.make_blastdb()
        self.blast()
        self.blast_parser()
        self.populate_unique()
        self.find_unique()
        self.create_lists()
        self.bin_fastq()
        self.assemble()
        self.target_creation()
        self.bowtie_build()
        self.bowtie_run()
        self.samtools_index()
        self.extract_overhangs()
        self.overhang_aligner()

    def combine_fastq(self):
        """
        Read in all the FASTQ files, and create a compressed output file
        """
        printtime('Combining and compressing FASTQ files', self.start)
        combine_command = 'cat {path}*.fastq | gzip > {combinedfastq} '\
            .format(path=os.path.join(self.fastq_path, ''),
                    combinedfastq=self.combined_reads)
        if not os.path.isfile(self.combined_reads):
            call(combine_command, shell=True, stdout=self.devnull, stderr=self.devnull)

    def trim_fastq(self):
        """
        Use porechop to trim the combined FASTQ file
        """
        printtime('Trimming FASTQ reads', self.start)
        trim_command = 'porechop --check_reads 1000 -i {combinedreads} -o {fastqpath} --threads {threads}'\
            .format(combinedreads=self.combined_reads,
                    fastqpath=os.path.join(self.trimmed_reads),
                    threads=self.cpus)
        if not os.path.isfile(self.trimmed_reads):
            # Let the stdout and stderr be printed to the terminal
            call(trim_command, shell=True)

    def clump_reads(self):
        """
        Use clumpify to group reads together for more efficient compression
        """
        printtime('Clumping FASTQ reads', self.start)
        clump_command = 'clumpify.sh in={trimmedreads} out={clumpedreads} reorder ziplevel=9'\
            .format(trimmedreads=self.trimmed_reads,
                    clumpedreads=self.reads)
        if not os.path.isfile(self.reads):
            call(clump_command, shell=True, stdout=self.devnull, stderr=self.devnull)

    def bait(self):
        """
        Bait the FASTQ reads with the combined target file
        """
        printtime('Baiting FASTQ reads with targets', self.start)
        # Create the system call to bbduk - use the user-supplied kmer length, and hdist values
        bait_command = 'bbduk.sh -overwrite=true -in={reads} minlength={minreadlength} ' \
                       'threads={cpus} outm={filteredbaits} k={kmer} maskmiddle=t hdist={hdist} ref={reference}'\
            .format(reads=self.reads,
                    minreadlength=self.minreadlength,
                    cpus=self.cpus,
                    kmer=self.kmer,
                    hdist=self.hdist,
                    reference=self.no_wheat if self.no_wheat else self.target,
                    filteredbaits=self.filteredfastq
                    )
        # Run the system call if the baited FASTQ file hasn't already been created
        if not os.path.isfile(self.filteredfastq):
            call(bait_command, shell=True, stdout=self.devnull, stderr=self.devnull)

    def fastq_to_fasta(self):
        """
        Convert baited FASTQ reads to FASTA format for future BLAST analyses
        """
        printtime('Converting FASTQ to FASTA', self.start)
        convert_command = 'fastq_to_fasta -i {filteredbaits} -o {fastareads}'\
            .format(filteredbaits=self.filteredfastq,
                    fastareads=self.filteredfasta)
        if not os.path.isfile(self.filteredfasta):
            call(convert_command, shell=True, stdout=self.devnull, stderr=self.devnull)

    def make_blastdb(self):
        """
        Create the BLAST database using the combined target file
        """
        printtime('Creating BLAST databases', self.start)
        # Create a variable to store the name of one of the database files
        nhr = '{}.nhr'.format(self.db)
        if not os.path.isfile(nhr):
            makedb_command = 'makeblastdb -dbtype nucl -in {target} -out {reference}'\
                .format(target=self.target,
                        reference=self.db)
            call(makedb_command, shell=True, stdout=self.devnull, stderr=self.devnull)

    def blast(self):
        """
        Run BLAST analyses between the FASTA-formatted reads and the target database
        """
        printtime('Performing BLAST analyses', self.start)
        try:
            size = os.path.getsize(self.report)
            # If a report was created, but no results entered - program crashed, or no sequences passed thresholds,
            # remove the report, and run the blast analyses again
            if size == 0:
                os.remove(self.report)
        except FileNotFoundError:
            pass
        # BLAST command line call. Note the high number of alignments.
        # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
        # alignments are reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
        blastn = NcbiblastnCommandline(query=self.filteredfasta,
                                       db=self.db,
                                       num_alignments=1000000,
                                       num_threads=self.cpus,
                                       outfmt="'6 qseqid sseqid positive mismatch gaps "
                                              "evalue bitscore slen length qstart qend sstart send'",
                                       out=self.report)
        # Only run blast if the report doesn't exist
        if not os.path.isfile(self.report):
            blastn()

    def blast_parser(self):
        """
        Parse the BLAST reports, and populate dictionaries of the length of the match is greater than the user-supplied
        value for minimum match length
        """
        printtime('Parsing BLAST outputs', self.start)
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(open(self.report), fieldnames=self.fieldnames, dialect='excel-tab')
        # Go through each BLAST result
        for row in blastdict:
            # Create variables to store the subject and query names, as well as the length of the alignment
            subject = row['subject_id']
            query = row['query_id']
            length = int(row['alignment_length'])
            # Only save the name of the read if the hit is longer than the user-supplied value
            if length >= self.length:
                try:
                    self.resultdict[subject].add(query)
                except KeyError:
                    self.resultdict[subject] = set()
                    self.resultdict[subject].add(query)

    def populate_unique(self):
        """
        Create a dictionary of the number of times each read is present over all the targets
        """
        printtime('Determining read frequency', self.start)
        for subject, queryset in self.resultdict.items():
            for query in queryset:
                # Find the number of times a read is present in the complete analysis
                try:
                    self.read_dict[query] += 1
                except KeyError:
                    self.read_dict[query] = 1

    def find_unique(self):
        """
        Parse dictionaries to determine which reads correspond to a single target
        """
        printtime('Finding reads unique to each subject', self.start)
        for subject, queryset in self.resultdict.items():
            # Initialise the list to store
            self.unique_dict[subject] = list()
            for query in queryset:
                # If a read is only present once, add it to the dictionary of unique reads
                if self.read_dict[query] == 1:
                    self.unique_dict[subject].append(query)

    def create_lists(self):
        """
        Create a text file of the read names for each target, and a file containing only the unique reads
        """
        printtime('Outputting lists of reads', self.start)
        for subject, queryset in self.resultdict.items():
            output = os.path.join(self.outpath, '{}.txt'.format(subject))
            self.listdict[subject] = output
            # Output the names of all the reads
            with open(output, 'w') as read_list:
                read_list.write('\n'.join([query for query in queryset]))
            # Output the names of the unique reads
            unique_output = os.path.join(self.outpath, '{}_unique_reads.txt'.format(subject))
            with open(unique_output, 'w') as unique:
                unique.write('\n'.join(self.unique_dict[subject]))

    def bin_fastq(self):
        """
        Use filterbyname.sh to extract the binned reads from the original FASTQ file
        """
        printtime('Extracting binned reads', self.start)
        for subject, report in sorted(self.listdict.items()):
            fastq_file = os.path.join(self.sequencepath, '{}.fastq.gz'.format(subject))
            self.sequencedict[subject] = fastq_file
            # Note that only reads longer than the user-supplied minimum read length parameter
            bin_command = 'filterbyname.sh -in={reads}' \
                          ' -out={out} names={names} threads={cpu} include=true ow=t minlen={minlen}'\
                .format(reads=self.reads,
                        out=fastq_file,
                        names=report,
                        cpu=self.cpus,
                        minlen=self.minreadlength)
            if not os.path.isfile(fastq_file):
                printtime('Binning {numreads} ({unique} unique) reads for {subject}'
                          .format(numreads=len(self.resultdict[subject]),
                                  unique=len(self.unique_dict[subject]),
                                  subject=subject), self.start)
                call(bin_command, shell=True, stdout=self.devnull, stderr=self.devnull)

    def assemble(self):
        """
        Use canu to assemble the binned FASTQ reads
        """
        printtime('Assembling binned FASTQ reads', self.start)
        for subject, sequence in sorted(self.sequencedict.items()):
            outdir = os.path.join(self.assemblypath, subject)
            contigs_file = os.path.join(outdir, '{sub}.contigs.fasta'
                                        .format(sub=subject))
            canu_command = 'canu -p {alias} -d {outdir} genomeSize=15000 -nanopore-raw {binnedreads} ' \
                           'minReadLength={minlength} stopOnReadQuality=false'\
                .format(alias=subject,
                        outdir=outdir,
                        minlength=self.minreadlength,
                        binnedreads=sequence)
            if not os.path.isfile(contigs_file):
                printtime('Assembling {numreads} ({unique} unique) reads for {subject}'
                          .format(numreads=len(self.resultdict[subject]),
                                  unique=len(self.unique_dict[subject]),
                                  subject=subject), self.start)
                call(canu_command, shell=True, stdout=self.devnull, stderr=self.devnull)

    def target_creation(self):
        """
        Split the multi-FASTA target file into FASTA files for each target in the file
        """
        # Create a generator of all the FASTA records
        records = SeqIO.parse(self.target, 'fasta')
        for record in records:
            # Extract the path information from the target file name/path variable
            target_base = os.path.dirname(self.target)
            # Create the name of the new target-specific FASTA file
            target_file = os.path.join(target_base, '{}.tfa'.format(record.id))
            # Populate the dictionary of subject: target file
            self.target_files[record.id] = target_file
            # Populate the dictionary of subject: record
            self.record_dict[record.id] = record
            if not os.path.isfile(target_file):
                with open(target_file, 'w') as target:
                    SeqIO.write(record, target, 'fasta')

    def bowtie_build(self):
        """
        Use bowtie2-build to index each target file
        """
        printtime('Preparing targets for reference mapping', self.start)
        for subject, target_file in sorted(self.target_files.items()):
            base_name = os.path.splitext(target_file)[0]
            self.basedict[subject] = base_name
            self.targetdict[subject] = target_file
            build_command = 'bowtie2-build {target} {base}'.format(target=target_file,
                                                                   base=base_name)
            if not os.path.isfile('{}.1.bt2'.format(base_name)):
                call(build_command, shell=True, stdout=self.devnull, stderr=self.devnull)

    def bowtie_run(self):
        """
        Map the binned FASTQ reads against the appropriate target file
        """
        printtime('Performing reference mapping', self.start)
        for subject, fastq in sorted(self.sequencedict.items()):
            # Set the name and create the folder to store the mapping outputs
            outpath = os.path.join(self.sequencepath, subject)
            make_path(outpath)
            sorted_bam = os.path.join(outpath, '{}_sorted.bam'.format(subject))
            # Populate the dictionary of subject: sorted BAM file name/path
            self.bamdict[subject] = sorted_bam
            # Bowtie2 command piped to samtools view to convert data to BAM format piped to samtools sort to
            # sort the BAM file
            map_command = 'bowtie2 -x {base} -U {fastq} --very-sensitive-local --local -p {threads} | ' \
                          'samtools view -@ {threads} -h -F 4 -bT {target} - | ' \
                          'samtools sort - -@ {threads} -o {sortedbam}'\
                .format(base=self.basedict[subject],
                        fastq=fastq,
                        threads=self.cpus,
                        target=self.targetdict[subject],
                        sortedbam=sorted_bam)
            if not os.path.isfile(sorted_bam):
                printtime('Mapping {numreads} ({unique} unique) reads to {subject}'
                          .format(numreads=len(self.resultdict[subject]),
                                  unique=len(self.unique_dict[subject]),
                                  subject=subject), self.start)
                call(map_command, shell=True, stdout=self.devnull, stderr=self.devnull)

    def samtools_index(self):
        """
        Use samtools index to index the sorted BAM files
        """
        for subject, bamfile in self.bamdict.items():
            index_command = 'samtools index {bamfile}'.format(bamfile=bamfile)
            if not os.path.isfile('{bamfile}.bai'.format(bamfile=bamfile)):
                call(index_command, shell=True, stdout=self.devnull, stderr=self.devnull)

    def extract_overhangs(self):
        """
        Find the clipped off ends from the reads, concatenate the reference sequence to the overhang
        """
        printtime('Extracting 5\' and 3\' overhangs', self.start)
        for subject, bamfile in self.bamdict.items():
            # Initialise dictionaries to hold lists
            self.overhang_dict[subject] = list()
            self.leftrecords[subject] = list()
            self.rightrecords[subject] = list()
            # Extract the reference sequence from the dictionary
            refseq = str(self.record_dict[subject].seq.lower())
            # Create a pysam alignment file object of the sorted bam file
            bamfile = pysam.AlignmentFile(bamfile, 'rb')
            # Iterate through each read in the bamfile
            for record in bamfile.fetch():
                # Examine the first entry in the first tuple in the cigartuples attribute. If the entry is 4, then
                # the tuple refers to a soft-clip at the beginning of the reference sequence
                if record.cigartuples[0][0] == 4:
                    # Ensure that the length of the extracted overhang is worth considering
                    if record.cigartuples[0][1] >= self.overhang_length:
                        # Create a sequence consisting of the 5' overlap joined to the reference sequence
                        seq = Seq('{}{}'.format(
                            record.query_sequence[0: int('{}'.format(record.cigartuples[0][1]))],
                            refseq),
                            IUPAC.unambiguous_dna)
                        # Create a sequence record of the overhang sequence
                        seqrecord = SeqRecord(seq,
                                              id=record.qname,
                                              description=str())
                        # Append the sequence record to the list
                        self.leftrecords[subject].append(seqrecord)
                # Examine the first entry in the final tuple in the cigartuples attribute. As above the entry must be 4
                if record.cigartuples[-1][0] == 4:
                    # Ensure long overhangs
                    if record.cigartuples[0][1] >= self.overhang_length:
                        seq = Seq('{}{}'.format(
                            refseq,
                            record.query_sequence[-int('{}'.format(record.cigartuples[-1][1])): -1]),
                            IUPAC.unambiguous_dna)
                        seqrecord = SeqRecord(seq,
                                              id=record.qname,
                                              description=str())
                        # Append the record to the list of records
                        self.rightrecords[subject].append(seqrecord)
            # Ensure that records exist before attempting to create an output file
            if self.leftrecords[subject]:
                # Set the name of the file
                left_file = os.path.join(self.overhang_path, '{}_left_overhang.fasta'.format(subject))
                # Populate the dictionary with subject: file name
                self.overhang_dict[subject].append(left_file)
                printtime('Saving {numreads} 5\' overhangs for {sub}'
                          .format(numreads=len(self.leftrecords[subject]),
                                  sub=subject), self.start)
                # Don't overwrite the output file
                if not os.path.isfile(left_file):
                    # Write the records to file
                    with open(left_file, 'w') as left:
                        SeqIO.write(self.leftrecords[subject], left, 'fasta')
            # Same as above, but with the 3' overhangs
            if self.rightrecords[subject]:
                right_file = os.path.join(self.overhang_path, '{}_right_overhang.fasta'.format(subject))
                self.overhang_dict[subject].append(right_file)
                printtime('Saving {numreads} 3\' overhangs for {sub}'
                          .format(numreads=len(self.rightrecords[subject]),
                                  sub=subject), self.start)
                if not os.path.isfile(right_file):
                    with open(right_file, 'w') as right:
                        SeqIO.write(self.rightrecords[subject], right, 'fasta')

    def overhang_aligner(self):
        """
        Perform a multiple sequence alignment of the overhang sequences
        """
        from Bio.Align.Applications import ClustalOmegaCommandline
        printtime('Aligning overhangs', self.start)
        for subject, file_list in self.overhang_dict.items():
            for outputfile in file_list:
                aligned = os.path.join(self.aligned_overhangs, '{}_aligned.fasta'
                                       .format(os.path.splitext(os.path.basename(outputfile))[0]))
                # Create the command line call
                clustalomega = ClustalOmegaCommandline(infile=outputfile,
                                                       outfile=aligned,
                                                       threads=self.cpus,
                                                       auto=True)
                if not os.path.isfile(aligned):
                    printtime('Aligning overhangs for {}'.format(subject), self.start)
                    # Perform the alignments
                    try:
                        clustalomega()
                    # Files with a single sequence cannot be aligned. Copy the original file over to the aligned folder
                    except ApplicationError:
                        shutil.copyfile(outputfile, aligned)

    def __init__(self, args):
        """
        Initialises the variables required for this class
        :param args: list of arguments passed to the script
        """
        # Define variables from the arguments - there may be a more streamlined way to do this
        self.args = args
        self.path = os.path.join(args.path)
        make_path(self.path)
        self.length = args.length
        self.start = args.startingtime
        if os.path.isdir(args.fastq):
            self.trim = True
            self.fastq_path = os.path.join(args.fastq)
            self.combined_reads = os.path.join(self.fastq_path, 'combined.fastq.gz')
            self.trimmed_reads = os.path.join(self.fastq_path, 'trimmed.fastq.gz')
            self.reads = os.path.join(self.fastq_path, 'trimmed_clumped.fastq.gz')
        elif os.path.isfile(args.fastq):
            self.trim = False
            self.reads = os.path.join(args.fastq)
            assert os.path.isfile(self.reads), 'Cannot find the FASTQ file supplied in the arguments: {}. ' \
                                               'Please ensure that the file name and path are correct'\
                .format(self.reads)
        self.kmer = args.kmersize
        self.hdist = args.hdist
        self.minreadlength = args.readlength
        self.overhang_length = args.overhang
        self.baitpath = os.path.join(self.path, 'bait')
        make_path(self.baitpath)
        self.blastpath = os.path.join(self.path, 'blast')
        make_path(self.blastpath)
        self.outpath = os.path.join(self.path, 'outputs')
        make_path(self.outpath)
        self.sequencepath = os.path.join(self.path, 'sequences')
        make_path(self.sequencepath)
        self.assemblypath = os.path.join(self.path, 'assemblies')
        make_path(self.assemblypath)
        self.overhang_path = os.path.join(self.path, 'overhangs')
        make_path(self.overhang_path)
        self.aligned_overhangs = os.path.join(self.path, 'aligned_overhangs')
        make_path(self.aligned_overhangs)
        self.target = os.path.join(args.targetfile)
        assert os.path.isfile(self.target), 'Cannot find the target file supplied in the arguments {}. Please ensure ' \
                                            'that the file name and path are correct'.format(self.target)
        if os.path.isfile('{}_no_wheat.fasta'.format(os.path.splitext(self.target)[0])):
            self.no_wheat = '{}_no_wheat.fasta'.format(os.path.splitext(self.target)[0])
        else:
            self.no_wheat = str()
        self.target_files = dict()
        # Remove the path and the file extension for the target file for use in BLASTing
        self.db = os.path.splitext(self.target)[0]
        self.filteredfastq = os.path.join(self.baitpath, 'filteredfastq.fastq')
        self.filteredfasta = os.path.join(self.baitpath, 'filteredfastq.fasta')
        self.report = os.path.join(self.blastpath, 'blast_report.csv')
        self.resultdict = dict()
        self.read_dict = dict()
        self.unique_dict = dict()
        self.listdict = dict()
        self.record_dict = dict()
        self.sequencedict = dict()
        self.targetdict = dict()
        self.basedict = dict()
        self.bamdict = dict()
        self.leftrecords = dict()
        self.rightrecords = dict()
        self.overhang_dict = dict()
        self.devnull = open(os.devnull, 'wb')
        self.cpus = multiprocessing.cpu_count()
        # Fields used for custom outfmt 6 BLAST output:
        self.fieldnames = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps',
                           'evalue', 'bit_score', 'subject_length', 'alignment_length',
                           'query_start', 'query_end', 'subject_start', 'subject_end']

# If the script is called from the command line, then call the argument parser
if __name__ == '__main__':
    # Parser for arguments
    parser = ArgumentParser(description='Pipeline for transgene analyses')
    parser.add_argument('path',
                        help='Specify path')
    parser.add_argument('-l', '--length',
                        default=50,
                        help='Minimum length of match. Default is 50')
    parser.add_argument('-t', '--targetfile',
                        required=True,
                        help='Name and path of multi-FASTA file of targets')
    parser.add_argument('-k', '--kmersize',
                        default=27,
                        help='Size of kmer to use when baiting reads against targets')
    parser.add_argument('--hdist',
                        default=2,
                        help='Number of mismatches to allow for baiting. Default is 2.')
    parser.add_argument('-r', '--readlength',
                        default=500,
                        help='Minimum read length for canu')
    parser.add_argument('-f', '--fastq',
                        help='Supply either the mame and path of FASTQ file to process, or the path of raw FASTQ files '
                             'to combine and trim prior to the rest of the analyses')
    parser.add_argument('-o', '--overhang',
                        default=100,
                        help='The minimum length a 5\' or 3\' overhang must be to be considered')
    arguments = parser.parse_args()
    arguments.startingtime = time()
    # Run the pipeline
    pipeline = MinionPipeline(arguments)
    pipeline.main()
