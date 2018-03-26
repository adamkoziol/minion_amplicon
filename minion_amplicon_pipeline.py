#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import make_path, printtime
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
from argparse import ArgumentParser
from subprocess import call
from csv import DictReader
import multiprocessing
from time import time
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
        self.reformat_reads()
        self.make_blastdb(self.target, self.db)
        self.blast()
        self.blast_parser()
        self.populate_unique()
        self.find_unique()
        self.create_lists()
        self.bin_fastq()
        if self.assemble_reads:
            self.assemble()
        self.target_creation(self.target)
        self.bowtie_build(self.target)
        self.bowtie_run()
        self.samtools_index()
        self.extract_overhangs()
        # If the target file with wheat sequences exists, use this file for the subsequent analyses. Otherwise use
        # the original targets file
        target = self.wheat if os.path.isfile(self.wheat) else self.target
        db = self.wheat_db if os.path.isfile(self.wheat) else self.db
        self.make_blastdb(target, db)
        self.target_creation(target)
        self.bowtie_build(target)
        self.blast_overhangs(db)
        self.parse_overhang_blast()
        self.output_overhang_results()
        self.overhang_bowtie_run()

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
        bait_command = 'bbduk.sh -overwrite=true -in={reads} ' \
                       'threads={cpus} outm={filteredbaits} k={kmer} maskmiddle=t hdist={hdist} ref={reference}'\
            .format(reads=self.reads,
                    cpus=self.cpus,
                    kmer=self.kmer,
                    hdist=self.hdist,
                    reference=self.target,
                    filteredbaits=self.filteredfasta
                    )
        # Run the system call if the baited FASTQ file hasn't already been created
        if not os.path.isfile(self.filteredfasta):
            call(bait_command, shell=True, stdout=self.devnull, stderr=self.devnull)

    def reformat_reads(self):
        """
        Remove any reads below the user-inputted minimum read length
        """
        printtime('Filtering FASTQ reads below {} bp'.format(self.minreadlength), self.start)
        # Create the system call to reformat.sh - use the user-supplied minimum read length value
        filter_command = 'reformat.sh -ow=true -in={reads} -out={reads} -minlength={minlength}'\
            .format(reads=self.filteredfasta,
                    minlength=self.minreadlength)
        # Run the system call
        call(filter_command, shell=True, stdout=self.devnull, stderr=self.devnull)

    def make_blastdb(self, target, db):
        """
        Create the BLAST database using the combined target file
        :param target:
        :param db:
        """
        printtime('Creating BLAST databases', self.start)
        # Create a variable to store the name of one of the database files
        nhr = '{}.nhr'.format(db)
        if not os.path.isfile(nhr):
            makedb_command = 'makeblastdb -dbtype nucl -in {target} -out {reference}'\
                .format(target=target,
                        reference=db)
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

    def target_creation(self, fasta):
        """
        Split the multi-FASTA target file into FASTA files for each target in the file
        """
        # Create a generator of all the FASTA records
        records = SeqIO.parse(fasta, 'fasta')
        for record in records:
            # Extract the path information from the target file name/path variable
            target_base = os.path.dirname(fasta)
            # Create the name of the new target-specific FASTA file
            target_file = os.path.join(target_base, '{}.tfa'.format(record.id))
            # Populate the dictionary of subject: target file
            self.target_files[record.id] = target_file
            # Populate the dictionary of subject: record
            self.record_dict[record.id] = record
            if not os.path.isfile(target_file):
                with open(target_file, 'w') as target:
                    SeqIO.write(record, target, 'fasta')

    def bowtie_build(self, fasta_file):
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
        Find the clipped off ends from the reads, and either concatenate the reference sequence to the overhang and
        export, or skip the concatenation step and export
        """
        printtime('Extracting 5\' and 3\' overhangs', self.start)
        left_read_ids = dict()
        right_read_ids = dict()
        for subject, bamfile in self.bamdict.items():
            # Initialise dictionaries to hold lists
            self.overhang_fasta[subject] = list()
            self.overhang_fastq[subject] = list()
            self.chimeric_overhang_dict[subject] = list()
            self.leftrecords[subject] = list()
            self.leftfastqrecords[subject] = list()
            self.chimericleftrecords[subject] = list()
            self.rightrecords[subject] = list()
            self.rightfastqrecords[subject] = list()
            self.chimericrightrecords[subject] = list()
            left_read_ids[subject] = list()
            right_read_ids[subject] = list()
            # Extract the reference sequence from the dictionary - set it to lowercase for quick identification
            # in chimeric files
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
                        # Create a sequence consisting of the 5' overlap
                        seq = Seq(
                            record.query_sequence[0: int('{}'.format(record.cigartuples[0][1]))],
                            IUPAC.unambiguous_dna)
                        # Create a sequence record of the overhang sequence
                        seqrecord = SeqRecord(seq,
                                              id=record.qname,
                                              description=str())
                        fastq_seqrecord = SeqRecord(seq,
                                                    id=record.qname,
                                                    description=str(),
                                                    letter_annotations={
                                                        'phred_quality':
                                                            list(record.query_qualities)
                                                            [0: int('{}'.format(record.cigartuples[0][1]))]})
                        # Append the sequence record to the list
                        if record.qname not in left_read_ids[subject]:
                            self.leftrecords[subject].append(seqrecord)
                            self.leftfastqrecords[subject].append(fastq_seqrecord)
                            left_read_ids[subject].append(record.qname)

                        # Create a sequence consisting of the 5' overlap joined to the reference sequence
                        chimera_seq = Seq('{}{}'.format(
                            record.query_sequence[0: int('{}'.format(record.cigartuples[0][1]))],
                            refseq),
                            IUPAC.unambiguous_dna)

                        # Create a sequence record of the overhang sequence
                        chimera_seqrecord = SeqRecord(chimera_seq,
                                                      id=record.qname,
                                                      description=str())
                        # Append the sequence record to the list
                        self.chimericleftrecords[subject].append(chimera_seqrecord)
                # Examine the first entry in the final tuple in the cigartuples attribute. As above the entry must be 4
                if record.cigartuples[-1][0] == 4:
                    # Ensure long overhangs
                    if record.cigartuples[-1][1] >= self.overhang_length:
                        seq = Seq(
                            record.query_sequence[-int('{}'.format(record.cigartuples[-1][1])): -1],
                            IUPAC.unambiguous_dna)
                        seqrecord = SeqRecord(seq,
                                              id=record.qname,
                                              description=str())
                        fastq_seqrecord = SeqRecord(seq,
                                                    id=record.qname,
                                                    description=str(),
                                                    letter_annotations={
                                                        'phred_quality':
                                                            list(record.query_qualities)
                                                            [-int('{}'.format(record.cigartuples[-1][1])): -1]})
                        # Append the record to the list of records
                        # Append the sequence record to the list
                        if record.qname not in right_read_ids[subject]:
                            self.rightrecords[subject].append(seqrecord)
                            self.rightfastqrecords[subject].append(fastq_seqrecord)
                            right_read_ids[subject].append(record.qname)
                        chimera_seq = Seq('{}{}'.format(
                            refseq,
                            record.query_sequence[-int('{}'.format(record.cigartuples[-1][1])): -1]),
                            IUPAC.unambiguous_dna)
                        chimera_seqrecord = SeqRecord(chimera_seq,
                                                      id=record.qname,
                                                      description=str())
                        # Append the record to the list of records
                        self.chimericrightrecords[subject].append(chimera_seqrecord)

            # Ensure that records exist before attempting to create an output file
            if self.leftrecords[subject]:
                # Set the name of the file
                left_file = os.path.join(self.overhang_path, '{}_left_overhang'.format(subject))
                printtime('Saving {numreads} 5\' overhangs for {sub}'
                          .format(numreads=len(self.leftrecords[subject]),
                                  sub=subject), self.start)
                # Don't overwrite the output file
                fasta = left_file + '.fasta'
                fastq = left_file + '.fastq'
                # Populate the dictionary with subject: file name
                self.overhang_fasta[subject].append(fasta)
                self.overhang_fastq[subject].append(fastq)
                if not os.path.isfile(fasta):
                    # Write the records to file
                    with open(fasta, 'w') as left:
                        SeqIO.write(self.leftrecords[subject], left, 'fasta')
                if not os.path.isfile(fastq):
                    # Write the records to file
                    with open(fastq, 'w') as left:
                        SeqIO.write(self.leftfastqrecords[subject], left, 'fastq')

                # Set the name of the chimeric file
                chimeric_left_file = os.path.join(self.overhang_path,
                                                  '{}_chimeric_left_overhang.fasta'.format(subject))
                # Populate the dictionary with subject: file name
                self.chimeric_overhang_dict[subject].append(left_file)
                # Don't overwrite the output file
                if not os.path.isfile(chimeric_left_file):
                    # Write the records to file
                    with open(chimeric_left_file, 'w') as left:
                        SeqIO.write(self.chimericleftrecords[subject], left, 'fasta')
            # Same as above, but with the 3' overhangs
            if self.rightrecords[subject]:
                right_file = os.path.join(self.overhang_path, '{}_right_overhang'.format(subject))

                printtime('Saving {numreads} 3\' overhangs for {sub}'
                          .format(numreads=len(self.rightrecords[subject]),
                                  sub=subject), self.start)
                fasta = right_file + '.fasta'
                fastq = right_file + '.fastq'
                self.overhang_fasta[subject].append(fasta)
                self.overhang_fastq[subject].append(fastq)
                if not os.path.isfile(fasta):
                    # Write the records to file
                    with open(fasta, 'w') as right:
                        SeqIO.write(self.rightrecords[subject], right, 'fasta')
                if not os.path.isfile(fastq):
                    # Write the records to file
                    with open(fastq, 'w') as right:
                        SeqIO.write(self.rightfastqrecords[subject], right, 'fastq')

                chimeric_right_file = os.path.join(self.overhang_path,
                                                   '{}_chimeric_right_overhang.fasta'.format(subject))
                self.chimeric_overhang_dict[subject].append(right_file)
                if not os.path.isfile(chimeric_right_file):
                    with open(chimeric_right_file, 'w') as right:
                        SeqIO.write(self.chimericrightrecords[subject], right, 'fasta')

    def blast_overhangs(self, db):
        """
        BLAST overhangs against the original target file to determine the neighbours of each gene
        """
        printtime('Performing BLAST analyses', self.start)
        for subject, overhang_files in self.overhang_fasta.items():
            self.blast_report_dict[subject] = list()
            out_path = os.path.join(self.overhang_path, subject)
            make_path(out_path)
            for overhang_file in overhang_files:
                direction = 'left' if 'left' in overhang_file else 'right'
                report = os.path.join(out_path, '{sub}_{dir}_blast_results.csv'
                                      .format(sub=subject,
                                              path=out_path,
                                              dir=direction))
                self.blast_report_dict[subject].append(report)
                try:
                    size = os.path.getsize(report)
                    # If a report was created, but no results entered - program crashed, or no sequences passed
                    # thresholds, remove the report, and run the blast analyses again
                    if size == 0:
                        os.remove(report)
                except FileNotFoundError:
                    pass
                # BLAST command line call. Note the high number of alignments. Due to the fact that all the
                # targets are combined into one database, this is to ensure that all potential alignments are
                # reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
                blastn = NcbiblastnCommandline(query=overhang_file,
                                               db=db,
                                               num_alignments=1000000,
                                               num_threads=self.cpus,
                                               outfmt="'6 qseqid sseqid positive mismatch gaps "
                                                      "evalue bitscore slen length qstart qend sstart send'",
                                               out=report)
                # Only run blast if the report doesn't exist
                if not os.path.isfile(report):
                    blastn()

    def parse_overhang_blast(self):
        """
        Parse the BLAST reports, and populate dictionaries of the length of the match is greater than the user-supplied
        value for minimum match length
        """
        printtime('Parsing overhang BLAST outputs', self.start)
        for subject, blast_reports in self.blast_report_dict.items():
            self.overhang_blast_results[subject] = dict()
            for report in blast_reports:
                direction = 'left' if 'left' in report else 'right'
                self.overhang_blast_results[subject][direction] = dict()
                # Open the sequence profile file as a dictionary
                blastdict = DictReader(open(report), fieldnames=self.fieldnames, dialect='excel-tab')
                # Go through each BLAST result
                for row in blastdict:
                    # Create variables to store the subject and query names, as well as the length of the alignment
                    gene = row['subject_id']
                    query = row['query_id']
                    length = int(row['alignment_length'])
                    # Only save the name of the read if the hit is longer than the user-supplied value
                    if length >= self.length:
                        # print(subject, direction, gene, query, length)
                        try:
                            self.overhang_blast_results[subject][direction][gene].add(query)
                        except KeyError:
                            self.overhang_blast_results[subject][direction][gene] = set()
                            self.overhang_blast_results[subject][direction][gene].add(query)

    def output_overhang_results(self):
        """
        Parse the BLAST outputs to determine with which other genes in the cassette the overhang sequences share
        sequence identity
        """
        printtime('Binning overhang BLAST hits', self.start)
        fastq_bins = dict()
        fasta_bins = dict()
        # Start iterating through the deeply nested dictionary
        for subject, nested_dict in self.overhang_blast_results.items():
            fastq_bins[subject] = dict()
            fasta_bins[subject] = dict()
            for direction, nested_gene_dict in nested_dict.items():
                fastq_bins[subject][direction] = dict()
                fasta_bins[subject][direction] = dict()
                for gene, query_set in nested_gene_dict.items():
                    for file_type in ['fastq', 'fasta']:
                        if file_type == 'fastq':
                            self.overhang_fastq_bins = self.output_file(subject, direction, self.overhang_fastq,
                                                                        gene, query_set, file_type, fastq_bins)
                        else:
                            self.overhang_fasta_bins = self.output_file(subject, direction, self.overhang_fasta,
                                                                        gene, query_set, file_type, fasta_bins)

    def output_file(self, subject, direction, dictionary, gene, query_set, file_type, outdict):
        """

        """
        # Determine which overhang sequence-containing FASTA file to use depending on the direction.
        fastafile = dictionary[subject][0] if direction in dictionary[subject][0] \
            else dictionary[subject][1]
        # Load the records from the FASTA file to a dictionary
        record_dict = SeqIO.to_dict(SeqIO.parse(fastafile, file_type))
        # Set the name of the output path for the BLAST-binned reads
        output_path = os.path.join(self.bin_path, subject)
        make_path(output_path)
        # Create the name of the file. It must contain the subject, the matching gene, and the direction
        output_file = os.path.join(output_path, '{subject}_{gene}_{dir}.{filetype}'
                                   .format(subject=subject,
                                           gene=gene,
                                           dir=direction,
                                           filetype=file_type))
        outdict[subject][direction][gene] = output_file
        # Create a text file to store only the read names
        text_file = '{}.txt'.format(os.path.splitext(output_file)[0])
        # Write the read names to the file
        with open(text_file, 'w') as text:
            text.write('\n'.join(query_set))
        # Create a list to store all the reads for the individual subject:gene match:direction combinations
        read_list = list()
        for read in query_set:
            read_list.append(record_dict[read])
        # Use SeqIO to output the reads to the output file
        with open(output_file, 'w') as output:
            SeqIO.write(read_list, output, file_type)
        return outdict

    def overhang_bowtie_run(self):
        """
        Map the binned FASTQ reads against the appropriate target file
        """
        printtime('Performing reference mapping', self.start)
        for subject, nested_dict in self.overhang_fastq_bins.items():
            for direction, nested_gene_dict in nested_dict.items():
                for gene, fastq in nested_gene_dict.items():
                    # Set the name and create the folder to store the mapping outputs
                    outpath = os.path.join(self.overhang_path, subject)
                    sorted_bam = os.path.join(outpath, '{subject}_{gene}_{dir}_sorted.bam'
                                              .format(subject=subject,
                                                      gene=gene,
                                                      dir=direction))
                    # Populate the dictionary of subject: sorted BAM file name/path
                    self.bamdict[subject] = sorted_bam
                    # Bowtie2 command piped to samtools view to convert data to BAM format piped to samtools sort to
                    # sort the BAM file
                    map_command = 'bowtie2 -x {base} -U {fastq} --very-sensitive-local --local -p {threads} | ' \
                                  'samtools view -@ {threads} -h -F 4 -bT {target} - | ' \
                                  'samtools sort - -@ {threads} -o {sortedbam}'\
                        .format(base=self.basedict[gene],
                                fastq=fastq,
                                threads=self.cpus,
                                target=self.targetdict[gene],
                                sortedbam=sorted_bam)
                    # Perform reference mapping
                    if not os.path.isfile(sorted_bam):
                        call(map_command, shell=True, stdout=self.devnull, stderr=self.devnull)
                    index_command = 'samtools index {bamfile}'.format(bamfile=sorted_bam)
                    # Index mapped file
                    if not os.path.isfile('{bamfile}.bai'.format(bamfile=sorted_bam)):
                        call(index_command, shell=True, stdout=self.devnull, stderr=self.devnull)

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
        self.assemble_reads = args.assemble
        self.baitpath = os.path.join(self.path, 'bait')
        make_path(self.baitpath)
        self.blastpath = os.path.join(self.path, 'blast')
        make_path(self.blastpath)
        self.outpath = os.path.join(self.path, 'outputs')
        make_path(self.outpath)
        self.sequencepath = os.path.join(self.path, 'sequences')
        make_path(self.sequencepath)
        if self.assemble_reads:
            self.assemblypath = os.path.join(self.path, 'assemblies')
            make_path(self.assemblypath)
        self.overhang_path = os.path.join(self.path, 'overhangs')
        make_path(self.overhang_path)
        self.bin_path = os.path.join(self.path, 'overhang_bins')
        make_path(self.bin_path)
        self.target = os.path.join(args.targetfile)
        assert os.path.isfile(self.target), 'Cannot find the target file supplied in the arguments {}. Please ensure ' \
                                            'that the file name and path are correct'.format(self.target)
        self.wheat = os.path.join(os.path.dirname(self.target), 'genesofinterest_wheat.fasta')
        self.target_files = dict()
        # Remove the path and the file extension for the target file for use in BLASTing
        self.db = os.path.splitext(self.target)[0]
        self.wheat_db = os.path.splitext(self.wheat)[0]
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
        self.leftfastqrecords = dict()
        self.chimericleftrecords = dict()
        self.rightrecords = dict()
        self.rightfastqrecords = dict()
        self.chimericrightrecords = dict()
        self.overhang_fasta = dict()
        self.overhang_fastq = dict()
        self.overhang_fastq_bins = dict()
        self.overhang_fasta_bins = dict()
        self.chimeric_overhang_dict = dict()
        self.blast_report_dict = dict()
        self.overhang_blast_results = dict()
        self.overhang_lists = dict()
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
    parser.add_argument('-a', '--assemble',
                        action='store_true',
                        help='Optionally assemble binned reads')
    arguments = parser.parse_args()
    arguments.startingtime = time()
    # Run the pipeline
    pipeline = MinionPipeline(arguments)
    pipeline.main()

'''
/home/bioinfo/Bioinformatics/180323/gwalk2_update
-t
/home/bioinfo/Bioinformatics/targets/genesofinterest.fasta
-f
/home/bioinfo/Bioinformatics/data/gwalk2_update/gwalk_combined_trimmed.fastq.gz
'''