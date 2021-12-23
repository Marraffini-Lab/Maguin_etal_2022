from Bio import SeqIO
import operator
import mmap


def spacer_finder(deep_seq_file, undetermined_deep_seq_file, name_output_file, seq_primer , genome_for, genome_revcomp):
    #open MiSeq file for a specific sample
    handle1 = open(deep_seq_file, "r")
    recs1 = SeqIO.parse(handle1, 'fastq')

    #open file with the reads for which the index was not resolved
    handle2 = open(undetermined_deep_seq_file, "r")
    recs2 = SeqIO.parse(handle2, 'fastq')

    nophage_seq = open(name_output_file, 'w')

    # open the phage genome txt file needed to match the spacer at the end of the script
    forward = open(genome_for)
    reverse = open(genome_revcomp)

    spacers_dictionary = {}
    reads = 0
    # extract spacers from determined index file
    for r in recs1:
        reads = reads + 1
        if reads % 100000 == 0:
            print reads
        repeat1 = 'GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC'
        repeat2 = 'GTTTTAGAGCTATGCTGTTTT'

        sequence = r.seq

        # skip reads longer than 146 bp as the single adapted reads should be under 146 bp in pur PCR
        if len(sequence) > 146:
            pass

        else:

            # nophage experiment
            # this search for the forward barcoded primer used seq_primer
            posF = sequence.find(seq_primer)

            if (posF != -1):
                rep1 = sequence.find(repeat1)

                if (rep1 != -1):
                    sequence = sequence[rep1 + 36:]

                    rep2 = sequence.find(repeat2)
                    if rep2 != -1:
                        spacer = sequence[rep2 - 20:rep2]

                        if spacer in spacers_dictionary:
                            spacers_dictionary[spacer] += 1
                        else:
                            spacers_dictionary[spacer] = 1
            else:
                rev_comp_seq = sequence.reverse_complement()

                # search for the barcoded primer in the reverse complement sequence
                posR = rev_comp_seq.find(seq_primer)
                if (posR != -1):
                    rep1 = rev_comp_seq.find(repeat1)

                    if (rep1 != -1):
                        rev_comp_seq = rev_comp_seq[rep1 + 36:]

                        rep2 = rev_comp_seq.find(repeat2)
                        if rep2 != -1:
                            spacer = rev_comp_seq[rep2 - 20:rep2]

                            if spacer in spacers_dictionary:
                                spacers_dictionary[spacer] += 1
                            else:
                                spacers_dictionary[spacer] = 1

    # extract spacers from undetermined indexes file
    reads_undeter = 0
    for r in recs2:
        reads_undeter = reads_undeter + 1
        if reads_undeter % 100000 == 0:
            print "undet" + " " + str(reads_undeter)

        repeat1 = 'GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC'
        repeat2 = 'GTTTTAGAGCTATGCTGTTTT'

        sequence = r.seq

        # nophage experiment
        # this search for the forward primer PM169
        # skip reads longer than 146 bp
        if len(sequence) > 146:
            pass

        else:
            posF = sequence.find(seq_primer)

            if (posF != -1):
                rep1 = sequence.find(repeat1)

                if (rep1 != -1):
                    sequence = sequence[rep1 + 36:]

                    rep2 = sequence.find(repeat2)
                    if rep2 != -1:
                        spacer = sequence[rep2 - 20:rep2]

                        if spacer in spacers_dictionary:
                            spacers_dictionary[spacer] += 1
                        else:
                            spacers_dicactionary[spacer] = 1
            else:
                rev_comp_seq = sequence.reverse_complement()

                # search for the barcoded primer in the reverse complement sequence
                posR = rev_comp_seq.find(seq_primer)
                if (posR != -1):
                    rep1 = rev_comp_seq.find(repeat1)

                    if (rep1 != -1):
                        rev_comp_seq = rev_comp_seq[rep1 + 36:]

                        rep2 = rev_comp_seq.find(repeat2)
                        if rep2 != -1:
                            spacer = rev_comp_seq[rep2 - 20:rep2]

                            if spacer in spacers_dictionary:
                                spacers_dictionary[spacer] += 1
                            else:
                                spacers_dictionary[spacer] = 1

    spacers_dictionary = sorted(spacers_dictionary.items(), key=operator.itemgetter(1), reverse=True)

    # creates a file with spacer, reads, location
    spacer_count = 0
    for spacer in spacers_dictionary:

        f = mmap.mmap(forward.fileno(), 0, access=mmap.ACCESS_READ)
        r = mmap.mmap(reverse.fileno(), 0, access=mmap.ACCESS_READ)

        # search in the forward orientation
        location = f.find(str(spacer[0]))
        if (location != -1):
            location = location + 21
            nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(location) + "\t" + "top_phage" + "\t" + "\n")


        else:
            # search in reverse orientation
            location = r.find(str(spacer[0]))
            if (location != -1):
                location_reverse = (len(r) - location - 20)
                nophage_seq.write(
                    str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + str(location_reverse) + "\t" + "bottom_phage" + "\t" + "\n")

            else:
                nophage_seq.write(str(spacer[0]) + "\t" + str(spacer[1]) + "\t" + '0' + "\t" + "not found" + "\t" + "\n")

    nophage_seq.close()
