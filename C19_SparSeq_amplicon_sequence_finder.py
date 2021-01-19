import re
fileoutSrbd = open("R1R2_Srbd_top3_PCR_hits.txt", "w")
fileoutSrbd.write(">Refseq_Srbd_pst_strand" + "\n" + "ATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGT" + "\n")

fileoutSpoly = open("R1R2_Spbs_top3_PCR_hits.txt", "w")
fileoutSpoly.write(">Refseq_Spoly_pst_strand" + "\n" + "TATGCGCTAGTTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGCTAGTCAATCCATCATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTAC" + "\n")

fileoutRdRP = open("R1R2_RdRP_top3_PCR_hits.txt", "w")
fileoutRdRP.write(">Refseq_RdRP_pst_strand" + "\n" + "GATGCCACAACTGCTTATGCTAATAGTGTTTTTAACATTTGTCAAGCTGTCACGGCCAATGTTAATGCACTTTTATCTACTGATGGTAACAAAATTGCCGATAAGTATGTCCGCAA" + "\n")

with open('stitched_files_selected.txt') as f:
    for line in f:
        filename = line.strip()
        #Part-1 binning sequences and creating "sequence - readcount" info for reads longer than 90nt
        #Minimum required read length can be changed from 90 to desired value at line 23.
        long_read_dict = dict()
        it=2
        with open(filename) as f:
            for line in f:
                it=it+1
                if it%4==0 :
                    seq = line.strip()
                    if ">" not in seq and len(seq)>90 :
                        if seq not in long_read_dict : long_read_dict[seq] = 1
                        else : long_read_dict[seq] = 1+long_read_dict[seq]

        #Part-2 binning sequences with read counts into gene-specific primer pools, filtering out reads with adaptor sequence
        rdrp_ = list()
        Srbd_ = list()
        Spoly_ = list()
        rdrp_c_ = 0
        Srbd_c_ = 0
        Spoly_c_ = 0

        long_read_t = long_read_dict.items()
        for key,val in long_read_t :
            if "ACTGCTTATG" in key and "AGATCGGAAG" not in key :
                rdrp_.append((val,key))
                rdrp_c_ = rdrp_c_ + val
            elif "GGTAGCACACCT" in key and "AGATCGGAAG" not in key :
                Srbd_.append((val, key))
                Srbd_c_ = Srbd_c_+ val
            elif "TCAGACTCAGAC" in key and "AGATCGGAAG" not in key :
                Spoly_.append((val, key))
                Spoly_c_ = Spoly_c_+ val
        #Sort amplicon read count and write the top enriched amplicon sequences.
        #Desired number of top enriched amplicons can be reported by chenging the parameter as in example:
        # "rdrp_ltop = rdrp_[:1]" -> "rdrp_ltop = rdrp_[:3]" Second version outputs top3 enriched amplicon sequences.
        rdrp_.sort(reverse=True)
        rdrp_ltop = rdrp_[:1]
        for val,item in rdrp_ltop :
            if val>50:
                fileoutRdRP.write(">" + filename + "_RdRP_" + str(val) + "of" + str(rdrp_c_) + "\n" + item + "\n")
        # Run the two lines below if output of all detected amplicon read counts are desired.
        #for val,item in rdrp_ :
        #    fileout1.write(">" + filename[:18] + "_RdRP_" + str("{:.2f}".format(val/rdrp_c_*100)) + "%" + " " + item + "\n")

        Srbd_.sort(reverse=True)
        Srbd_ltop = Srbd_[:1]
        for val,item in Srbd_ltop :
            if val>50:
                fileoutSrbd.write(">" + filename + "_Srbd_" + str(val) + "of" + str(rdrp_c_) + "\n" + item + "\n")
        #for val,item in Srbd_ :
        #    fileout4.write(">" + filename[:18] + "_Srbd_" + str("{:.2f}".format(val/Srbd_c_*100)) + "%" + " " + item + "\n")

        Spoly_.sort(reverse=True)
        Spoly_ltop = Spoly_[:1]
        for val,item in Spoly_ltop :
            if val>50:
                fileoutSpoly.write(">" + filename + "_Spbs_" + str(val) + "of" + str(rdrp_c_) + "\n" + item + "\n")
        #for val,item in Spoly_ :
        #    fileout8.write(">" + filename[:18] + "_Spbs_" + str("{:.2f}".format(val/Spoly_c_*100)) + "%" + " " + item + "\n")
