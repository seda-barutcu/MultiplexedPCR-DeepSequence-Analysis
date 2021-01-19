fileout_R1_dimers = open("pilot_R1_dimers_per_primer.txt", "w")

RdRP_R1_Amplicon = "ACCGTTTCTATAGATTAGCTAATGAGTGTGCTCAAGTATTGAGTGAAATGGTCATGTGTGGCGGTTCACTATATGTTAAACCAGGTGGAACCTCATCAGGAGATGCCACAACTGCTTATGCTAATAGTGTTTTTAACATTTGTCAAGCTGTCACGGCCAATGTTAATGCACTT"
RdRP_R2_Amplicon = "AAGTGCATTAACATTGGCCGTGACAGCTTGACAAATGTTAAAAACACTATTAGCATAAGCAGTTGTGGCATCTCCTGATGAGGTTCCACCTGGTTTAACATATAGTGAACCGCCACACATGACCATTTCACTCAATACTTGAGCACACTCATTAGCTAATCTATAGAAACGGT"

E_R1_Amplicon = "AGACAGGTACGTTAATAGTTAATAGCGTACTTCTTTTTCTTGCTTTCGTGGTATTCTTGCTAGTTACACTAGCCATCCTTACTGCGCTTCGATTGTGTGCGTACTGCTGCAATATTGTTAACGTG"
E_R2_Amplicon = "CACGTTAACAATATTGCAGCAGTACGCACACAATCGAAGCGCAGTAAGGATGGCTAGTGTAACTAGCAAGAATACCACGAAAGCAAGAAAAAGAAGTACGCTATTAACTATTAACGTACCTGTCT"

N_R1_Amplicon = "CCAGGCAGCAGTAGGGGAACTTCTCCTGCTAGAATGGCTGGCAATGGCGGTGATGCTGCTCTTGCTTTGCTGCTGCTTGACAGATTGAACCAGCTTGAGAGCAAAATGTCTGGTAAAGGCCAA"
N_R2_Amplicon = "TTGGCCTTTACCAGACATTTTGCTCTCAAGCTGGTTCAATCTGTCAAGCAGCAGCAAAGCAAGAGCAGCATCACCGCCATTGCCAGCCATTCTAGCAGGAGAAGTTCCCCTACTGCTGCCTGG"

Srbd_R1_Amplicon = "ATCAGGCCGGTAGCACACCTTGTAATGGTGTTGAAGGTTTTAATTGTTACTTTCCTTTACAATCATATGGTTTCCAACCCACTAATGGTGTTGGTTACCAACCATACAGAGT"
Srbd_R2_Amplicon = "ACTCTGTATGGTTGGTAACCAACACCATTAGTGGGTTGGAAACCATATGATTGTAAAGGAAAGTAACAATTAAAACCTTCAACACCATTACAAGGTGTGCTACCGGCCTGAT"

Spoly_R1_Amplicon = "TATGCGCTAGTTATCAGACTCAGACTAATTCTCCTCGGCGGGCACGTAGTGTAGCTAGTCAATCCATCATTGCCTACACTATGTCACTTGGTGCAGAAAATTCAGTTGCTTAC"
Spoly_R2_Amplicon = "GTAAGCAACTGAATTTTCTGCACCAAGTGACATAGTGTAGGCAATGATGGATTGACTAGCTACACTACGTGCCCGCCGAGGAGAATTAGTCTGAGTCTGATAACTAGCGCATA"

PPIB_R1_Amplicon = "TGCTGCTGCCGGGACCTTCTGCGGCCGATGAGAAGAAGAAGGGGCCCAAAGTCACCGTCAAGGTGTATTTTGACCTACGAATTGGAGATGAAGATGTAGGCCGGGTGATCTTTGGTCTCTTCGGAA"
PPIB_R2_Amplicon = "TTCCGAAGAGACCAAAGATCACCCGGCCTACATCTTCATCTCCAATTCGTAGGTCAAAATACACCTTGACGGTGACTTTGGGCCCCTTCTTCTTCTCATCGGCCGCAGAAGGTCCCGGCAGCAGCA"

PPIB_F = "TGCTGCTGCCGGGACCTTC"
Srbd_F = "ATCAGGCCGGTAGCACACCT"
Spoly_F = "TATGCGCTAGTTATCAGACTCAGAC"
N_F = "CCAGGCAGCAGTAGGGGAAC,CCAGGCAGCAATAGGGGAAC"
E_F = "AGACAGGTACGTTAATAGTTAATAGCG"
RdRp_F = "ACCGTTTCTATAGATTAGCTAATGAGT"
Read2_revcomp = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

H2O_R1_dimer_dic = dict()
HEK_R1_dimer_dic = dict()
Twist_R1_dimer_dic = dict()
LTRI_R1_dimer_dic = dict()

Dimers = dict()

Dimers["H2O_RdRp_R1_dim"] = 0
Dimers["HEK_RdRp_R1_dim"] = 0
Dimers["Twist_RdRp_R1_dim"] = 0
Dimers["LTRI_RdRp_R1_dim"] = 0

Dimers["H2O_E_R1_dim"] = 0
Dimers["HEK_E_R1_dim"] = 0
Dimers["Twist_E_R1_dim"] = 0
Dimers["LTRI_E_R1_dim"] = 0

Dimers["H2O_N_R1_dim"] = 0
Dimers["HEK_N_R1_dim"] = 0
Dimers["Twist_N_R1_dim"] = 0
Dimers["LTRI_N_R1_dim"] = 0

Dimers["H2O_Srbd_R1_dim"] = 0
Dimers["HEK_Srbd_R1_dim"] = 0
Dimers["Twist_Srbd_R1_dim"] = 0
Dimers["LTRI_Srbd_R1_dim"] = 0

Dimers["H2O_Spoly_R1_dim"] = 0
Dimers["HEK_Spoly_R1_dim"] = 0
Dimers["Twist_Spoly_R1_dim"] = 0
Dimers["LTRI_Spoly_R1_dim"] = 0

Dimers["H2O_PPIB_R1_dim"] = 0
Dimers["HEK_PPIB_R1_dim"] = 0
Dimers["Twist_PPIB_R1_dim"] = 0
Dimers["LTRI_PPIB_R1_dim"] = 0

Dimers["H2O_unmatched_R1_dim"] = 0
Dimers["HEK_unmatched_R1_dim"] = 0
Dimers["Twist_unmatched_R1_dim"] = 0
Dimers["LTRI_unmatched_R1_dim"] = 0

pilot_H2O_raw_read_num = 0
pilot_Twist_raw_read_num = 0
pilot_HEK_raw_read_num = 0
pilot_LTRI_raw_read_num = 0
with open('filelist_R1.txt') as f:
    for line in f:
        filename = line.strip()
        it=2
        with open(filename) as f:
            for line in f:
                it = it+1
                if it%4==0 :
                    seq = line.strip()
                    if Read2_revcomp[:10] in seq:
                        if "H2O" in filename :
                            pilot_H2O_raw_read_num = pilot_H2O_raw_read_num +1
                            if seq not in H2O_R1_dimer_dic : H2O_R1_dimer_dic[seq] = 1
                            else : H2O_R1_dimer_dic[seq] = 1+H2O_R1_dimer_dic[seq]
                        elif "HEK" in filename :
                            pilot_HEK_raw_read_num = pilot_HEK_raw_read_num +1
                            if seq not in HEK_R1_dimer_dic : HEK_R1_dimer_dic[seq] = 1
                            else : HEK_R1_dimer_dic[seq] = 1+HEK_R1_dimer_dic[seq]
                        elif "Twist" in filename :
                            pilot_Twist_raw_read_num = pilot_Twist_raw_read_num +1
                            if seq not in Twist_R1_dimer_dic : Twist_R1_dimer_dic[seq] = 1
                            else : Twist_R1_dimer_dic[seq] = 1+Twist_R1_dimer_dic[seq]
                        elif "LTRI" in filename :
                            pilot_LTRI_raw_read_num = pilot_LTRI_raw_read_num +1
                            if seq not in LTRI_R1_dimer_dic : LTRI_R1_dimer_dic[seq] = 1
                            else : LTRI_R1_dimer_dic[seq] = 1+LTRI_R1_dimer_dic[seq]

for seqs in H2O_R1_dimer_dic :
    if seqs[2:12] in PPIB_F : Dimers["H2O_PPIB_R1_dim"] = Dimers["H2O_PPIB_R1_dim"]+H2O_R1_dimer_dic[seqs]
    elif seqs[2:12] in Srbd_F : Dimers["H2O_Srbd_R1_dim"] = Dimers["H2O_Srbd_R1_dim"]+H2O_R1_dimer_dic[seqs]
    elif seqs[2:12] in Spoly_F : Dimers["H2O_Spoly_R1_dim"] = Dimers["H2O_Spoly_R1_dim"]+H2O_R1_dimer_dic[seqs]
    elif seqs[2:12] in E_F : Dimers["H2O_E_R1_dim"] = Dimers["H2O_E_R1_dim"]+H2O_R1_dimer_dic[seqs]
    elif seqs[2:12] in N_F : Dimers["H2O_N_R1_dim"] = Dimers["H2O_N_R1_dim"]+H2O_R1_dimer_dic[seqs]
    elif seqs[2:12] in RdRp_F : Dimers["H2O_RdRp_R1_dim"] = Dimers["H2O_RdRp_R1_dim"]+H2O_R1_dimer_dic[seqs]
for seqs in Twist_R1_dimer_dic :
    if seqs[2:12] in PPIB_F : Dimers["Twist_PPIB_R1_dim"] = Dimers["Twist_PPIB_R1_dim"]+Twist_R1_dimer_dic[seqs]
    elif seqs[2:12] in Srbd_F : Dimers["Twist_Srbd_R1_dim"] = Dimers["Twist_Srbd_R1_dim"]+Twist_R1_dimer_dic[seqs]
    elif seqs[2:12] in Spoly_F : Dimers["Twist_Spoly_R1_dim"] = Dimers["Twist_Spoly_R1_dim"]+Twist_R1_dimer_dic[seqs]
    elif seqs[2:12] in E_F : Dimers["Twist_E_R1_dim"] = Dimers["Twist_E_R1_dim"]+Twist_R1_dimer_dic[seqs]
    elif seqs[2:12] in N_F : Dimers["Twist_N_R1_dim"] = Dimers["Twist_N_R1_dim"]+Twist_R1_dimer_dic[seqs]
    elif seqs[2:12] in RdRp_F : Dimers["Twist_RdRp_R1_dim"] = Dimers["Twist_RdRp_R1_dim"]+Twist_R1_dimer_dic[seqs]
for seqs in HEK_R1_dimer_dic :
    if seqs[2:12] in PPIB_F : Dimers["HEK_PPIB_R1_dim"] = Dimers["HEK_PPIB_R1_dim"]+HEK_R1_dimer_dic[seqs]
    elif seqs[2:12] in Srbd_F : Dimers["HEK_Srbd_R1_dim"] = Dimers["HEK_Srbd_R1_dim"]+HEK_R1_dimer_dic[seqs]
    elif seqs[2:12] in Spoly_F : Dimers["HEK_Spoly_R1_dim"] = Dimers["HEK_Spoly_R1_dim"]+HEK_R1_dimer_dic[seqs]
    elif seqs[2:12] in E_F : Dimers["HEK_E_R1_dim"] = Dimers["HEK_E_R1_dim"]+HEK_R1_dimer_dic[seqs]
    elif seqs[2:12] in N_F : Dimers["HEK_N_R1_dim"] = Dimers["HEK_N_R1_dim"]+HEK_R1_dimer_dic[seqs]
    elif seqs[2:12] in RdRp_F : Dimers["HEK_RdRp_R1_dim"] = Dimers["HEK_RdRp_R1_dim"]+HEK_R1_dimer_dic[seqs]
for seqs in LTRI_R1_dimer_dic :
    if seqs[2:12] in PPIB_F : Dimers["LTRI_PPIB_R1_dim"] = Dimers["LTRI_PPIB_R1_dim"]+LTRI_R1_dimer_dic[seqs]
    elif seqs[2:12] in Srbd_F : Dimers["LTRI_Srbd_R1_dim"] = Dimers["LTRI_Srbd_R1_dim"]+LTRI_R1_dimer_dic[seqs]
    elif seqs[2:12] in Spoly_F : Dimers["LTRI_Spoly_R1_dim"] = Dimers["LTRI_Spoly_R1_dim"]+LTRI_R1_dimer_dic[seqs]
    elif seqs[2:12] in E_F : Dimers["LTRI_E_R1_dim"] = Dimers["LTRI_E_R1_dim"]+LTRI_R1_dimer_dic[seqs]
    elif seqs[2:12] in N_F : Dimers["LTRI_N_R1_dim"] = Dimers["LTRI_N_R1_dim"]+LTRI_R1_dimer_dic[seqs]
    elif seqs[2:12] in RdRp_F : Dimers["LTRI_RdRp_R1_dim"] = Dimers["LTRI_RdRp_R1_dim"]+LTRI_R1_dimer_dic[seqs]

for items in Dimers :
    fileout_R1_dimers.write(items + "," + str(Dimers[items]) + "\n")
fileout_R1_dimers.write("pilot_H2O_raw_read_num," + str(pilot_H2O_raw_read_num) + "\n")
fileout_R1_dimers.write("pilot_HEK_raw_read_num," + str(pilot_HEK_raw_read_num) + "\n")
fileout_R1_dimers.write("pilot_LTRI_raw_read_num," + str(pilot_LTRI_raw_read_num) + "\n")
fileout_R1_dimers.write("pilot_Twist_raw_read_num," + str(pilot_Twist_raw_read_num) + "\n")

fileout_H2O = open("pilot_H2O_dimers.txt", "w")
H2O_dimer_l = list()
H2O_dimer_t = H2O_R1_dimer_dic.items()
for key,val in H2O_dimer_t :
    H2O_dimer_l.append((val,key))
H2O_dimer_l.sort(reverse=True)
for val,item in H2O_dimer_l :
    fileout_H2O.write(item + "," + str(val) + "\n")

fileout_HEK = open("pilot_HEK_dimers.txt", "w")
HEK_dimer_l = list()
HEK_dimer_t = HEK_R1_dimer_dic.items()
for key,val in HEK_dimer_t :
    HEK_dimer_l.append((val,key))
HEK_dimer_l.sort(reverse=True)
for val,item in HEK_dimer_l :
    fileout_HEK.write(item + "," + str(val) + "\n")

fileout_Twist = open("pilot_Twist_dimers.txt", "w")
Twist_dimer_l = list()
Twist_dimer_t = Twist_R1_dimer_dic.items()
for key,val in Twist_dimer_t :
    Twist_dimer_l.append((val,key))
Twist_dimer_l.sort(reverse=True)
for val,item in Twist_dimer_l :
    fileout_Twist.write(item + "," + str(val) + "\n")

fileout_LTRI = open("pilot_LTRI_dimers.txt", "w")
LTRI_dimer_l = list()
LTRI_dimer_t = LTRI_R1_dimer_dic.items()
for key,val in LTRI_dimer_t :
    LTRI_dimer_l.append((val,key))
LTRI_dimer_l.sort(reverse=True)
for val,item in LTRI_dimer_l :
    fileout_LTRI.write(item + "," + str(val) + "\n")
