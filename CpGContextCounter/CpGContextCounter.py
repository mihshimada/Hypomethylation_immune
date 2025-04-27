f_1 = open("EPIC_annotation.txt", "r", encoding = "utf-8")
w_1 = open("EPIC_context.txt", "w", encoding = "utf-8")


filename = []

tmp = ["ProbeID", "CpG_count", "CpG_category", "seq4", "seq4_category"]
list_w = '\t'.join(tmp)
list_w2 = list_w + "\n"
w_1.writelines(list_w2)

i = 1
for file_1 in f_1:
    line_1 = file_1.strip().split('\t')

    if "cg" in line_1[0] and len(line_1) >=2:
        seq_F = line_1[1].split("[CG]")[0][-35:]
        seq_R = line_1[1].split("[CG]")[1][:35]

        CpG_count = seq_F.count('CG') + seq_R.count('CG')
        if CpG_count == 0:
            CpG_category = 0
        elif CpG_count == 1:
            CpG_category = 1
        elif CpG_count == 2:
            CpG_category = 2
        elif CpG_count >= 3:
            CpG_category = 3

        Pre_CG = line_1[1].split("[CG]")[0][-1]
        Post_CG = line_1[1].split("[CG]")[1][0]

        seq4 = Pre_CG + "CG" + Post_CG

        if Pre_CG in ["A","T"] and Post_CG in ["A","T"]:
            seq4_category = "WCGW"
        elif Pre_CG in ["C","G"] and Post_CG in ["C","G"]:
            seq4_category = "SCGS"
        elif Pre_CG in ["A","T"] and Post_CG in ["C","G"]:
            seq4_category = "SCGW"
        elif Pre_CG in ["C","G"] and Post_CG in ["A","T"]:
            seq4_category = "SCGW"

        tmp = [line_1[0], str(CpG_count), str(CpG_category), seq4, seq4_category]
        list_w = '\t'.join(tmp)
        list_w2 = list_w + "\n"
        w_1.writelines(list_w2)


    print(i)
    i += 1
