# Change the working directory
import os
import pandas as pd
import re

#feature list file
feature_file = "feature_file.txt" # <-- Put your own file name here

# Open input and output files
f_0 = open("files.txt", "r", encoding="utf-8")
o_0 = open("count_result.txt", "w", encoding="utf-8")
list_w = ["feature", "inside", "outside", "\n"]
list_w_2 = '\t'.join(list_w)
o_0.writelines(list_w_2)

for file_1 in f_0:
    file_2 = file_1.strip().split('\t')
    file_name = file_2[0]
    item_name = file_1.strip().split('_')[0]
    out_name = item_name + "_inside_features.txt"
    print(item_name)

    #--------------------------------------------------------- Load ENCODE annotation file
    f_1 = pd.read_table(file_name, header=None, sep=r'\s+')
    f_1.columns = ["Chromosome_int", "start", "end"]

    # Sort the ENCODE file
    f_1_s = f_1.sort_values(['Chromosome_int', 'start'])
    chr_list_ENCODE = set(f_1_s['Chromosome_int'].values.tolist())

    #--------------------------------------------------------- Load and process feature data
    f_2 = pd.read_table(feature_file, header=None)
    f_2.columns = ["featureID"]

    f_2_2 = pd.concat([f_2, f_2['featureID'].str.split(':', expand=True)], axis=1)
    f_2_3 = pd.concat([f_2_2, f_2_2[1].str.split('_', expand=True)], axis=1).drop(1, axis=1)
    f_2_3.columns = ['featureID', 'chr', 'position']
    f_2_3['chr'] = f_2_3['chr'].astype(int)
    f_2_3['position'] = f_2_3['position'].astype(int)

    f_2_s = f_2_3.sort_values(['chr', 'position'])

    # Store the features into a multi-dimensional list by chromosome
    feature_list = []
    chr_list_feature = set(f_2_s['chr'].values.tolist())
    for i in range(1, 23):
        if i in chr_list_feature:
            df_tmp = f_2_s.groupby("chr").get_group(i)
            list_tmp = df_tmp.values.tolist()
            feature_list.append(list_tmp)
        else:
            feature_list.append([[0, 0, 0]])

    #--------------------------------------------------------- Compare features with ENCODE regions
    o_1 = open(out_name, 'w')

    feature_list_no = 1
    feature_outside = 0
    feature_inside = 0

    for i in range(22):
        if i + 1 in chr_list_ENCODE:
            f_1_s2 = f_1_s.groupby("Chromosome_int").get_group(i + 1)
        else:
            f_1_s2 = pd.DataFrame({'Chromosome': [0], 'start': [0], 'end': [0]})

        num = 0
        encode_line = f_1_s2[num:num + 1]
        pos_1 = int(encode_line['start'].iloc[0])
        pos_2 = int(encode_line['end'].iloc[0])

        for item_feature in feature_list[i]:
            if item_feature[0] != 0:
                pos_feature = int(item_feature[2])

                if pos_feature < pos_1:
                    feature_outside += 1
                elif pos_1 <= pos_feature <= pos_2:
                    feature_inside += 1
                    print(item_feature[0], file=o_1)
                elif pos_feature > pos_2:
                    if num < len(f_1_s2) - 1:
                        while pos_feature > pos_2 and num < len(f_1_s2) - 1:
                            num += 1
                            encode_line = f_1_s2[num:num + 1]
                            pos_1 = int(encode_line['start'].iloc[0])
                            pos_2 = int(encode_line['end'].iloc[0])
                        if pos_feature < pos_1:
                            feature_outside += 1
                        elif pos_1 <= pos_feature <= pos_2:
                            feature_inside += 1
                            print(item_feature[0], file=o_1)
                        else:
                            feature_outside += 1
                    else:
                        if pos_feature < pos_1:
                            feature_outside += 1
                        elif pos_1 <= pos_feature <= pos_2:
                            feature_inside += 1
                            print(item_feature[0], file=o_1)
                        else:
                            feature_outside += 1

                feature_list_no += 1

    print("feature_outside:", feature_outside)
    print("feature_inside:", feature_inside)

    list_w = [item_name, str(feature_inside), str(feature_outside), "\n"]
    list_w_2 = '\t'.join(list_w)
    o_0.writelines(list_w_2)

# Close all files
o_0.close()
f_0.close()
o_1.close()
