import plink_data_analysis as pda
import settings as set
import os, csv

if __name__ == "__main__":
    protein_list = ["BSA", "BQ", "LBM"]
    crosslink_report_root_list = set.crosslink_report_withdis_root_list
    looplink_report_root_list = set.looplink_report_withdis_root_list

    fasta_root = set.fasta_root
    fasta_name_list = ["bsa", "conalbumin", "Lactoferrin"]
    pdb_root = set.pdb_root

    for i in [0]:
        with open(fasta_root + fasta_name_list[i] + ".fasta", "r") as f:
            fasta = "".join([i.strip() for i in f][1:])
        for j in [[looplink_report_root_list]]:
            all_filename_list = os.listdir(j[0][i])
            all_filename_list.sort(
                key=lambda x: int(x.split("_")[1]) * 100 + int(
                    x.split(protein_list[i])[-1].strip("_").strip(".csv")))

            # output_list = [all_filename_list, [], []]
            for k in all_filename_list:
                data = []
                with open(j[0][i] + k, "r") as file:
                    csv_reader = csv.reader(file)
                    data = [row for row in csv_reader]
                    # crosslink_pos_values = [
                    #     pda.find_link_pos(row[4], fasta) for row in data
                    #     if "-" in row[4]  #CrossLink
                    # ]
                    crosslink_pos_values = [
                        pda.find_link_pos(
                            str(")-" + row[4].split("(")[0] + "(").join(
                                row[4].split(")(")), fasta) for row in data
                        if ")(" in row[4]  #LoopLink
                    ]
                    for l in range(len(data)):
                        if l == 0:
                            for s in [
                                    "Pep1_fasta", "Pep1_num", "Cross_pos_1",
                                    "Pep1_start", "Pep1_end", "CL1_fasta_num",
                                    "Pep2_fasta", "Pep2_num", "Cross_pos_2",
                                    "Pep2_start", "Pep2_end", "CL2_fasta_num"
                            ]:
                                data[l].append(s)
                        else:
                            for s in crosslink_pos_values[
                                    l - 1][0] + crosslink_pos_values[l - 1][1]:
                                data[l].append(s)
                with open(j[0][i] + k, "w", newline='') as file:
                    writer = csv.writer(file)
                    for l in range(len(data)):
                        if l == 0:
                            writer.writerow(data[l])
                        elif int(data[l][-2]) != 1 and int(data[l][-8]) != 1:
                            # if str(data[l][-4]) == "K" or str(
                            #         data[l][-10]) == "K":
                            writer.writerow(data[l])
                        else:
                            print("Error Pep")
