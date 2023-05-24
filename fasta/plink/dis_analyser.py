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

    for pdb_name in ["LF-AF-Q6LBN7-F1-model_v4"]:
        for i in [2]:
            with open(fasta_root + fasta_name_list[i] + ".fasta", "r") as f:
                fasta = "".join([i.strip() for i in f][1:])
            for j in [[crosslink_report_root_list],
                      [looplink_report_root_list]]:
                all_filename_list = os.listdir(j[0][i])
                all_filename_list.sort(
                    key=lambda x: int(x.split("_")[1]) * 100 + int(
                        x.split(protein_list[i])[-1].strip("_").strip(".csv")))

                for k in all_filename_list:
                    data = []
                    with open(j[0][i] + k, "r") as file:
                        csv_reader = csv.reader(file)
                        data = [row for row in csv_reader]
                    with open(j[0][i] + k, "w", newline='') as file:
                        writer = csv.writer(file)
                        for l in range(len(data)):
                            if l == 0:
                                data[l].append(pdb_name)
                                writer.writerow(data[l])
                            else:
                                try:
                                    data[l].append(
                                        pda.cal_distance_pos(int(data[l][26]),
                                                             int(data[l][32]),
                                                             pdb_root +
                                                             pdb_name + ".pdb",
                                                             autoprint=True))
                                except:
                                    data[l].append("")
                                writer.writerow(data[l])
