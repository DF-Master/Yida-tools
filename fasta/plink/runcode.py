import plink_data_analysis as pda
import settings as set
import os

if __name__ == "__main__":
    protein_list = ["BSA", "BQ", "LBM"]
    crosslink_report_root_list = set.crosslink_report_withdis_root_list
    monolink_report_root_list = set.monolink_report_root_list
    looplink_report_root_list = set.looplink_report_withdis_root_list
    fasta_root = set.fasta_root
    fasta_name_list = ["bsa", "conalbumin", "Lactoferrin"]
    pdb_root = set.pdb_root
    for i in [0]:
        with open(fasta_root + fasta_name_list[i] + ".fasta", "r") as f:
            fasta = "".join([i.strip() for i in f][1:])
        for j in [[crosslink_report_root_list], [looplink_report_root_list]]:
            all_filename_list = os.listdir(j[0][i])
            all_filename_list.sort()
            # key=lambda x: int(x.split("_")[1]) * 100 + int(
            #     x.split(protein_list[i])[-1].strip("_").strip(".csv")))

            output_list = [all_filename_list, [], []]
            for k in all_filename_list:
                output_list[1].append(
                    pda.plink_report_valid_link_sfDR(j[0][i] + k,
                                                     fasta=fasta,
                                                     least_ssn=1,
                                                     ms_peptide_list=4,
                                                     type="looplink",
                                                     score_count_mode=False)[1]
                )  # Same for CL and LL,[1] stands for report_list, and [0] stands for hard_loc_list

            print(all_filename_list)
            print(
                '|-A-|-C-|-D-|-E-|-F-|-G-|-H-|-I-|-K-|-L-|-M-|-N-|-P-|-Q-|-R-|-S-|-T-|-V-|-W-|-Y-|ALL|'
            )
            print('|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|---|')
            for k in output_list[1]:
                pda.report_animo_ratio(k,
                                       pos_count_form=False,
                                       pos_form=True,
                                       fasta_file=fasta,
                                       tableform=True,
                                       ratio_output=True,
                                       threshold=10000)

                # link_list = []
                # report_list = []
                # dis_list = []
                # for p in k:
                #     link_list.append(
                #         [fasta[int(p[0]) - 1], fasta[int(p[1]) - 1]])
                # for amino in ["A", "L", "V"]:
                #     for cross_dimer in range(len(link_list)):
                #         if link_list[cross_dimer][0] == "K" or link_list[
                #                 cross_dimer][1] == "K":
                #             if link_list[cross_dimer][0] == "K" and link_list[
                #                     cross_dimer][1] == amino:
                #                 report_list.append(k[cross_dimer])
                #             elif link_list[cross_dimer][1] == "K" and link_list[
                #                     cross_dimer][0] == amino:
                #                 report_list.append(k[cross_dimer][::-1])
                # print(
                #     pda.cal_distance_pos_list(
                #         report_list,
                #         fasta,
                #         pdb_file=pdb_root + fasta_name_list[i] + ".pdb",
                #         animo_core1='CA',
                #         animo_core2='CA',
                #         print_all=False,
                #     )[2:])
