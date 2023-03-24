import plink_data_analysis as pda
import os

if __name__ == "__main__":
    protein_list = ["BSA", "BQ", "LBM"]
    crosslink_report_root_list = [
        "G:/MSdata/230317BQMC/BSA/reports/crosslink/",
        "G:/MSdata/230306BQM/BQ3/reports/crosslink/",
        "G:/MSdata/230306BQM/LBM3/reports/crosslink/"
    ]
    monolink_report_root_list = [
        "G:/MSdata/230306BQM/BSAmono/reports/monolink/",
        "G:/MSdata/230306BQM/BQmono/reports/monolink/",
        "G:/MSdata/230306BQM/LBMmono/reports/monolink/"
    ]
    looplink_report_root_list = [
        "G:/MSdata/230317BQMC/BSA/reports/looplink/",
        "G:/MSdata/230306BQM/BQ3/reports/looplink/",
        "G:/MSdata/230306BQM/LBM3/reports/looplink/"
    ]
    crosslink_report_file_list = [
        "G:/MSdata/230317BQMC/BSA/reports/bsa_con_2023.03.17.filtered_cross-linked_spectra.csv",
        "G:/MSdata/230306BQM/BQ3/reports/conalbumin_con_2023.03.08.filtered_cross-linked_spectra.csv",
        "G:/MSdata/230306BQM/LBM3/reports/Lactoferrin_con_2023.03.09.filtered_cross-linked_spectra.csv"
    ]
    monolink_report_file_list = [
        "G:/MSdata/230306BQM/BSAmono/reports/bsa_con_2023.03.07.filtered_mono-linked_spectra.csv",
        "G:/MSdata/230306BQM/BQmono/reports/conalbumin_con_2023.03.07.filtered_mono-linked_spectra.csv",
        "G:/MSdata/230306BQM/LBMmono/reports/Lactoferrin_con_2023.03.07.filtered_mono-linked_spectra.csv"
    ]
    looplink_report_file_list = [
        "G:/MSdata/230317BQMC/BSA/reports/bsa_con_2023.03.17.filtered_loop-linked_spectra.csv",
        "G:/MSdata/230306BQM/BQ3/reports/conalbumin_con_2023.03.08.filtered_loop-linked_spectra.csv",
        "G:/MSdata/230306BQM/LBM3/reports/Lactoferrin_con_2023.03.09.filtered_loop-linked_spectra.csv"
    ]
    fasta_root = "G:/MSdata/fasta/"
    fasta_name_list = ["bsa", "conalbumin", "Lactoferrin"]
    pdb_root = "G:/MSdata/pdb/"
    for i in [0, 1, 2]:
        with open(fasta_root + fasta_name_list[i] + ".fasta", "r") as f:
            fasta = "".join([i.strip() for i in f][1:])

        # pda.calc_fasta(fasta, table_form=True)

        for j in [
            [crosslink_report_file_list, crosslink_report_root_list],
                # [monolink_report_file_list, monolink_report_root_list],
                # [looplink_report_file_list, looplink_report_root_list]
        ]:
            # pda.plink2normalform(j[0][i], j[1][i])
            all_filename_list = os.listdir(j[1][i])
            all_filename_list.sort(
                key=lambda x: int(x.split("_")[1]) * 100 + int(
                    x.split(protein_list[i])[-1].strip("_").strip(".csv")))

            output_list = [all_filename_list, [], []]
            for k in all_filename_list:
                output_list[1].append(
                    pda.plink_report_valid_link_sfDR(
                        j[1][i] + k,
                        fasta=fasta,
                        least_ssn=1,
                        type="crosslink",
                        score_count_mode=False)[0])
            # routine = int(len(output_list[0]) / 2)
            # for k in range(routine * 2):
            #     output_list[2].append(
            #         pda.cal_multi_list([
            #             output_list[1][k],
            #             output_list[1][k],
            #         ]))

            print(all_filename_list)
            print(
                '|-A-|-C-|-D-|-E-|-F-|-G-|-H-|-I-|-K-|-L-|-M-|-N-|-P-|-Q-|-R-|-S-|-T-|-V-|-W-|-Y-|ALL|'
            )
            print('|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|---|')
            for k in output_list[1]:
                # pda.report_animo_ratio(k,
                #                        pos_count_form=False,
                #                        pos_form=True,
                #                        fasta_file=fasta,
                #                        tableform=True,
                #                        ratio_output=False,
                #                        threshold=10000)

                link_list = []
                report_list = []
                dis_list = []
                for p in k:
                    link_list.append(
                        [fasta[int(p[0]) - 1], fasta[int(p[1]) - 1]])
                for amino in ["A", "L", "V"]:
                    for cross_dimer in range(len(link_list)):
                        if link_list[cross_dimer][0] == "K" or link_list[
                                cross_dimer][1] == "K":
                            if link_list[cross_dimer][0] == "K" and link_list[
                                    cross_dimer][1] == amino:
                                report_list.append(k[cross_dimer])
                            elif link_list[cross_dimer][1] == "K" and link_list[
                                    cross_dimer][0] == amino:
                                report_list.append(k[cross_dimer][::-1])
                print(
                    pda.cal_distance_pos_list(
                        report_list,
                        fasta,
                        pdb_file=pdb_root + fasta_name_list[i] + ".pdb",
                        animo_core1='CA',
                        animo_core2='CA',
                        print_all=False,
                    )[2:])
                # print(dis_list)
