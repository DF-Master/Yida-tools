bsa_fasta = "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLTSSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"


def analyze_neighborhood(fasta, target='K', distance=7, percent=False):
    animo_list = list("ACDEFGHIKLMNPQRSTVWY")
    cross_dic = {}
    for animo_tag in animo_list:
        cross_dic[animo_tag] = 0
    cross_dic["all"] = 0
    cross_dic["nofind"] = 0

    fasta_list = list(fasta)
    for i in range(len(fasta)):
        if fasta[i] == target:
            for j in range(i - distance, i + distance + 1):
                try:
                    if fasta_list[j] in animo_list:
                        cross_dic[fasta_list[j]] += 1
                        cross_dic["all"] += 1
                    else:
                        cross_dic['nofind'] += 1
                        cross_dic["all"] += 1
                except:
                    continue
    if percent == False:
        return (cross_dic)
    else:
        for animo_tag in animo_list:
            cross_dic[animo_tag] = cross_dic[animo_tag] / cross_dic['all']
        return (cross_dic)


if __name__ == '__main__':
    print(analyze_neighborhood(bsa_fasta, percent=True))
