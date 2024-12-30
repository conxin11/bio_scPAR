


import os
import pandas as pd
import re
'''
go富集的结果包括：本文方法，celltype_go_enrich.tsv; deg方法，celltype_gsea_go_deg.tsv
reactome富集的结果包括：本文方法，celltype_reactome_enrich.tsv; deg方法，celltype_gsea_reactome_deg.tsv
分别统计两种富集方法中，Term列涉及到的条目中，各自的数量，交集数量，交集占各自的比例，并集的数量，各自特异性的数量，以及总的数量，输出的结果以pandas的数据框的形式保存。
'''

class GSEASummary:
    def __init__(self, ar_fp, deg_fp):
        self.ar_fp = ar_fp
        self.deg_fp = deg_fp
        self.ar = None
        self.deg = None
        self.terms_ar = None
        self.terms_deg = None
        self.result = None

    @staticmethod
    def read_data(fp):
        data = pd.read_csv(fp, header=0, index_col=False, sep='\t')
        term_set = data.loc[:, 'Term']
        term_set = set(term_set)
        return term_set

    def statistics(self,prefix=None):
        # ar = pd.read_csv(self.ar_fp, sep='\t')
        # deg = pd.read_csv(self.deg_fp, sep='\t')
        terms_ar = self.read_data(self.ar_fp)
        terms_deg = self.read_data(self.deg_fp)
        terms_ar_num = len(terms_ar)
        terms_deg_num = len(terms_deg)
        terms_ar_deg = terms_ar & terms_deg
        terms_ar_deg_num = len(terms_ar_deg)
        terms_ar_specific = terms_ar - terms_ar_deg
        terms_deg_specific = terms_deg - terms_ar_deg
        terms_ar_specific_num = len(terms_ar_specific)
        terms_deg_specific_num = len(terms_deg_specific)
        terms_union = terms_ar | terms_deg
        terms_union_num = len(terms_union)
        if terms_ar_num > 0:
            terms_rate_intersecton_in_ar = terms_ar_deg_num / terms_ar_num
        else:
            terms_rate_intersecton_in_ar = 0
        if terms_deg_num > 0:
            terms_rate_intersecton_in_deg = terms_ar_deg_num / terms_deg_num
        else:
            terms_rate_intersecton_in_deg = 0
        result = pd.DataFrame(
            columns=['AR', 'DEG', 'Intersection', 'Union', 'AR_Specific', 'DEG_Specific', 'AR_in_DEG_Rate', 'DEG_in_AR_Rate']
        )
        result.loc[0, 'AR'] = terms_ar_num
        result.loc[0, 'DEG'] = terms_deg_num
        result.loc[0, 'Intersection'] = terms_ar_deg_num
        result.loc[0, 'Union'] = terms_union_num
        result.loc[0, 'AR_Specific'] = terms_ar_specific_num
        result.loc[0, 'DEG_Specific'] = terms_deg_specific_num
        result.loc[0, 'Intersecton_in_ar_Rate'] = terms_rate_intersecton_in_ar
        result.loc[0, 'Intersecton_in_deg_Rate'] = terms_rate_intersecton_in_deg
        if prefix is not None:
            result.columns = [prefix + '_' + col for col in result.columns]
        self.result = result

'''
首先从结果文件中，分离出以tsv为后缀的文件，
这些文件的文件名都是以"_"作为分割的，
根据分割后的list中第一个元素的值，作为细胞类型，保存在self.result中，以细胞类型为键，文件名为值
根据分割后的list中是否含有“reactome”或“go”作为富集类型 enrich_type，将文件分为两类，分别存在self.result[cell_type]的“go”和“reactome”
然后分别以“ar”和“deg”作为键，将文件分为两类，分别存在self.result[cell_type][enrich_type]的“ar”和“deg”
'''
class ExtractFile:
    def __init__(self, result_dir, statistics_fp="statistics.tsv"):
        self.result_dir = result_dir
        self.result = {}
        self.resutl_all_cell_pd = {}
        self.statistics_fp = os.path.join(result_dir, statistics_fp)
    def extract(self):
        files = os.listdir(self.result_dir)
        for file in files:
            if file.endswith('.tsv'):
                fp = os.path.join(self.result_dir, file)
                file = file.strip('.tsv')
                strs = file.split('_')
                cell_type = strs[0]
                enrich_type = 'go' if 'go' in strs else 'reactome'
                method = 'deg' if 'deg' in strs else 'ar'
                if cell_type not in self.result:
                    self.result[cell_type] = {}
                if enrich_type not in self.result[cell_type]:
                    self.result[cell_type][enrich_type] = {}
                if method not in self.result[cell_type][enrich_type]:
                    self.result[cell_type][enrich_type][method] = [fp]
                else:
                    self.result[cell_type][enrich_type][method].append(fp)
    def pair_files_and_run_statistics(self):
        for cell_type, value in self.result.items():
            cell_pd = None
            for enrich_type, value1 in value.items():
                ar_fps = value1['ar']
                deg_fps = value1['deg']
                assert isinstance(deg_fps, list) and len(deg_fps) == 1
                for ar_fp in ar_fps:
                    for deg_fp in deg_fps:
                        gseasummary = GSEASummary(ar_fp, deg_fp)
                        prefix = self.extract_prefix(enrich_type, ar_fp)
                        gseasummary.statistics(prefix=prefix)
                        gseasummary.result.loc[:, 'cell_type'] = cell_type
                        if cell_pd is None:
                            cell_pd = gseasummary.result
                        else:
                            # concat dataframe by column
                            cell_pd = pd.concat([cell_pd, gseasummary.result], axis=1)
            cell_pd['cell_type'] = cell_type
            self.resutl_all_cell_pd[cell_type] = cell_pd
            # if self.resutl_all_cell_pd is None:
            #     self.resutl_all_cell_pd = cell_pd
            # else:
            #     self.resutl_all_cell_pd = pd.concat([self.resutl_all_cell_pd, cell_pd])

    def extract_prefix(self, enrich_type, file_path):
        fn = os.path.basename(file_path)
        if "community" in fn:
            match = re.search(r"(community_id_\d+)\.tsv", fn)
            prefix = match.group(1) if match else None
        else:
            prefix = "all"
        assert prefix is not None
        prefix = enrich_type + '_' + prefix
        return prefix

    # pandas dataframe to csv string
    def to_csv_string(self, df):
        output = df.to_csv(sep='\t', index=False)
        return output

    def save(self):
        all_cell_string = ""
        for cell_type, pd in self.resutl_all_cell_pd.items():
            df_str  = self.to_csv_string(pd)
            all_cell_string = all_cell_string + '\n' + df_str
        with open(self.statistics_fp, 'w') as f:
            f.write(all_cell_string)


if __name__ == '__main__':
    result_dir = r"/mnt/sda/liuyq/tem/gsea"
    extractfile = ExtractFile(result_dir)
    extractfile.extract()
    extractfile.pair_files_and_run_statistics()
    extractfile.save()















