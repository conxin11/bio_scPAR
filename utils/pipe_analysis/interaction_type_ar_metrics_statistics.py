import numpy as np
import pandas as pd
from scipy import stats

# read network data
from utils.data_process.read_regnetworks import RegulationNetwork
from utils.data_process.read_ppi import PPI
from utils.data_process.read_reactome import ReactomeNetwork

from config import Config

# read ar result
from utils.data_process.ar_metrics_process import LoadArMetrics, FilterArMetrics, SaveArMetrics

class Statistics:
    '''
    data:{"group1":[],"group2":[]}
    '''
    def __init__(self, data):
        self.data = data
        self.result = None

    def t_test(self, group1, group2):
        tstat, pval = stats.ttest_ind(group1, group2, alternative="two-sided")
        return pval

    def ranksum_test(self, group1, group2):
        stat, pval = stats.ranksums(group1, group2, alternative="two-sided")
        return pval

    def calculate_statistics(self, method="t_test"):
        group1 = self.data["group1"]
        group2 = self.data["group2"]
        if method == "t_test":
            pval = self.t_test(group1, group2)
        elif method == "ranksum_test":
            pval = self.ranksum_test(group1, group2)
        else:
            raise ValueError("method should be t_test or ranksum_test")
        return pval


class NetworkIntegration:
    def __init__(self):
        self.networks = {}

    def read_regulation_network(self, tf_fp, regnetwork_fp):
        regnetworks = RegulationNetwork(tf_fp=tf_fp, interaction_fp=regnetwork_fp)
        regnetworks.create_network()
        return regnetworks.data_full_columns

    def read_reactome_network(self, gmt_fp, reactome_fp):
        reactome = ReactomeNetwork(gmt_fp=gmt_fp, interaction_fp=reactome_fp)
        reactome.prepare()
        data = reactome.get_data()
        return data

    def read_ppi_network(self, ppi_fp, ppi_mapping_fp):
        ppi = PPI(mapping_fp=ppi_mapping_fp, ppi_fp=ppi_fp)
        ppi.read_mapper()
        ppi.read_ppi()
        ppi.map_id_to_gene()
        return ppi.interactions

    def read_all_networks(self):
        # read tf network
        tf_fp = Config.tf_fp
        regnetwork_fp = Config.regnetwork_fp
        tf_network = self.read_regulation_network(tf_fp, regnetwork_fp)

        # read reactome network
        gmt_fp = Config.gmt_fp
        reactome_fp = Config.reactome_fp
        reactome_network = self.read_reactome_network(gmt_fp, reactome_fp)

        # read ppi network
        # ppi_fp = Config.ppi_fp
        # ppi_mapping_fp = Config.ppi_mapping_fp
        # ppi_network = self.read_ppi_network(ppi_fp, ppi_mapping_fp)

        self.networks["tf"] = tf_network
        self.networks["reactome"] = reactome_network
        # self.networks["ppi"] = ppi_network
        return self.networks


class CalculateInteractionTypeARMetrics:
    def __init__(self, threshold:dict, ar_metrics_fp:str):
        self.network_integration = NetworkIntegration()
        self.all_networks = self.network_integration.read_all_networks()
        self.threshold = threshold
        self.ar_metrics_fp = ar_metrics_fp
        self.all_metrics = [k for k in Config.threshold_dict.keys()]
        self.result = {}
        self.result_df = pd.DataFrame(columns=["cell_type", "interaction_database", "metrics","interaction_type", "pval"])

    def filter_ar_metrics(self):
        ar_metrics = LoadArMetrics(self.ar_metrics_fp).load_result()
        filter_ar_metrics = FilterArMetrics(ar_metrics, threshold_dict=self.threshold)
        result = filter_ar_metrics.filter_pairs(results=filter_ar_metrics.results,
                                                 threshold_dict_tem=self.threshold)
        return result

    # def ar_result_set_index(self, ar_metrics):
    #     ar_metrics = ar_metrics.set_index(["antecedent", "consequent"],drop=False)
    #     return ar_metrics

    def calculate_ar_metrics(self, cell_type, ar_metrics):

        ar_metrics = ar_metrics.set_index(["antecedent", "consequent"],drop=False)

        # reactome network
        reactome_network = self.all_networks["reactome"]
        reactome_network = reactome_network.set_index(["Gene1", "Gene2"], drop=False)
        index_intersection = reactome_network.index.intersection(ar_metrics.index)
        ar_metrics_reactome = ar_metrics.loc[index_intersection]
        # add direction column
        Direction = reactome_network.loc[index_intersection,["Direction"]]
        Direction = Direction.reset_index(drop=False)
        #remove duplicate rows
        Direction = Direction.drop_duplicates(subset=["Gene1", "Gene2"])
        Direction = Direction.set_index(["Gene1", "Gene2"])

        try:
            ar_metrics_reactome["Direction"] = Direction.loc[ar_metrics_reactome.index, "Direction"]
        except:
            print("ar_metrics_reactome",len(ar_metrics_reactome), "reactome", len(Direction))
        else:
            interaction_types = set( ar_metrics_reactome["Direction"])
            self.result[cell_type] = {}
            self.result[cell_type]["reactome"] = {}
            for metrics in self.all_metrics:
                for interaction_type in interaction_types:

                    _tem1 = ar_metrics_reactome[ar_metrics_reactome["Direction"] == interaction_type]
                    ar_metrics_reactome_metrics1 = _tem1[metrics]
                    group1 = ar_metrics_reactome_metrics1.values.tolist()

                    _tem2 = ar_metrics_reactome[ar_metrics_reactome["Direction"] != interaction_type]
                    ar_metrics_reactome_metrics2 = _tem2[metrics]
                    group2 = ar_metrics_reactome_metrics2.values.tolist()

                    if group1 and group2:
                        statistics = Statistics(data={"group1":group1, "group2":group2})
                        pval = statistics.calculate_statistics(method="ranksum_test")
                        self.result[cell_type]["reactome"][metrics + "_" + interaction_type] = pval

                        #add a row to result_df
                        # self.result_df = self.result_df.append({"cell_type":cell_type, "interaction_database":"reactome", "metrics":metrics , "interaction_type":interaction_type, "pval":pval}, ignore_index=True)
                        self.result_df = pd.concat([self.result_df, pd.DataFrame([{"cell_type": cell_type,
                                                                                   "interaction_database": "reactome",
                                                                                   "metrics": metrics,
                                                                                   "interaction_type": interaction_type,
                                                                                   "pval": pval}])], ignore_index=True)


        # tf network
        # TF	Target	Type Up_or_Down_or_Unknown
        tf_network = self.all_networks["tf"]
        tf_network = tf_network.set_index(["TF", "Target"], drop=False)
        index_intersection = tf_network.index.intersection(ar_metrics.index)
        ar_metrics_tf = ar_metrics.loc[index_intersection]
        # add direction column
        UpDown_tf = tf_network.loc[index_intersection,["Up_or_Down_or_Unknown"]]
        UpDown_tf = UpDown_tf.reset_index(drop=False)
        #remove duplicate rows
        UpDown_tf = UpDown_tf.drop_duplicates(subset=["TF", "Target"])
        UpDown_tf = UpDown_tf.set_index(["TF", "Target"])

        try:
            ar_metrics_tf["Up_or_Down_or_Unknown"] = UpDown_tf.loc[ar_metrics_tf.index, "Up_or_Down_or_Unknown"]
        except:
            print("ar_metrics_tf",len(ar_metrics_tf), "UpDown_tf", len(UpDown_tf))
        else:
            interaction_types = set(ar_metrics_tf["Up_or_Down_or_Unknown"])
            self.result[cell_type]["tf_up_down"] = {}
            for metrics in self.all_metrics:
                for interaction_type in interaction_types:
                    _tem1 = ar_metrics_tf[ar_metrics_tf["Up_or_Down_or_Unknown"] == interaction_type]
                    ar_metrics_tf_metrics1 = _tem1[metrics]
                    group1 = ar_metrics_tf_metrics1.values.tolist()

                    _tem2 = ar_metrics_tf[ar_metrics_tf["Up_or_Down_or_Unknown"] != interaction_type]
                    ar_metrics_tf_metrics2 = _tem2[metrics]
                    group2 = ar_metrics_tf_metrics2.values.tolist()
                    if group1 and group2:
                        statistics = Statistics(data={"group1":group1, "group2":group2})
                        pval = statistics.calculate_statistics(method="ranksum_test")
                        self.result[cell_type]["tf_up_down"][metrics + "_" + interaction_type] = pval

                        #add a row to result_df
                        # self.result_df = self.result_df.append({"cell_type":cell_type, "interaction_database":"tf_up_down", "metrics":metrics , "interaction_type":interaction_type, "pval":pval}, ignore_index=True)

                        self.result_df = pd.concat([self.result_df, pd.DataFrame([{"cell_type": cell_type,
                                                                                   "interaction_database": "tf_up_down",
                                                                                   "metrics": metrics,
                                                                                   "interaction_type": interaction_type,
                                                                                   "pval": pval}])], ignore_index=True)

        #add type column
        Type_tf = tf_network.loc[index_intersection,["Type"]]
        Type_tf = Type_tf.reset_index(drop=False)
        #remove duplicate rows
        Type_tf = Type_tf.drop_duplicates(subset=["TF", "Target"])
        Type_tf = Type_tf.set_index(["TF", "Target"])

        try:
            ar_metrics_tf["Type"] = Type_tf.loc[ar_metrics_tf.index, "Type"]
        except:
            print("ar_metrics_tf",len(ar_metrics_tf), "Type_tf", len(Type_tf))
        else:
            interaction_types = set(ar_metrics_tf["Type"])
            self.result[cell_type]["tf_type"] = {}
            for metrics in self.all_metrics:
                for interaction_type in interaction_types:
                    _tem1 = ar_metrics_tf[ar_metrics_tf["Type"] == interaction_type]
                    ar_metrics_tf_metrics1 = _tem1[metrics]
                    group1 = ar_metrics_tf_metrics1.values.tolist()

                    _tem2 = ar_metrics_tf[ar_metrics_tf["Type"] != interaction_type]
                    ar_metrics_tf_metrics2 = _tem2[metrics]
                    group2 = ar_metrics_tf_metrics2.values.tolist()
                    if group1 and group2:
                        statistics = Statistics(data={"group1":group1, "group2":group2})
                        pval = statistics.calculate_statistics(method="ranksum_test")
                        self.result[cell_type]["tf_type"][metrics + "_" + interaction_type] = pval

                        #add a row to result_df
                        # self.result_df = self.result_df.append(
                        #     {"cell_type": cell_type, "interaction_database": "tf_type", "metrics": metrics,
                        #      "interaction_type": interaction_type, "pval": pval}, ignore_index=True)
                        self.result_df = pd.concat([self.result_df, pd.DataFrame([{"cell_type": cell_type,
                                                                                   "interaction_database": "tf_type",
                                                                                   "metrics": metrics,
                                                                                   "interaction_type": interaction_type,
                                                                                   "pval": pval}])], ignore_index=True)

    # save result
    def save_result(self, fp="result.txt"):
        # write self.result to text file
        with open(fp, "w") as f:
            for cell_type, cell_type_result in self.result.items():
                f.write(cell_type)
                f.write("\n")
                f.write("%%%"*30)
                f.write("\n")
                f.write("%%%" * 30)
                f.write("\n")
                for network, network_result in cell_type_result.items():
                    f.write(network)
                    f.write("\n")
                    f.write("==="*30)
                    f.write("\n")
                    for metrics, pval in network_result.items():
                        f.write(metrics)
                        f.write("\n")
                        f.write(str(pval))
                        f.write("\n")
                    f.write("\n")
                f.write("\n")

    # save self.result_df to tsv file
    def save_result_df(self, fp="result_df.tsv"):
        self.result_df.to_csv(fp, sep="\t", index=False)




    def run_cell_types_ar_metrics(self, results):
        # cell types
        for ct, ct_metrics in results.items():
            print(ct)
            print("==="*30)
            self.calculate_ar_metrics(cell_type=ct, ar_metrics=ct_metrics)



if __name__ == '__main__':
    # calculate interaction type ar metrics
    threshold_dict = {"support":[0.05, None], "confidence":[0.5, None], "lift":[1, None], "leverage":[None, None], "conviction":[1, None]}
    ar_metrics_fp = "/home/liuyq/data/ar_result/ar/result_edge_style_filted"
    calcu_interaction_type_ar_metrics = CalculateInteractionTypeARMetrics(threshold=threshold_dict, ar_metrics_fp=ar_metrics_fp)
    results = calcu_interaction_type_ar_metrics.filter_ar_metrics()
    calcu_interaction_type_ar_metrics.run_cell_types_ar_metrics(results)
    print(calcu_interaction_type_ar_metrics.result)
    output_fp = "/home/liuyq/data/ar_result/ar/test_inter_type_statistics.txt"
    calcu_interaction_type_ar_metrics.save_result(output_fp)
    calcu_interaction_type_ar_metrics.save_result_df(output_fp.replace(".txt", ".tsv"))















