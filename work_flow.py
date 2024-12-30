import config
import os
from config import Config
from utils.data_process.ar_metrics_process import SaveArMetrics, LoadArMetrics
from utils.pipe_analysis.cell_type_pipe import CellTypePipe
from utils.pipe_analysis.interaction_type_ar_metrics_statistics import CalculateInteractionTypeARMetrics, Statistics, NetworkIntegration

thread_dict = Config.threshold_dict
cell_type_pipe = CellTypePipe(data_path=Config.data_fp,
                              ar_result_path=Config.ar_fp,
                              cell_type_fp=Config.cell_type_fp,
                              threshold_dict=Config.threshold_dict)

if os.path.exists(Config.cahe_result_edge_style_fp) and os.path.exists(Config.cahe_result_edge_style_filted_fp):
    print("load from cache...")
    result_edge_style = LoadArMetrics(Config.cahe_result_edge_style_fp).load_result()
    result_edge_style_filted = LoadArMetrics(Config.cahe_result_edge_style_filted_fp).load_result()
else:

    # import data
    adata_subs, data_obj = cell_type_pipe.read_data(min_support=Config.min_support)

    # select genes
    bio_net_fp_dict = {"ppi_fp": Config.ppi_fp,
                       "ppi_mapping_fp": Config.ppi_mapping_fp,
                       "gmt_fp": Config.gmt_fp,
                       "reactome_fp": Config.reactome_fp,
                       "tf_fp": Config.tf_fp,
                       "regnetwork_fp": Config.regnetwork_fp}
    adata_subs_gene_selected = cell_type_pipe.subset_by_genes(adata_subs, bio_net_fp_dict, deg_sel=True,
                                                              network_sel=True)

    # run ar
    results = cell_type_pipe.run_ar(adata_subs=adata_subs_gene_selected)

    # filter results
    # results_filted = cell_type_pipe.filter_result(results)

    result_edge_style = cell_type_pipe.transform_to_edge_style(results)
    res_fp = os.path.join(config.Config.ar_fp, "result_edge_style")
    save_obj = SaveArMetrics(results=result_edge_style, out_path=res_fp)
    save_obj.write_result()
    print("finish writing result_edge_style")

    result_edge_style_filted = cell_type_pipe.filter_transform_to_edge_style(result_edge_style)
    res_filtered_fp = os.path.join(config.Config.ar_fp, "result_edge_style_filted")
    save_obj = SaveArMetrics(results=result_edge_style_filted, out_path=res_filtered_fp)
    save_obj.write_result()
    print("finish writing result_edge_style_filtered")

    # export edge style results to file
    output_path = Config.ar_fp
    unfilterd_output_path = os.path.join(output_path, "edge_style", "unfilterd")
    try:
        os.makedirs(unfilterd_output_path)
    except FileExistsError:
        print("FileExists, edgle style output will be overwrited...")
    except Exception as e:
        print(F"Make dir error: {e}")

    cell_type_pipe.export_edge_style(result_edge_style, unfilterd_output_path)

    filterd_output_path = os.path.join(output_path, "edge_style", "filterd")
    try:
        os.makedirs(filterd_output_path)
    except FileExistsError:
        print("FileExists, filtered edgle style output will be overwrited...")
    except Exception as e:
        print(F"Make dir error: {e}")

    cell_type_pipe.export_edge_style(result_edge_style_filted, filterd_output_path)

    # cache result
    cache_result_fp = os.path.join(Config.cahe_result_edge_style_fp)
    save_obj = SaveArMetrics(results=result_edge_style, out_path=cache_result_fp)
    save_obj.write_result()
    print("result_edge_style has been cached...")

    cache_result_filted_fp = os.path.join(Config.cahe_result_edge_style_filted_fp)
    save_filterd_obj = SaveArMetrics(results=result_edge_style_filted, out_path=cache_result_filted_fp)
    save_filterd_obj.write_result()
    print("result_edge_style_filted has been cached...")

# get all result_edge_style_filted genes and run gsea
# gene_dict = cell_type_pipe.get_all_genes(result_edge_style_filted)


# community_detection
graph_meta_dict = {
    "graph_outdir": Config.graph_fp,
    "gsea_outdir": Config.gsea_fp,
    "gsea_meta": Config.gesea_meta,
    "gsea_threshold_dict": Config.gsea_threshold_dict
}

cell_type_pipe.graph_and_gsea_steps(results=result_edge_style_filted,
                                    meta_dict=graph_meta_dict)

# calculate interaction type ar metrics
threshold_dict = {"support": [0.05, None], "confidence": [0.5, None], "lift": [1, None], "leverage": [None, None],
                  "conviction": [1, None]}
ar_metrics_fp = res_filtered_fp
calcu_interaction_type_ar_metrics = CalculateInteractionTypeARMetrics(threshold=threshold_dict,
                                                                      ar_metrics_fp=ar_metrics_fp)
results = calcu_interaction_type_ar_metrics.filter_ar_metrics()
calcu_interaction_type_ar_metrics.run_cell_types_ar_metrics(results)
print(calcu_interaction_type_ar_metrics.result)
output_fp = Config.interaction_type_statis_fp
calcu_interaction_type_ar_metrics.save_result(output_fp)
calcu_interaction_type_ar_metrics.save_result_df(output_fp.replace(".txt", ".tsv"))



print("all task finish...")
