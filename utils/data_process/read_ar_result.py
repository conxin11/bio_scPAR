import numpy as np
import pandas as pd
import os
from utils.data_process.ar_metrics_process import FilterArMetrics

'''
threshold_dict = {"support":[0.05, None], "confidence":[0.5, None], "lift":[1, None], "leverage":[None, None], "conviction":[1, None]}

'''


class ArMetrics:
    def __init__(self, fp=""):
        self.fp = fp
        self.metrics = None

    def load_result(self):
        self.metrics = pd.read_csv(self.fp, index_col=0, header=0)

    def filter_result(self, threshold_dict):
        metrics = self.metrics
        for key, value in threshold_dict.items():
            if value is not None:
                metrics = metrics.loc[(metrics[key] >= value[0]) & (metrics[key] <= value[1])]
        return metrics


