import sys
import os
sys.path.append("C:\\Users\\Admin\\Desktop\\matteo\\git_repo\\microstate_based_on_stats\\")

from micros_based_stats.plot_stat_per_k import plot_microstate_Fstatistics
from micros_based_stats.aggregate_peak_backfits import process_microstate_results_Fstat


if __name__ == "__main__":
    base_folder= os.path.join(os.getcwd(), 'results', '1000')
    process_microstate_results_Fstat(base_folder=base_folder)
    plot_microstate_Fstatistics(base_folder=base_folder,alpha=0.05)
