from utils.data_process.prepare_envs import DirManager
if __name__ == '__main__':
    dir_manager = DirManager()
    dir_manager.clean_cache_dir()
    dir_manager.mk_results_dir()
    dir_manager.clean_results_dir()
    print("done")