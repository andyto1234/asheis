import pandas as pd
import glob
import os
import multiprocessing
from contextlib import contextmanager
import signal
# from multiprocessing import Pool
from functools import partial
from time import sleep, time  # Import time
from tqdm import tqdm
import numpy as np
import astropy.units as u
from ashmcmc import ashmcmc, interp_emis_temp
import argparse
import platform
from mcmc.mcmc_utils import calc_chi2
from demregpy import dn2dem
import demregpy
import shutil
from mcmc_para import *
import matplotlib
matplotlib.use('Agg')  # Use a non-interactive backend
import sys
from asheis.core import asheis
import gc  # Add garbage collection
import resource  # Add this import for resource usage statistics
import traceback  # Import for detailed error tracking
from datetime import datetime
# Add a global to keep track of process pools
global_pools = []

@contextmanager
def managed_pool(processes, **kwargs):
    """
    Context manager for properly handling process pools.
    Ensures processes are terminated even if an exception occurs.
    
    Parameters
    ----------
    processes : int
        Number of processes to use
    **kwargs
        Additional arguments for Pool
        
    Yields
    ------
    multiprocessing.Pool
        The process pool
    """
    pool = multiprocessing.Pool(processes, **kwargs)
    global_pools.append(pool)
    try:
        yield pool
    finally:
        pool.close()
        pool.join()
        if pool in global_pools:
            global_pools.remove(pool)

def terminate_all_pools():
    """Terminate all active pools"""
    global global_pools
    for pool in global_pools:
        try:
            pool.terminate()
            pool.join()
        except:
            pass
    global_pools = []

def signal_handler(signum, frame):
    """Handle termination signals by cleaning up pools"""
    if multiprocessing.current_process().name == 'MainProcess':
        print(f"\nReceived signal {signum} in MAIN process {os.getpid()}. Cleaning up and exiting...")
    else:
        print(f"\nReceived signal {signum} in WORKER process {os.getpid()} (name: {multiprocessing.current_process().name}). Cleaning up and exiting...")
    
    resources = resource.getrusage(resource.RUSAGE_SELF)
    print(f"Current resource usage: {resources}")
    print(f"Memory usage: {resources.ru_maxrss/1024:.2f} MB, CPU time: {resources.ru_utime:.2f}s")
    
    if multiprocessing.current_process().name == 'MainProcess':
        terminate_all_pools()
    sys.exit(1)# Register signal handlers
signal.signal(signal.SIGINT, signal_handler)   # Ctrl+C
signal.signal(signal.SIGTERM, signal_handler)  # Termination signal

def clean_directory_contents(directory, makedir=True, all=False):
    """
    Remove all contents within SO_EIS_data_{cluster_id} while keeping the directory
    If all=True, clean everything in the directory, otherwise only clean images and dem_columns
    """
    # directory = f'SO_EIS_data_{cluster_id}'
    try:
        if all:
            # Remove all contents while keeping the directory itself
            for item in os.listdir(directory):
                item_path = os.path.join(directory, item)
                if os.path.isdir(item_path):
                    shutil.rmtree(item_path)
                else:
                    os.remove(item_path)
            print(f"Successfully cleaned all contents in {directory}")
        else:
            # Only remove specific directories
            if os.path.exists(os.path.join(directory, 'images')):
                shutil.rmtree(os.path.join(directory, 'images'))
                print(f"Successfully cleaned images/ in {directory}")
            
            if os.path.exists(os.path.join(directory, 'dem_columns')):
                shutil.rmtree(os.path.join(directory, 'dem_columns'))
                print(f"Successfully cleaned dem_columns/ in {directory}")
            
        if makedir==True:
            os.makedirs(os.path.join(directory, 'images'), exist_ok=True)
            os.makedirs(os.path.join(directory, 'dem_columns'), exist_ok=True)
            print(f"Successfully recreated directories in {directory}")
    except FileNotFoundError:
        print(f"Directory not found: {directory}")
    except Exception as e:
        print(f"Error cleaning directory {directory}: {str(e)}")

def check_dens_exists(filename):
    """
    Check if any of the density diagnostic line pairs exist in the EIS file.
    
    Parameters
    ----------
    filename : str
        The EIS file to check
        
    Returns
    -------
    bool
        True if at least one complete density diagnostic line pair exists, False otherwise
    """
    from asheis.eis_density.density_config import DENSITY_DIAGNOSTICS
    import eispac
    
    a = asheis(filename)
    
    # Check each density diagnostic
    for diagnostic in DENSITY_DIAGNOSTICS:
        diag_info = DENSITY_DIAGNOSTICS[diagnostic]
        
        # Check if both lines exist in the file by trying to read their windows
        try:
            nom_template_name = a.dict.get(diag_info['nom_line'], [None])[0]
            denom_template_name = a.dict.get(diag_info['denom_line'], [None])[0]
            
            if nom_template_name and denom_template_name:
                nom_template = eispac.read_template(eispac.data.get_fit_template_filepath(nom_template_name))
                denom_template = eispac.read_template(eispac.data.get_fit_template_filepath(denom_template_name))
                
                # Try to read cubes for both lines
                nom_cube = eispac.read_cube(filename, window=nom_template.central_wave)
                denom_cube = eispac.read_cube(filename, window=denom_template.central_wave)
                # Explicitly delete large objects
                del nom_template, denom_template
                del nom_cube, denom_cube
                gc.collect()  # Force garbage collection

                
                # If we get here, both lines exist
                print(f"DENS FOUND: Found density diagnostic pair: {diag_info['description']}")
                return True
        except Exception as e:
            # If there's an error reading either line, this pair doesn't exist
            continue
    
    
    # No complete pairs found
    print("No complete density diagnostic line pairs found")
    # Explicitly delete asheis instance
    del a
    gc.collect()
    return False

    
def main(args):
    try:
        # Read the dataframe
        df = pd.read_csv('/disk/solar16/st3/python/eis_flare_project_esa_computer/dataframes/s8_filtered_flares_in_fov_AIAEIS_valid_shifts.csv')
        
        # Drop duplicates and reset index
        df = df.drop_duplicates().reset_index(drop=True)    
        
        # Filter dataframe to only include files for this cluster
        cluster_df = df[df.index % args.total_clusters == args.cluster_id].copy()
        
        print(f"Cluster {args.cluster_id}/{args.total_clusters-1} processing {len(cluster_df)} files out of {len(df)} total files")

        # Create cluster-specific completed files tracking
        completed_files_path = f"completed_files_cluster_{args.cluster_id}.txt"
        failed_files_path = f"failed_files_cluster_{args.cluster_id}.txt"
        
        # Create the completed files tracking file if it doesn't exist
        if not os.path.exists(completed_files_path):
            with open(completed_files_path, 'w') as f:
                pass  # Create an empty file
                
        # Create the failed files tracking file if it doesn't exist
        if not os.path.exists(failed_files_path):
            with open(failed_files_path, 'w') as f:
                pass  # Create an empty file

        # Load already completed files
        with open(completed_files_path, 'r') as f:
            completed_files = f.read().splitlines()
            
        # Load already failed files
        with open(failed_files_path, 'r') as f:
            failed_files = f.read().splitlines()

        filenames = cluster_df['python_filename'].tolist()

        for file_num, filename_full in enumerate(filenames):
            # Check if the previous file's results folder is empty and remove it if it is
            if file_num > 0:
                prev_filename = filenames[file_num - 1]
                prev_results_folder = os.path.join('results', os.path.basename(prev_filename).replace('.data.h5', ''))
                if os.path.exists(prev_results_folder) and len(os.listdir(prev_results_folder)) == 0:
                    print(f"Removing empty results folder from previous file: {prev_results_folder}")
                    shutil.rmtree(prev_results_folder)
            
            # Set a 40-minute timeout for each file processing
            file_start_time = time()
            file_timeout = 40 * 60  # 40 minutes in seconds
            
            filedir = f'SO_EIS_data_{args.cluster_id}/'
            filename = filedir+filename_full.replace(" [processing]", '')

            # Skip if already processed or failed
            if filename in completed_files:
                print(f"Skipping already processed file: {filename}")
                continue
                
            if filename in failed_files:
                print(f"Skipping previously failed file: {filename}")
                continue

            try:
                print(f"Processing: {filename}")

                # Check if the results folder exists and contains the combined NPZ file
                results_folder = os.path.join('results', os.path.basename(filename).replace('.data.h5', ''))
                combined_npz_exists = False
                if os.path.exists(results_folder):
                    for file in os.listdir(results_folder):
                        if file.endswith('_combined.npz'):
                            combined_npz_exists = True
                            break

                if combined_npz_exists:
                    print(f"Skipping already processed file (combined NPZ found): {filename}")
                    with open(completed_files_path, "a") as completed_file:
                        completed_file.write(f"{filename}\n")
                    print(f"Added {filename} to completed files list")
                    continue

                file_dir = download_data(filename)
                
                # Check for timeout after download
                if time() - file_start_time > file_timeout:
                    print(f"Timeout of {file_timeout} seconds reached for file {filename}. Moving to next file...")
                    with open(failed_files_path, "a") as failed_file:
                        failed_file.write(f"{filename}\n")
                    print(f"Added {filename} to failed files list (timeout)")
                    clean_directory_contents(f'SO_EIS_data_{args.cluster_id}', all=True)
                    gc.collect()
                    continue
                
                if not check_dens_exists(file_dir):
                    print(f"No density diagnostic pairs found in {filename}, skipping...")
                    # Append the filename to completed files to avoid reprocessing
                    with open(completed_files_path, "a") as completed_file:
                        completed_file.write(f"{filename}\n")
                    print(f"Added {filename} to completed files list")
                    continue

                # Check for timeout before process_data
                if time() - file_start_time > file_timeout:
                    print(f"Timeout of {file_timeout} seconds reached for file {filename}. Moving to next file...")
                    with open(failed_files_path, "a") as failed_file:
                        failed_file.write(f"{filename}\n")
                    print(f"Added {filename} to failed files list (timeout)")
                    clean_directory_contents(f'SO_EIS_data_{args.cluster_id}', all=True)
                    gc.collect()
                    continue
                
                np_file = process_data(filename, 60)
                outdir = os.path.dirname(np_file)
                print(f"Processed: {filename}")

                # Check for timeout before calc_composition
                if time() - file_start_time > file_timeout:
                    print(f"Timeout of {file_timeout} seconds reached for file {filename}. Moving to next file...")
                    with open(failed_files_path, "a") as failed_file:
                        failed_file.write(f"{filename}\n")
                    print(f"Added {filename} to failed files list (timeout)")
                    clean_directory_contents(f'SO_EIS_data_{args.cluster_id}', all=True)
                    gc.collect()
                    continue
                
                line_databases = {
                    "sis": ['si_10_258.37', 's_10_264.23', 'SiX_SX'],
                    "CaAr": ['ca_14_193.87', 'ar_14_194.40', 'CaXIV_ArXIV'],
                    "FeS": ['fe_16_262.98', 's_13_256.69', 'FeXVI_SXIII'],
                }
                calc_composition(filename, np_file, line_databases, args.cores)
                plt.close('all')

                # Append the completed filename to the cluster-specific completed files list *IMMEDIATELY*
                with open(completed_files_path, "a") as completed_file:
                    completed_file.write(f"{filename}\n")
                print(f"Added {filename} to completed files list")

                # Remove both dem_columns and images folders if specified
                if args.unsave:
                    # Remove dem_columns folder
                    dem_column_dir = os.path.join(outdir, 'dem_columns')
                    if os.path.exists(dem_column_dir):
                        shutil.rmtree(dem_column_dir)
                        print(f"Deleted folder: {dem_column_dir}")
                    else:
                        print(f"Folder not found: {dem_column_dir}, not deleted")

                    # Remove images folder
                    images_dir = os.path.join(outdir, 'images')
                    if os.path.exists(images_dir):
                        shutil.rmtree(images_dir)
                        print(f"Deleted folder: {images_dir}")
                    else:
                        print(f"Folder not found: {images_dir}, not deleted")

                    # filedir_path = os.path.dirname(np_file)
                    # if os.path.exists(filedir_path):
                    #     shutil.rmtree(filedir_path)
                    #     print(f"Deleted folder: {filedir_path}")
                    # else:
                    #     print(f"Folder not found: {filedir_path}, not deleted")
                    gc.collect() # Force garbage collection

            except KeyboardInterrupt:
                print("\nProcess interrupted by user. Exiting cleanly...")
                terminate_all_pools()
                sys.exit(1)  # Exit with a non-zero status to indicate interruption

            except Exception as e:
                error_details = traceback.format_exc()
                print(f"An error occurred processing {filename}:")
                print(error_details)
                
                # Record the failed file with error details
                with open(failed_files_path, "a") as failed_file:
                    failed_file.write(f"{filename}\n")
                
                # Create a detailed error log for this specific file
                error_log_dir = "error_logs"
                os.makedirs(error_log_dir, exist_ok=True)
                error_log_file = os.path.join(error_log_dir, f"error_{os.path.basename(filename)}_{int(time())}.log")
                
                with open(error_log_file, "w") as error_file:
                    error_file.write(f"Error processing file: {filename}\n")
                    error_file.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                    error_file.write(f"Error message: {str(e)}\n\n")
                    error_file.write("Traceback:\n")
                    error_file.write(error_details)
                    error_file.write("\n\nSystem info:\n")
                    error_file.write(f"Python version: {sys.version}\n")
                    resources = resource.getrusage(resource.RUSAGE_SELF)
                    error_file.write(f"Memory usage: {resources.ru_maxrss/1024:.2f} MB\n")
                    error_file.write(f"CPU time: {resources.ru_utime:.2f}s\n")
                
                print(f"Error details saved to {error_log_file}")
                
                dir_to_remove = 'results/'+ os.path.basename(filename).replace('.data.h5','')+'/'
                clean_directory_contents(dir_to_remove, makedir=False)
                # Ensure we clean up any running processes
                terminate_all_pools()
                gc.collect() # Force garbage collection

            # Clean up all files within SO_EIS_data directory after processing is complete
            clean_directory_contents(f'SO_EIS_data_{args.cluster_id}', all=True)

            gc.collect() # Force garbage collection
            
            # Final check for timeout at the end of the loop
            if time() - file_start_time > file_timeout:
                print(f"File processing took longer than {file_timeout} seconds. This will be noted for future runs.")
    except Exception as e:
        print(f"Main process error: {e}")
        terminate_all_pools()
        raise
    finally:
        terminate_all_pools()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process EIS data with multiple clusters')
    parser.add_argument('--cluster-id', type=int, required=True, 
                        help='ID of this cluster (e.g., 0, 1, 2, etc.)')
    parser.add_argument('--total-clusters', type=int, required=True,
                        help='Total number of clusters to divide the work')
    parser.add_argument('--cores', type=int, default=1,
                        help='Number of CPU cores to use for processing (default: 1)')
    parser.add_argument('--unsave', action='store_true',
                        help='Delete dem_columns and images folders after processing')
    parser.add_argument('--retry-failed', action='store_true',
                        help='Retry processing files that previously failed')

    args = parser.parse_args()

    # Validate arguments
    if args.cluster_id >= args.total_clusters:
        raise ValueError(f"cluster-id ({args.cluster_id}) must be less than total-clusters ({args.total_clusters})")
    
    if args.cores < 1:
        raise ValueError(f"cores ({args.cores}) must be greater than 0")

    # If retry-failed flag is set, clear the failed files list
    if args.retry_failed:
        failed_files_path = f"failed_files_cluster_{args.cluster_id}.txt"
        if os.path.exists(failed_files_path):
            print(f"Clearing failed files list to retry all previously failed files")
            os.remove(failed_files_path)
            with open(failed_files_path, 'w') as f:
                pass  # Create an empty file

    # Print configuration
    print(f"Configuration:")
    print(f"  Cluster ID: {args.cluster_id}")
    print(f"  Total Clusters: {args.total_clusters}")
    print(f"  CPU Cores: {args.cores}")
    print(f"  Delete intermediate files: {args.unsave}")
    print(f"  Retry failed files: {args.retry_failed}")

    try:
        # Execute main function
        main(args)
    except KeyboardInterrupt:
        print("\nProcess interrupted by user. Terminating all worker processes...")
        terminate_all_pools()
        sys.exit(1)
    except Exception as e:
        print(f"Unhandled exception: {e}")
        terminate_all_pools()
        raise
    finally:
        terminate_all_pools()
        print("All worker processes terminated.")

    # CALL LINE: python eis_mass_composition.py --cluster-id 0 --total-clusters 2 --cores 80 --unsave