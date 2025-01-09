import requests
from tqdm import tqdm
from pathlib import Path

def download_hdf5(filename, output_dir=None):
    filename = filename.replace('.data.h5', '').replace('.head.h5', '')
    year = filename[4:8]
    month = filename[8:10]
    date = filename[10:12]
    time = filename.split('_')[2][:6]

    head_file = filename+'.head.h5'
    data_file = filename+'.data.h5'
    
    # Define base URLs
    primary_base = f"https://vsolar.mssl.ucl.ac.uk/eispac/hdf5/{year}/{month}/{date}/"
    backup_base = f"https://eis.nrl.navy.mil/level1/hdf5/{year}/{month}/{date}/"
    
    # Set default output directory if none provided
    if output_dir is None:
        output_dir = Path('data_eis')
    else:
        output_dir = Path(output_dir)
        
    # Try primary URL first, then backup
    for base_url in [primary_base, backup_base]:
        data_url = base_url + data_file
        head_url = base_url + head_file
        
        # Create directory for downloads if it doesn't exist
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Download files
        data_path = output_dir / data_file
        head_path = output_dir / head_file
        print(f"Attempting download from {base_url}")
        
        data_file_downloaded = download_request(data_url, data_path)
        head_file_downloaded = download_request(head_url, head_path)
        
        # If both files downloaded successfully, we're done
        if data_file_downloaded and head_file_downloaded:
            return True
        
        print(f"Failed to download from {base_url}, trying next source...")
    
    # If we get here, both sources failed
    print("Failed to download from all sources")
    return False

def download_request(url, save_path):
    max_retries = 2
    retry_delay = 10  # seconds
    
    for attempt in range(max_retries):
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            block_size = 1024  # 1 KB
            
            with open(save_path, 'wb') as f, tqdm(
                desc=save_path.name,
                total=total_size,
                unit='iB',
                unit_scale=True,
                unit_divisor=1024,
            ) as progress_bar:
                for data in response.iter_content(block_size):
                    size = f.write(data)
                    progress_bar.update(size)
            
            print(f"Downloaded: {save_path}")
            return True  # Success, return True
        
        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 404:
                print(f"File not found (404 error): {url}")
                return False  # Skip this file
            elif attempt < max_retries - 1:
                print(f"Error occurred: {str(e)}. Retrying in {retry_delay} seconds...")
                import time
                time.sleep(retry_delay)
            else:
                print(f"Failed to download after {max_retries} attempts: {url}")
                print(f"Error: {str(e)}")
                return False  # Skip this file after max retries
        
        except (requests.exceptions.RequestException, OSError) as e:
            if attempt < max_retries - 1:
                print(f"Error occurred: {str(e)}. Retrying in {retry_delay} seconds...")
                import time
                time.sleep(retry_delay)
            else:
                print(f"Failed to download after {max_retries} attempts: {url}")
                print(f"Error: {str(e)}")
                return False  # Skip this file after max retries
