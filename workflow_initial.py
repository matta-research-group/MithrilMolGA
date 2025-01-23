import os
import time

def wait_for_file(file_path, sleep_time=120, timeout=10):
    """
    Waits for a file to be produced within a specified timeout period.

    Parameters
    ----------
    file_path (str): The path to the file to wait for.
    sleep_time (int): The time to wait between checking for the file, in seconds. Default is 120 seconds.
    timeout (int): The maximum time to wait for the file, in seconds. Default is 10 seconds.

    Returns
    -------
    bool: True if the file is produced within the timeout period, False otherwise.
    """
    start_time = time.time()
    while not os.path.exists(file_path):
        if time.time() - start_time > timeout:
            return False
        time.sleep(sleep_time)
    # Perform the desired action once the file is present
    return True
