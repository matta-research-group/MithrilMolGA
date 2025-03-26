import os
import time

def wait_for_file(file_path, sleep_time=180, timeout=21600):
    """
    Waits for a file to be produced within a specified timeout period and checks for errors.

    Parameters
    ----------
    file_path (str): The path to the file to wait for.
    sleep_time (int): The time to wait between checking for the file, in seconds. Default is 180 seconds.
    timeout (int): The maximum time to wait for the file, in seconds. Default is 21600 seconds (6 hours).

    Returns
    -------
    bool: True if the file is produced within the timeout period and no errors are found, False otherwise.
    """
    start_time = time.time()
    err_file_path = file_path.replace('.txt', '.err')
    
    while not os.path.exists(file_path):
        if time.time() - start_time > timeout:
            return False
        if os.path.exists(err_file_path):
            with open(err_file_path, 'r') as err_file:
                if err_file.read().strip():
                    return 'Calculation failed'
        time.sleep(sleep_time)
    
    # Perform the desired action once the file is present
    return True

def is_file_present(file_path, reorganisation=None):
    """
    Checks if a file is present at the given file path. If the file is not present,
    it checks for an error file with the same name but with a '.err' extension.

    Parameters
    ----------
    file_path (str): The path to the file to check.
    
    Returns
    -------
    str: 'Success' if the file is present.
            'Calculation failed' if the error file is present and contains any content.
            'Waiting' if neither the file nor the error file is present.
    """
    if reorganisation == None:
        err_file_path = file_path.replace('_energy_and_gap.txt', '.err')
        py_file_path = file_path.replace('_energy_and_gap.txt', '.py')
        print(err_file_path)

    if reorganisation == 'n_c_geo':
        print('n_c_geo')
        err_file_path = file_path.replace('_n_c_geo_energy_and_gap.txt', '_cation.err')
        py_file_path = file_path.replace('_n_c_geo_energy_and_gap.txt', '_cation.py')
        print(err_file_path)
        print(py_file_path)
    
    if reorganisation == 'opt_c':
        print('opt_c')
        err_file_path = file_path.replace('_opt_c_energy_and_gap.txt', '_cation.err')
        py_file_path = file_path.replace('_opt_c_energy_and_gap.txt', '_cation.py')
        print(err_file_path)
        print(py_file_path)

    if reorganisation == 'opt_a':
        print('opt_a')
        err_file_path = file_path.replace('_opt_a_energy_and_gap.txt', '_anion.err')
        py_file_path = file_path.replace('_opt_a_energy_and_gap.txt', '_anion.py')
        print(err_file_path)
        print(py_file_path)

    if reorganisation == 'n_a_geo':
        print('n_a_geo')
        err_file_path = file_path.replace('_n_a_geo_energy_and_gap.txt', '_anion.err')
        py_file_path = file_path.replace('_n_a_geo_energy_and_gap.txt', '_anion.py')
        print(err_file_path)
        print(py_file_path)
    
    while not os.path.exists(file_path):
        if os.path.exists(err_file_path):
            with open(err_file_path, 'r') as err_file:
                error_content = err_file.read().strip()
                if any(word in error_content for word in ['OptError', 'failed', 'Could not converge SCF', 'ValueError', 'Could not converge geometry optimization']):
                    return 'Calculation failed'
                else:
                    return 'Waiting'
        #if the error file and the txt file don't exist then the calculations has submitted but isn't running yet
        elif os.path.exists(py_file_path):
            return 'Waiting Start'

    return 'Success'


def calculations_status(tasks, sleep_time=10):
    """
    Processes a list of tasks, moving on to the next task and coming back to unfinished tasks later.

    Parameters
    ----------
    tasks (List[Tuple[str, Callable]]): A list of tuples containing task name and callable tasks to be processed.
    sleep_time (int): The time to wait between each loop iteration, in seconds. Default is 10 seconds.

    Returns
    -------
    data_dict (dict): A dictionary with task names as keys and 'Success' or 'Failed' as values.
    failed_dict (dict): A dictionary with task names as keys and the number of failed attempts as values.
    attempts (dict): A dictionary showing how many attempts each calculation took to retreave

    Example Tasks
    -------------
    tasks = [
    ('Task 1', lambda: is_file_present('file1.txt')),
    ('Task 2', lambda: is_file_present('file2.txt')),
    ('Task 3', lambda: is_file_present('file3.txt'))
    ]
    """

    data_dict = {}
    failed_dict = {}
    attempts = {task_name: 0 for task_name, _ in tasks}

    while tasks:
        task_name = (tasks[0][0])
        result = (tasks[0][1])
        tasks.pop(0)  # Remove the task from the list
        print(f'Current directory: {os.getcwd()}')
        #result = task()  # Execute the task
        print(f' {result} {task_name} running')

        if result == 'Success':  # If the task is completed successfully
            print(f'Task {task_name} completed successfully')
            data_dict[task_name] = 'Success'

        elif 'Waiting' in result:
            print(f'Task {task_name} is still waiting')
            attempts[task_name] += 1
            task = lambda: is_file_present(f'{task_name}/{task_name}_opt_energy_and_gap.txt')
            tasks.append((task_name, task()))  # Re-add the task to the end of the list
            print(f'Task {task_name} waiting, will retry (attempt {attempts[task_name]})')

        elif 'Waiting Start' in result:
            print(f'Task {task_name} is awaiting to start')
            attempts[task_name] += 1
            task = lambda: is_file_present(f'{task_name}/{task_name}_opt_energy_and_gap.txt')
            tasks.append((task_name, task()))  # Re-add the task to the end of the list
            print(f'Task {task_name} waiting to start, will retry (attempt {attempts[task_name]})')

        elif 'Calculation failed' in result:
            print(f'Task {task_name} failed')
            failed_dict[task_name] = 'Failed'
        else:
            data_dict[task_name] = 'In Progress'
        
        print(f'sleeping for {sleep_time} seconds, end {task_name}')
        time.sleep(sleep_time)  # Sleep between each loop iteration

    return data_dict, failed_dict, attempts

def reorder_dict(reference_dict, dict_to_reorder):
    reordered_dict = {key: dict_to_reorder[key] for key in reference_dict if key in dict_to_reorder}
    return reordered_dict