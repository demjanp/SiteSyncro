import multiprocessing as mp
import threading
import time

from typing import Union, Generator, Callable


def _worker(worker_fnc: Callable, params_mp: mp.Queue, collect_mp: mp.Queue, max_queue_size: int, args: list) -> None:
	while True:
		try:
			params = params_mp.get(timeout=10)
		except queue.Empty:
			return
		except:
			time.sleep(0.01)
			continue
		
		while collect_mp.qsize() > max_queue_size:
			time.sleep(0.01)
			pass
		collect_mp.put(worker_fnc(params, *args))


def process_mp(worker_fnc: Callable, params_list: Union[list, Generator], worker_args: list = [],
			   collect_fnc: Callable = None, collect_args: list = [], progress_fnc: Callable = None,
			   progress_args: list = [], max_cpus: int = -1, max_queue_size: int = -1) -> None:
	"""
	Process multiple tasks in parallel using multiprocessing.

	This function takes a worker function and a list of parameters, and applies the worker function to each set of parameters in parallel using multiprocessing. The results are collected and processed by a collector function.

	Parameters:
	worker_fnc (Callable): The worker function to apply to each set of parameters. It should take a set of parameters and any additional arguments, and return a result.
	params_list (list or generator): A list or generator of sets of parameters to apply the worker function to.
	worker_args (list, optional): Additional arguments to pass to the worker function. Default is an empty list.
	collect_fnc (Callable, optional): A function to process the results. It should take a result and any additional arguments. Default is None, which means the results are not processed.
	collect_args (list, optional): Additional arguments to pass to the collector function. Default is an empty list.
	progress_fnc (Callable, optional): A function to report progress. It should take the number of tasks done, the total number of tasks, and any additional arguments. Default is None, which means progress is not reported.
	progress_args (list, optional): Additional arguments to pass to the progress function. Default is an empty list.
	max_cpus (int, optional): The maximum number of CPUs to use for multiprocessing. Default is -1, which means all available CPUs will be used.
	max_queue_size (int, optional): The maximum size of the queue for multiprocessing. Default is -1, which means the queue size is unlimited.

	Returns:
	None
	"""
	
	def call_progress(progress_fnc, done, todo, progress_args):
		
		if progress_fnc == True:
			print("\r%d/%d             " % (done, todo), end="")
		elif callable(progress_fnc):
			progress_fnc(done, todo, *progress_args)
	
	def start_process(proc):
		proc.start()
	
	params_mp = mp.Queue()
	todo = 0
	for params in params_list:
		params_mp.put(params)
		todo += 1
	done = 0
	collect_mp = mp.Queue(todo)
	if max_cpus > 0:
		n_cpus = min(max_cpus, mp.cpu_count() - 1, todo)
	else:
		n_cpus = min(mp.cpu_count() - 1, todo)
	if max_queue_size < 0:
		max_queue_size = n_cpus * 10
	call_progress(progress_fnc, done, todo, progress_args)
	procs = []
	while len(procs) < n_cpus:
		procs.append(
			mp.Process(target=_worker, args=(worker_fnc, params_mp, collect_mp, max_queue_size, worker_args)))
		threading.Thread(target=start_process, args=(procs[-1],)).start()
	while done < todo:
		if not collect_mp.empty():
			data = collect_mp.get()
			done += 1
			call_progress(progress_fnc, done, todo, progress_args)
			if collect_fnc is not None:
				collect_fnc(data, *collect_args)
		else:
			time.sleep(0.5)
	for proc in procs:
		try:
			proc.terminate()
		except:
			pass
		proc = None
