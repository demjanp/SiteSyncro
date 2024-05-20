from queue import Empty as QueueEmpty
import multiprocessing as mp
from itertools import islice
import threading
import time

from typing import Union, Generator, Callable


def _worker(worker_fnc: Callable, params_mp: mp.Queue, collect_mp: mp.Queue, max_queue_size: int, batch_size: int, args: list) -> None:
	while True:
		try:
			params = params_mp.get(timeout=10)
		except QueueEmpty:
			return
		except:
			time.sleep(0.01)
			continue
		
		while collect_mp.qsize() > max_queue_size:
			time.sleep(0.01)
			pass
		if batch_size < 1:
			collect_mp.put(worker_fnc(params, *args))
		else:
			collect_mp.put([worker_fnc(params_, *args) for params_ in params])


def process_mp(worker_fnc: Callable, params_list: Union[list, Generator], worker_args: list = [],
				collect_fnc: Callable = None, collect_args: list = [], progress_fnc: Callable = None,
				progress_args: list = [], max_cpus: int = -1, max_queue_size: int = -1, batch_size: int = -1) -> None:
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
	batch_size (int, optional): If set to >0, process data in batches. collect_fnc will then be called with a list of results of length <= batch_size

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
	
	def fnc_generator(params_list):
		for params in params_list:
			yield params
	
	def add_to_params(params_mp, params_list, batch_size):
		if batch_size < 1:
			try:
				params = next(params_list)
			except StopIteration:
				return 0
			params_mp.put(params)
			return 1
		to_add = list(islice(params_list, batch_size))
		n_added = len(to_add)
		params_mp.put(to_add)
		return n_added
	
	if isinstance(params_list, list):
		params_list = fnc_generator(params_list)
	
	params_mp = mp.Queue()
	todo = 0
	if max_queue_size > 0:
		if batch_size < 1:
			for params in islice(params_list, max_queue_size):
				params_mp.put(params)
				todo += 1
		else:
			while todo < max_queue_size:
				todo += add_to_params(params_mp, params_list, batch_size)
	else:
		if batch_size < 1:
			for params in params_list:
				params_mp.put(params)
				todo += 1
		else:
			n_added = 1
			while n_added:
				n_added = add_to_params(params_mp, params_list, batch_size)
				todo += n_added
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
			mp.Process(target=_worker, args=(worker_fnc, params_mp, collect_mp, max_queue_size, batch_size, worker_args)))
		threading.Thread(target=start_process, args=(procs[-1],)).start()
	while done < todo:
		to_add = max_queue_size - params_mp.qsize()
		while to_add > 0:
			n_added = add_to_params(params_mp, params_list, batch_size)
			if not n_added:
				break
			todo += n_added
			to_add -= n_added
		if not collect_mp.empty():
			data = collect_mp.get()
			if batch_size < 1:
				done += 1
			else:
				done += len(data)
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
