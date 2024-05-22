import networkx as nx
import sys

if sys.platform == "win32":
	from sitesyncro import pygraphviz

def to_agraph(N):
	directed = N.is_directed()
	strict = nx.number_of_selfloops(N) == 0 and not N.is_multigraph()
	A = pygraphviz.AGraph(name=N.name, strict=strict, directed=directed)

	# default graph attributes
	A.graph_attr.update(N.graph.get("graph", {}))
	A.node_attr.update(N.graph.get("node", {}))
	A.edge_attr.update(N.graph.get("edge", {}))

	A.graph_attr.update(
		(k, v) for k, v in N.graph.items() if k not in ("graph", "node", "edge")
	)

	# add nodes
	for n, nodedata in N.nodes(data=True):
		A.add_node(n)
		# Add node data
		a = A.get_node(n)
		a.attr.update({k: str(v) for k, v in nodedata.items()})

	# loop over edges
	if N.is_multigraph():
		for u, v, key, edgedata in N.edges(data=True, keys=True):
			str_edgedata = {k: str(v) for k, v in edgedata.items() if k != "key"}
			A.add_edge(u, v, key=str(key))
			# Add edge data
			a = A.get_edge(u, v)
			a.attr.update(str_edgedata)

	else:
		for u, v, edgedata in N.edges(data=True):
			str_edgedata = {k: str(v) for k, v in edgedata.items()}
			A.add_edge(u, v)
			# Add edge data
			a = A.get_edge(u, v)
			a.attr.update(str_edgedata)

	return A


def pygraphviz_layout(G, prog="dot", root=None, args=""):
	if sys.platform != "win32":
		return {}
	if root is not None:
		args += f"-Groot={root}"
	A = to_agraph(G)
	A.layout(prog=prog, args=args)
	node_pos = {}
	for n in G:
		node = pygraphviz.Node(A, n)
		try:
			xs = node.attr["pos"].split(",")
			node_pos[n] = tuple(float(x) for x in xs)
		except:
			print("no position for node", n)
			node_pos[n] = (0.0, 0.0)
	return node_pos

