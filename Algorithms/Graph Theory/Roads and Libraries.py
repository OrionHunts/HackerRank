#!/usr/bin/env python3


class Vertex:
    def __init__(self, vertexName):
        self.name = str(vertexName)  # name of the vertex as a string
        # '_adjList' is a dictionary of the form (vertexObject, wt), which denotes that (self, vertexObject) is an edge with weight 'wt'
        self._adjList = dict()
        # '_adjSet' is the set of adjacent vertices of 'self'
        self._adjSet = set()

    # adds a vertexObject and corresponding edge weight to the dictionary
    # of self
    def addNeighbour(self, vertexObject, weight=0):
        self._adjList[vertexObject] = weight
        self._adjSet.add(vertexObject)

    # returns all the vertexObjects, and corresponding weights adjacent to 'self' as a list
    # Note: the list contains tuples of the form (vertexObjects, weight)
    @property
    def adjacencyDictionary(self):
        return list(self._adjList.items())

    # returns all the vertexObjects adjacent to self as a list
    # Note: the list contains vertexObjects
    @property
    def adjacencyList(self):
        return list(self._adjList.keys())

    # returns all the vertexObjects adjacent to self as a set
    # Note: the set contains vertexObjects
    @property
    def adjacencySet(self):
        return self._adjSet

    # returns the degree of the vertex 'self'
    @property
    def degree(self):
        return len(self._adjSet)

    # returns the weight of the edge (self, vertexObject)
    def getWeight(self, vertexObject):
        return self._adjList[vertexObject]

    @property
    def printAdjacencyList(self):
        print("{}: ".format(self.name), end='')
        print([vobj.name for vobj in self.adjacencyList])

    @property
    def printAdjacencyDictionary(self):
        print("{}: ".format(self.name), end='')
        print([(vobj.name, wt) for vobj, wt in self.adjacencyDictionary])


class SimpleGraph:
    def __init__(self, directed=False):
        self._directed = directed  # is the graph directed?
        # a dictionary as (key, value), where key is vertexName and
        # value is the vertexObject
        self._vertexDictionary = dict()
        # the vertex set of the graph, that contains vertexNames (as strings)
        self._vertexSet = set()

    # adds an edge to the graph
    # edge is a tuple (u, v[, w ])
    def addEdge(self, edge):
        if not isinstance(edge, tuple):
            print("An edge must a tuple: ", edge)
            exit(1)
        if len(edge) == 2 and not self._directed:
            uname, vname = str(edge[0]), str(edge[1])
            wt = 0  # default weight for undirected graphs
        elif len(edge) == 3 and self._directed:
            uname, vname, wt = str(t[0]), str(t[1]), float(t[2])
        else:
            print(
                "An edge must a tuple of length 2 (undirected) or 3 (directed): ",
                edge)
            exit(1)
        # get the object references for the two vertices
        uobj = self.addVertex(uname)
        vobj = self.addVertex(vname)
        # add the edge (uobj, vobj, wt)
        uobj.addNeighbour(vobj, wt)
        # if the graph is undirected add edge (vobj, uobj, wt)
        if not self._directed:
            vobj.addNeighbour(uobj, wt)

    # adds vertexName to the graph
    def addVertex(self, vertexName):
        vertexName = str(vertexName)
        # if vertexName is present in the graph, then return the corresponding
        # vertexObject
        if vertexName in self._vertexSet:
            vertexObject = self._vertexDictionary[vertexName]
        # else, create a vertexObject with vertexName as identifier, and return
        # vertexObject
        else:
            vertexObject = Vertex(vertexName)
            self._vertexDictionary[vertexName] = vertexObject
            self._vertexSet.add(vertexName)
        return vertexObject

    # returns the bfs sequence of traversal and the set of visited vertices
    # starting from the given source vertex
    def bfs(self, srcName):
        visited = set()  # contains the name of the visited vertices
        visitSeq = list()  # contains the name of the visited vertices in sequence
        queue = list()
        # add the src vertexObject vertex into the queue
        queue.append(self._vertexDictionary[srcName])
        # continue untill the queue is empty
        while queue:
            vertexObject = queue.pop(0)  # Queue is FIFO
            if vertexObject.name not in visited:
                visited.add(vertexObject.name)
                visitSeq.append(vertexObject.name)
                # add all the neighbours of vertexObject that have not been
                # visited to the queue
                queue.extend(list(vertexObject.adjacencySet - visited))
        return visitSeq, visited

    # returns the dfs sequence of traversal and the set of visited vertices
    # starting from the given source vertex
    def dfs(self, srcName):
        visited = set()  # contains the name of the visited vertices
        visitSeq = list()  # contains the name of the visited vertices in sequence
        stack = list()
        # add the src vertexObject into the stack
        stack.append(self._vertexDictionary[srcName])
        # continue untill the stack is empty
        while stack:
            vertexObject = stack.pop()  # Stack is LIFO
            if vertexObject.name not in visited:
                visited.add(vertexObject.name)
                visitSeq.append(vertexObject.name)
                # add all the neighbours of vertexObject that have not been
                # visited to the stack
                stack.extend(list(vertexObject.adjacencySet - visited))
        return visitSeq, visited

    # construct the graph from the given list of edges
    def fromEdgeList(self, edgeList):
        for edge in edgeList:
            edge = tuple(edge)
            self.addEdge(edge)

    # construct the graph from the given list of edges contained in a file
    # present in the working directory
    def fromFile(self, filename):
        with open(filename, 'r') as f:
            for line in f.readlines():
                line = str(line.strip('\r').strip('\n'))
                if self.directed:
                    u, v, wt = line.split()
                    edge = (u, v, wt)
                else:
                    u, v = line.split()
                    edge = (u, v)
                self.addEdge(edge)

    @property
    def connectedComponents(self):
        components = list()  # this stores the individual components of the graph as a list of sets
        visited = set()  # this keeps track of the visited vertices of the graph
        # a generator that keeps track of unvisited vertices in the graph
        unvisited = (
            vertexName for vertexName in self.vertexSet if vertexName not in visited)
        for vertexName in unvisited:
            visitSeq, component = self.dfs(vertexName)
            visited |= component
            components.append(component)
        # return the list of connected components of the graph
        # Note: number of connected components of the graph = len(components)
        return components  # return sorted(components, key=len) # OPTIONAL

    @property
    def printGraph(self):
        print("The adjacency list representation of the graph is:")
        for vname, vobj in self.vertexDictionary:
            vobj.printAdjacencyList

    @property
    def printGraphWithWeights(self):
        print("The adjacency list representation of the weighted graph is:")
        for vname, vobj in self.vertexDictionary:
            vobj.printAdjacencyDictionary

    def printAdjacencyList(self, vertexName):
        self._vertexDictionary[vertexName].printAdjacencyList

    def printAdjacencyDictionary(self, vertexName):
        self._vertexDictionary[vertexName].printAdjacencyDictionary

    @property
    def vertexDictionary(self):
        return list(self._vertexDictionary.items())

    # returns the list of vertices of the graph
    @property
    def vertexList(self):
        return list(self._vertexDictionary.keys())

    # returns the set of vertices of the graph
    @property
    def vertexSet(self):
        return self._vertexSet


def main():
    t = int(input())
    for i in range(t):
        n, m, clib, croad = map(int, input().split())
        g = SimpleGraph()
        for i in range(m):
            uname, vname = input().split()
            g.addEdge((uname, vname))
        if clib < croad:
            totalCost = n * clib
        else:
            components = g.connectedComponents
            totalComponents = len(components)
            totalComponents += len(set(str(i) for i in range(1, n + 1)) - g.vertexSet)
            totalCost = totalComponents * clib + (n - totalComponents) * croad
        print(totalCost)


if __name__ == '__main__':
    main()
