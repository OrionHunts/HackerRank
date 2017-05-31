#!/usr/bin/env python3


class Vertex:
    def __init__(self, vertexName):
        self.name = str(vertexName)
        self._adjList = dict()
        self._adjSet = set()

    def addNeighbour(self, vertexObject, weight=0):
        self._adjList[vertexObject] = weight
        self._adjSet.add(vertexObject)

    @property
    def adjacencyDictionary(self):
        return list(self._adjList.items())

    @property
    def adjacencyList(self):
        return list(self._adjList.keys())

    @property
    def adjacencySet(self):
        return self._adjSet

    @property
    def degree(self):
        return len(self._adjSet)

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
        self._directed = directed
        self._vertexDictionary = dict()
        self._vertexSet = set()

    def addEdge(self, edge):
        if not isinstance(edge, tuple):
            print("An edge must a tuple: ", edge)
            exit(1)
        if len(edge) == 2 and not self._directed:
            uname, vname = str(edge[0]), str(edge[1])
            wt = 0
        elif len(edge) == 3 and self._directed:
            uname, vname, wt = str(t[0]), str(t[1]), float(t[2])
        else:
            print(
                "An edge must a tuple of length 2 (undirected) or 3 (directed): ",
                edge)
            exit(1)
        uobj = self.addVertex(uname)
        vobj = self.addVertex(vname)
        uobj.addNeighbour(vobj, wt)
        if not self._directed:
            vobj.addNeighbour(uobj, wt)

    def addVertex(self, vertexName):
        vertexName = str(vertexName)
        if vertexName in self._vertexSet:
            vertexObject = self._vertexDictionary[vertexName]
        else:
            vertexObject = Vertex(vertexName)
            self._vertexDictionary[vertexName] = vertexObject
            self._vertexSet.add(vertexName)
        return vertexObject

    def bfs(self, srcName):
        visited = set()
        visitSeq = list()
        queue = list()
        queue.append(self._vertexDictionary[srcName])
        while queue:
            vertexObject = queue.pop(0)
            if vertexObject.name not in visited:
                visited.add(vertexObject.name)
                visitSeq.append(vertexObject.name)
                queue.extend(list(vertexObject.adjacencySet - visited))
        return visitSeq, visited

    @property
    def connectedComponents(self):
        components = list()
        visited = set()
        unvisited = (
            vertexName for vertexName in self.vertexSet if vertexName not in visited)
        for vertexName in unvisited:
            visitSeq, component = self.dfs(vertexName)
            visited |= component
            components.append(component)
        return components

    def dfs(self, srcName):
        visited = set()
        visitSeq = list()
        stack = list()
        stack.append(self._vertexDictionary[srcName])
        while stack:
            vertexObject = stack.pop()
            if vertexObject.name not in visited:
                visited.add(vertexObject.name)
                visitSeq.append(vertexObject.name)
                stack.extend(list(vertexObject.adjacencySet - visited))
        return visitSeq, visited

    def fromEdgeList(self, edgeList):
        for edge in edgeList:
            edge = tuple(edge)
            self.addEdge(edge)

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
    def printGraph(self):
        print("The adjacency list representation of the graph is:")
        for vname, vobj in self.vertexDictionary:
            vobj.printAdjacencyList

    @property
    def printGraphWithWeights(self):
        print("The adjacency list representation of the weighted graph is:")
        for vname in self._vertexDictionary:
            self._vertexDictionary[vname].printAdjacencyDictionary

    def printAdjacencyList(self, vertexName):
        self._vertexDictionary[vertexName].printAdjacencyList

    def printAdjacencyDictionary(self, vertexName):
        self._vertexDictionary[vertexName].printAdjacencyDictionary

    @property
    def vertexDictionary(self):
        return list(self._vertexDictionary.items())

    @property
    def vertexList(self):
        return list(self._vertexDictionary.keys())

    @property
    def vertexSet(self):
        return self._vertexSet

def main():
    g = SimpleGraph()
    n, p = map(int, input().split())
    edgeList = list()
    for i in range(p):
        u, v = input().split()
        edgeList.append((u, v))
    g.fromEdgeList(edgeList)
    components = g.connectedComponents
    for i in range(n):
        if str(i) not in g.vertexSet:
            components.append([1])
    total = 0
    s = len(components[0])
    for i in range(1, len(components)):
        total += s * len(components[i])
        s += len(components[i])
    print(total)


if __name__ == '__main__':
    main()
