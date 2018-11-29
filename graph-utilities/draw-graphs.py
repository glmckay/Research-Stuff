from graphics import *
import math

"""
    A small tool to draw graphs that are in g6 format

    Requires John Zelle's graphics.py module
"""

HEIGHT = 400
WIDTH = 300
TITLE_CENTER = Point(WIDTH // 2, 20)
TITLE_FONT_SIZE = 14
GRAPH_CENTER = Point(WIDTH // 2, 190)
VERTEX_RAD = 100
LABEL_SPACING = 15
VERTEX_SIZE = 5
VERTEX_COLOUR = 'black'


def g6ToGraph(s):
    numVertices = ord(s[0]) - 63
    adjMatrixSize = (numVertices * (numVertices - 1)) // 2
    adjMatrix = []
    pos = 0
    for c in s[1:]:
        cVal = ord(c) - 63
        for i in range(5, -1, -1):
            adjMatrix.append(cVal >> i & 1)
            pos += 1
            if pos >= adjMatrixSize:
                break
    return (numVertices, adjMatrix)


def areAdjacent(u1,u2, adjMatrix):
    if u1 < u2:
        u = u1
        w = u2
    else:
        u = u2
        w = u1
    return adjMatrix[((w * (w-1)) // 2) + u] == 1


def makeGraphWindow(name="Graph View"):
    return GraphWin(name, WIDTH, HEIGHT, autoflush=False)


def drawGraph(win, n, adj, graphName=""):
    V = []      # vertices
    elems = []  # elements drawn
    for i in range(n):
        directionX = -math.cos(2.0 * math.pi * (i + 0.5) / n)
        directionY = -math.sin(2.0 * math.pi * (i + 0.5) / n)
        pt = Point(GRAPH_CENTER.getX() + VERTEX_RAD * directionX,
                   GRAPH_CENTER.getY() + VERTEX_RAD * directionY)
        textAnchor = Point(GRAPH_CENTER.getX() + (VERTEX_RAD + LABEL_SPACING) * directionX,
                           GRAPH_CENTER.getY() + (VERTEX_RAD + LABEL_SPACING) * directionY)
        V.append(pt)
        vert = Circle(pt, VERTEX_SIZE)
        elems.append(vert)
        vert.setFill(VERTEX_COLOUR)
        vert.draw(win)
        txt = Text(textAnchor, str(i+1))
        elems.append(txt)
        txt.draw(win)

    for w in range(1, n):
        for u in range(w):
            if areAdjacent(u,w,adj):
                l = Line(V[u], V[w])
                l.draw(win)
                elems.append(l)

    title = Text(TITLE_CENTER, graphName)
    title.setSize(TITLE_FONT_SIZE)
    title.draw(win)
    elems.append(title)

    # Update window and wait for mouse click to close it
    win.update()
    win.getMouse()
    for e in elems:
        e.undraw()


# To just draw a graph given the string
def justDrawGraph(g6str):
    win = makeGraphWindow()
    drawGraph(win, *g6ToGraph(g6str), graphName=g6str)
    win.close()


# justDrawGraph("FCOf?")
# justDrawGraph("FCOf_")
# justDrawGraph("FCOfo")
# justDrawGraph("FCOfw")


# with open('graph_data/connected/graphs_4.g6', 'r') as f:
#     win = makeGraphWindow()
#     i = 1
#     for line in f:
#         g6string = line[:-1]
#         drawGraph(win ,*g6ToGraph(g6string), graphName=g6string)
#         i += 1
#         # if (i > 10):
#             # break
#     win.close()