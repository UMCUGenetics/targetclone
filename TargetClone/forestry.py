from bokeh.plotting import figure, output_file, show,ColumnDataSource
from bokeh.models.glyphs import Bezier,Segment,AnnularWedge,Annulus
from bokeh.models import Range1d,HoverTool,Circle,LabelSet
from numpy import exp, abs, angle
import numpy as np
from bokeh.models.widgets import Panel, Tabs

from bokeh.palettes import Set2 as colorPalette
colors = colorPalette[8]

from bokeh.layouts import column, row, widgetbox

from bokeh.models.widgets import Paragraph, Div

def plantTree(nameList,parentList,muList,weightList,descList,
	xNodeSpread=.775,lineMaxWidth=8, #.675
	linePullFactor=0.5,lineStartDist=0.1):
	
	if None not in parentList:
		nameList.insert(0,'Precursor')
		parentList.insert(0,None)
		weightList.insert(0,0)
		descList.insert(0,[])
		muList.insert(0,0)

	#maxWeight = max(weightList)
	#weightList = [x/maxWeight for x in weightList]

	dataTable = list(zip(*[nameList,parentList,muList,weightList,descList]))

	nameDict = {}
	parentDict = {}
	for index,node in enumerate(dataTable):
		if node[0] not in nameDict:
			nameDict[node[0]]=[node[1]]
		else:
			print 'ERROR Node specified twice: ' + node[0]
			exit()

		if node[1] not in parentDict:
			parentDict[node[1]]=[]
		parentDict[node[1]].append(node[0])

	byHeightList = []

	def makeTree(curChild,curDepth):
		#byHeightList
		while len(byHeightList) <= curDepth:
			byHeightList.append([])
		byHeightList[curDepth].append(curChild)
		if curChild in parentDict:
			for child in parentDict[curChild]:
				makeTree(child,curDepth+1)

	makeTree(None,0)
	byHeightList.pop(0)

	maxWidth = max([len(x) for x in byHeightList])
	#print 'Maximum amount of nodes on a single line:',maxWidth

	xList = []
	yList = []
	nList = []
	dList = []
	wList = []
	mList = []
	for height,nodes in enumerate(byHeightList):
		heightWidth = len(nodes)
		for horz,node in enumerate(nodes):
			#print height,horz,node,heightWidth
			xList.append((horz-heightWidth/2.)*xNodeSpread)
			yList.append(height)
			nList.append(node)
			dList.append(descList[nameList.index(node)])
			wList.append(weightList[nameList.index(node)])
			mList.append(muList[nameList.index(node)])
	
	yMax = max(yList)
	yList = [yMax-y for y in yList]


	hover = HoverTool(
			tooltips=[
				("Name","@name"),
				("Events", "@desc"),
				("Distance to parent","@weight"),
				("Tumor fraction", "@mu")
			]
		)

	plotSize = max(max(xList),len(byHeightList))

	mins=plotSize/2
	maxes=plotSize/2
	print mins,maxes
	margin = 1
	treePlot = figure(
		plot_width=800, 
		plot_height=800,
		x_range=Range1d(-plotSize*.5-margin, plotSize*.5+margin), 
		y_range=Range1d(-.5-margin, -.5+plotSize+margin),
		x_axis_location=None, 
		y_axis_location=None,
		min_border=64,
		v_symmetry=True,
		h_symmetry=True,
		active_scroll='wheel_zoom',
		tools=['pan','wheel_zoom','tap','save',hover])
	treePlot.xgrid.grid_line_color = None
	treePlot.ygrid.grid_line_color = None


	for parent in parentDict:
		if parent not in nList:
			continue
		parentIndex = nList.index(parent)
		fromX = xList[parentIndex]
		fromY = yList[parentIndex]


		for child in parentDict[parent]:
			childIndex = nList.index(child)
			toX = xList[childIndex]
			toY = yList[childIndex]

			maxD = max(len(x) for x in dList)
			dWidth = len(dList[childIndex])/float(maxD)*lineMaxWidth+1

			oldIndex = nameList.index(child)
			childWeight = weightList[oldIndex]			
		
			anglePull = linePullFactor
			startDist = lineStartDist
			bezGlyph = Bezier(
				x0=fromX,
				y0=fromY - startDist,
				x1=toX,
				y1=toY + startDist,
				cx0=fromX,
				cy0=fromY - anglePull,
				cx1=toX,
				cy1=toY + anglePull, 
				line_color='black', 
				line_width=dWidth,#childWeight*lineMaxWidth,
				line_alpha=1-childWeight*.75)
			treePlot.add_glyph(bezGlyph)

	source = ColumnDataSource(
			data=dict(
				x=xList,
				y=yList,
				desc=dList,
				name=nList,
				weight=wList,
				mu=mList
			)
		)

	renderer = treePlot.circle('x', 'y', 
		source=source,
		color=colors[1],
		line_color='#000000',
		line_width=2,
		radius=0.35) #0.25
	selected_circle = Circle(line_width=2,fill_color=colors[4],line_color='#000000',fill_alpha=1)
	nonselected_circle = Circle(line_width=2,fill_color=colors[2],line_color='#000000',fill_alpha=1)
	renderer.selection_glyph = selected_circle
	renderer.nonselection_glyph = nonselected_circle

	labels = LabelSet(x='x',
		y='y',
		text='name',
		text_align='center',
		text_font_size='10px',
		text_baseline='middle',
		text_alpha=0.9,
		x_offset=0,
		y_offset=0,
		source=source)
	treePlot.add_layout(labels)

	return treePlot


def plantForest(dataList,titles):
	treeList = []
	for dataSet in dataList:
		tree = plantTree(dataSet[0],dataSet[1],dataSet[2],dataSet[3], dataSet[4])
		message = dataSet[5]
		text = Div(text="""""",width=800)
		# We defineren wat text, kunnen we natuurlijk ook gewoon aan de functie doorgeven
		if message != '':
			text = Div(text="""<center><img src="warning_icon.png"></img>The Infinite Sites Assumption was not resolved, reporting the tree with the fewest assumption violations""",width=800)
		treeList.append(column(tree,text))
		
	return treeList


def growForest(dataList,titles):
	treeList = plantForest(dataList,titles)

	tabList = []
	for i,tree in enumerate(treeList):
		treeTab = Panel(child=tree, title=titles[i])
		tabList.append(treeTab)

	tabs = Tabs(tabs=tabList)
	return tabs

