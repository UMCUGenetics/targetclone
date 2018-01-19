import re

#Obtain a segmentation from the file. We use this to do a p/q arm segmentation. In the main tool, we report LOH events as being on a p/q arm. In the simulations, we divide the positions into p/q arms. 
class Segmentation:
	
	segments = None
	
	def __init__(self, segments = []):
		self.segments = segments
	
	def setSegmentationFromFile(self, segmentationFile):
		#read the segmentation file
		#assign the segmentation to the measurements
		#we store the start, end and segmentation name separately for convenient searching
		#for every segment, make a new object
		
		samples = []
		with open(segmentationFile) as f:
			content = f.readlines()
			#split newlines
			content = [x.strip() for x in content]
			for row in content:
				#split the line further by tab
				columns = re.split(r'\t+', row)
				segment = Segment(columns[0], columns[1], columns[2])
				self.segments.append(segment)
		
	def getSegmentName(self, chromosome, start, end):
		for segment in self.segments:
			if re.match(str(chromosome), segment.name): #if the segment name matches on the current chromosome (so for chr 1 we find 1p and 1q)
				if segment.start <= int(start) and segment.end >= int(end):
					return segment.name
		return None #TODO: properly catch errors here
			
class Segment:
	
	start = 0
	end = 0
	name = ''
	
	def __init__(self, name, start, end):
		self.start = int(start)
		self.end = int(end)
		self.name = name